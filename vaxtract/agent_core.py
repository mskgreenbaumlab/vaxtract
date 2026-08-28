"""SDK-free core logic for the extraction agent.

Loads the packet's schema+vocab (the contract), and provides the pure functions
the SDK tools wrap: validate a candidate record, save it (validate-then-write),
re-validate a written file (outer guard / quarantine), and read tables/PDF text.
Importing this module does NOT import claude_agent_sdk, so it stays unit-testable
without the SDK or a model.
"""
from __future__ import annotations

import base64
import importlib.resources
import io
import json
import pathlib
import re
import time

from . import schema  # the contract (schema.py, sibling module)
from . import vocab   # noqa: F401  (controlled vocabularies; re-exported for callers/tests)
from . import table_map  # pure mapping DSL
from .schema import ExtractedPaper

ROOT = pathlib.Path(__file__).resolve().parent

SCHEMA_VERSION: str = schema.SCHEMA_VERSION
# Layer-2 prompt deltas (shipped as package data). importlib.resources keeps this
# working under zipimport/frozen installs, not only loose-file directory installs.
DELTAS_TEXT: str = importlib.resources.files(__package__).joinpath(
    "layer2_prompt_deltas.md"
).read_text(encoding="utf-8")

# File readers are imported lazily (this module stays SDK-free and importable on a
# schema-only `pip install vaxtract`). The reader deps live in optional extras:
# pypdf/openpyxl/python-docx in [agent], Pillow/pymupdf in [figures]. `_require`
# turns a missing reader into an actionable "install this extra" message instead of
# a bare ModuleNotFoundError.
_READER_EXTRA = {
    "pypdf": "agent", "openpyxl": "agent", "docx": "agent",
    "PIL": "figures", "fitz": "figures",
}


def _require(modpath: str):
    """Import a reader dependency, or raise pointing at the extra that ships it."""
    import importlib
    try:
        return importlib.import_module(modpath)
    except ModuleNotFoundError as exc:
        top = modpath.split(".", 1)[0]
        extra = _READER_EXTRA.get(top, "agent")
        raise ModuleNotFoundError(
            f"vaxtract needs {top!r} for this reader, which ships in the "
            f"'{extra}' extra. Install it with:  pip install 'vaxtract[{extra}]'"
        ) from exc


def validate_record(candidate_json: str) -> tuple[bool, str]:
    """True/'VALID' if the candidate satisfies the schema, else False/errors."""
    try:
        ExtractedPaper(**json.loads(candidate_json))
        return True, "VALID"
    except Exception as e:  # pydantic ValidationError or json error
        return False, f"INVALID — fix these and call validate again:\n{e}"


def save_record(candidate_json: str, out_path: str) -> tuple[bool, str]:
    """Validate then write. Writes nothing unless the record is valid."""
    ok, msg = validate_record(candidate_json)
    if not ok:
        return False, f"NOT SAVED — still invalid:\n{msg}"
    paper = ExtractedPaper(**json.loads(candidate_json))
    p = pathlib.Path(out_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(paper.model_dump_json(indent=2, exclude_none=True))
    return True, f"SAVED {out_path} (schema v{SCHEMA_VERSION})"


def _orphaned_epitope_msg(rec: "ExtractedPaper") -> str | None:
    """Health check on a SCHEMA-VALID record: every epitope must link to a parent
    immunizing peptide via parent_peptide_ids. The schema already rejects parent ids
    that don't exist (dangling refs), but it permits an EMPTY list — which silently
    orphans the epitope so the DB join (and the report) can't attach its sequence.
    Returns a fix message if any epitope is orphaned, else None.
    """
    # Schema permits empty parents ("paper did not resolve a long peptide"). The
    # silent-orphan failure mode is an untargeted predicted epitope, not an
    # assayed sequence evidence already points at (Hu Supp11b).
    ev_epi = {e.epitope_paper_id for e in rec.evidence if getattr(e, "epitope_paper_id", None)}
    orphans = [e.paper_local_id for e in rec.epitopes
               if not e.parent_peptide_ids and e.paper_local_id not in ev_epi]
    if not orphans:
        return None
    return (
        f"[gate] epitope linkage broken: {len(orphans)}/{len(rec.epitopes)} epitopes have "
        f"EMPTY parent_peptide_ids (orphaned -> sequences won't join to peptides). finalize "
        f"already tried to DERIVE each parent by sequence containment against the loaded "
        f"immunizing_peptides and could not, so the peptide these epitopes came from is MISSING "
        f"from the record: load it (add_epitope_manifest on the manifest sheet, or add_table over "
        f"the long-peptide column). If you are transcribing epitopes by hand with add_table, map "
        f"parent_peptide_ids with template_list using the SAME id scheme as the immunizing_peptides' "
        f"paper_local_id. NEVER point at a peptide the source doesn't tie the epitope to. "
        f"First orphans: {orphans[:3]}"
    )


def _peptide_count_recon_msg(rec: "ExtractedPaper") -> str | None:
    """STRICT peptide-count reconciliation: the number of immunizing peptides in the
    record must equal the sum of every patient's declared n_peptides_synthesized. A
    silent drop (run 6 dropped 4 indel neoantigens) shows up as declared > present.

    STRICT was chosen over a SOFT (warn-only) variant on 2026-06-03 — see
    docs/superpowers/specs/2026-06-03-peptide-drop-guard-design.md. Revisit if a
    non-personalized / shared-peptide paper legitimately has peptides != sum(synthesized).
    """
    declared = sum(p.n_peptides_synthesized for p in rec.patients)
    actual = len(rec.immunizing_peptides)
    if declared == actual:
        return None
    if declared > actual:
        return (
            f"[gate] peptide count mismatch: patients declare {declared} synthesized "
            f"(sum of n_peptides_synthesized) but the record holds {actual} immunizing_peptides "
            f"(over by {declared - actual}). If the immunizing-peptide lane was loaded from the "
            f"vaccine-peptide table, you COUNTED PREDICTED EPTs as synthesized -- set each "
            f"vaccinated patient's n_peptides_synthesized so the SUM equals {actual} (LONG "
            f"vaccine peptides only). Do not add peptides. Do not count 8-11mer EPTs."
        )
    return (
        f"[gate] peptide count mismatch: patients declare {declared} synthesized "
        f"(sum of n_peptides_synthesized) but the record holds {actual} immunizing_peptides "
        f"(gap {actual - declared} missing). Peptides were dropped -- add the missing "
        f"ones (INCLUDING indel/frameshift neoantigens; the sequence is the short mutant "
        f"neoantigen sequence, not the long mRNA/indel context)."
    )


_PATIENT_POOL_RE = re.compile(r"\bpools?\b", re.I)


def _patients_needing_pool(rec: "ExtractedPaper") -> list[str]:
    """Patients whose evidence describes a POOLED response (quoted_text mentions a
    'pool') but who have no ExtractedPeptidePool entity.

    An in-pool / neoantigen-pool ELISpot response means the peptides were assayed as a
    pool, so the pool should be materialized (members = that patient's in-pool peptides).
    The per-peptide evidence stays stable run-to-run, but the SEPARATE pool entity was
    created ad hoc and got dropped (Rojas pool_P25: present run13, gone run14). This ties
    the two together. The signal is the word "pool" (covers both the table value
    "...in pool" and the main-text "...neoantigen pools"); the check only fires when NO
    pool entity exists for the patient, so explicit-pool papers (Keskin A-D) never trip.
    """
    pooled = {e.patient_paper_id for e in rec.evidence
              if e.quoted_text and _PATIENT_POOL_RE.search(e.quoted_text)}
    have = {p.patient_paper_id for p in rec.pools}
    return sorted(pooled - have)


def _pool_evidence_not_collapsed(rec: "ExtractedPaper") -> list[str]:
    """Patients whose POOLED responses are encoded as per-member rows instead of ONE pool row.

    The canonical rule (resolves the run-to-run evidence-count variance root-caused on Rojas P25):
    a response MEASURED at the pool level is ONE `pool`-target evidence row; an
    `immunizing_peptide`/`epitope` row is reserved for a response the source reports at the
    INDIVIDUAL peptide level (deconvolution). The deterministic signal that a member row is really a
    POOLED measurement is its own quoted_text — Rojas writes "De novo response IN POOL" on each
    member (vs a bare "De novo response" for the deconvoluted ones). So the anti-pattern is, for a
    patient that HAS a pool entity: >=2 non-pool-target evidence rows whose quoted_text says "pool"
    AND no `pool`-target row consolidating them. Those member rows should collapse to one pool row;
    the genuinely-deconvoluted (no "pool" wording) rows stay. Reproduces gold P25 (1 pool + 3 IMP).
    Counts the per-member-pool rows from data alone, so it converges every run to the same shape.
    SOFT nudge, overridable (allow_member_level_pool_evidence) when the source truly deconvolutes all
    members despite pool wording. Distinct from `_patients_needing_pool` (which fires when NO pool
    entity exists at all); this fires when the entity exists but the EVIDENCE wasn't consolidated.
    """
    have_pool_entity = {p.patient_paper_id for p in rec.pools}
    have_pool_evidence = {e.patient_paper_id for e in rec.evidence if e.target_kind == "pool"}
    pooled_member_rows: dict[str, int] = {}
    for e in rec.evidence:
        if (e.target_kind != "pool" and e.quoted_text
                and _PATIENT_POOL_RE.search(e.quoted_text)
                and e.patient_paper_id in have_pool_entity
                and e.patient_paper_id not in have_pool_evidence):
            pooled_member_rows[e.patient_paper_id] = pooled_member_rows.get(e.patient_paper_id, 0) + 1
    return sorted(p for p, n in pooled_member_rows.items() if n >= 2)


def _candidate_bridge_seq_mismatch(rec: "ExtractedPaper") -> list[str]:
    """Candidates (v2.8 P19) whose `selected_peptide_id` resolves to an IMP but
    whose `sequence` DIFFERS from that IMP's sequence.

    The paper-level cross-ref guard only proves the bridge RESOLVES, not that the
    candidate and its IMP are the same peptide. A valid-but-wrong bridge attaches
    candidate C's prioritization scores to the wrong IMP's outcome → label noise in
    the score→outcome training signal (easy to hit when candidates and IMPs are
    loaded in separate add_table passes and ids are matched by hand). SOFT by
    design: the review says flag/needs_review, NOT hard-reject. Returns the sorted
    list of offending candidate paper_local_ids.
    """
    imp_seq = {i.paper_local_id: i.sequence for i in rec.immunizing_peptides}
    bad = []
    for c in rec.candidates:
        if c.selected_peptide_id is not None and c.selected_peptide_id in imp_seq:
            if imp_seq[c.selected_peptide_id] != c.sequence:
                bad.append(c.paper_local_id)
    return sorted(bad)


def _funnel_size_unknown(rec: "ExtractedPaper") -> bool:
    """True when a candidate funnel exists but its true denominator is unrecorded
    (v2.8 P19): candidates non-empty AND n_predicted_reported is None.

    Without the paper-stated predicted count, a TRUNCATED funnel (50 of 322) is
    indistinguishable from a complete one, so any "selection rate" /
    "immunogenic-per-predicted" computed downstream has a silently wrong
    denominator. SOFT honesty nudge, not a hard-reject.
    """
    return bool(rec.candidates) and rec.n_predicted_reported is None


# v2.9 P21 / v2.16 recalibration: the TRIAL-CONSTANT, structural regimen fields of VaccineDelivery.
# The signature must compare only fields that are uniform across an arm by protocol, so divergence
# means agent inconsistency (or a genuine dose-escalation arm -> override). EXCLUDED:
#  - per-patient OUTCOME fields that legitimately vary: n_boost_doses (boosters a patient ACTUALLY
#    received, like n_doses_received / weeks_surgery_to_first_dose — already exempt). Including
#    n_boost_doses made every personalized-vaccine paper falsely diverge (40480654: identical regimen,
#    boosters 4..26 per patient).
#  - FREE-TEXT fields that vary by incidental wording or per-patient dosing: adjuvant_detail,
#    formulation_detail, dose_amount_raw, schedule_detail (39762422: '25 ug'/'38 ug', substudy phrasing).
# KEEPS the enum/numeric structural fields an agent should record identically per arm. Dose-escalation
# (dose_per_peptide_ug varies by cohort) still flags -> genuine, resolved via allow_regimen_divergence.
_REGIMEN_FIELDS = (
    "adjuvant", "dose_basis", "dose_per_peptide_ug", "n_priming_doses",
)


def _regimen_divergence(rec: "ExtractedPaper") -> list[str]:
    """Vaccine platforms (arms) whose patients carry DIVERGENT trial-constant delivery regimen
    (v2.9 P21). The regimen is normally identical for every patient of an arm; divergence usually
    means the agent re-derived it per patient and got it inconsistent. Per-patient delivery fields
    are exempt. Only patients with a vaccine_delivery set are compared; an arm with <2 such
    patients can't diverge. SOFT nudge, overridable for genuine per-arm regimens (dose escalation)."""
    by_arm: dict = {}
    for p in rec.patients:
        if p.vaccine_delivery is not None:
            by_arm.setdefault(p.vaccine_platform, []).append(p.vaccine_delivery)
    bad = []
    for arm, deliveries in by_arm.items():
        if len(deliveries) < 2:
            continue
        sigs = {tuple(getattr(vd, f) for f in _REGIMEN_FIELDS) for vd in deliveries}
        if len(sigs) > 1:
            bad.append(arm)
    return sorted(bad)


def outer_guard(out_path: str) -> tuple[bool, str]:
    """Re-validate a written file; rename to *.QUARANTINED if it does not validate.
    A schema-valid record also passes record-internal health checks (orphaned epitopes;
    strict peptide-count reconciliation)."""
    p = pathlib.Path(out_path)
    if not p.exists():
        return False, f"[guard] file not found: {out_path}"
    try:
        rec = ExtractedPaper(**json.loads(p.read_text()))
    except Exception as e:
        q = p.with_name(p.name + ".QUARANTINED")
        p.rename(q)
        return False, f"[quarantine] invalid record -> {q}: {e}"
    # schema-valid but possibly inconsistent: keep the file (it's valid), tell the agent
    # to fix. Report every failed health check at once so the agent fixes in one pass.
    problems = [m for m in (_peptide_count_recon_msg(rec), _orphaned_epitope_msg(rec)) if m]
    if problems:
        return False, " | ".join(problems)
    return True, f"[gate] {out_path} validates under schema v{SCHEMA_VERSION}"


# Oversized-sheet guards for read_table_rows. Supplementary "source data" sheets can
# carry thousands of rows with multi-kB neoORF-context cells (Keskin Table S5); a naive
# dump blows the CLI's ~64KB per-tool-result cap and spills to a file the agent is
# forbidden to read. We clip each cell and byte-cap the whole preview so the agent
# always gets a readable slice + the TRUE total, and can narrow with row_filter/columns.
_MAX_CELL_CHARS = 240
_MAX_RESULT_BYTES = 48_000


def _clip_cell(v):
    """Truncate an over-long string cell (keeps numbers/None untouched)."""
    if isinstance(v, str) and len(v) > _MAX_CELL_CHARS:
        return v[:_MAX_CELL_CHARS] + f"…(+{len(v) - _MAX_CELL_CHARS} chars)"
    return v


def _trim_width(matrix: list) -> list:
    """Drop trailing all-empty columns. openpyxl often reports a sheet's width as the
    spreadsheet maximum (Keskin S5 → thousands of phantom empty columns), which makes
    even a sparse row serialize past the byte cap and the preview collapse to nothing."""
    width = 0
    for row in matrix:
        for i in range(len(row) - 1, -1, -1):
            if row[i] is not None and str(row[i]).strip() != "":
                width = max(width, i + 1)
                break
    return [row[:width] for row in matrix] if width else matrix


def _header_index(matrix: list) -> int:
    """Index of the real header row, skipping leading TITLE rows (a single non-empty
    cell while the sheet is wider — e.g. 'Table S2. Selective neoantigens…'). These are
    ubiquitous in supplementary xlsx and otherwise make column projection silently fail.
    """
    width = max((len(r) for r in matrix), default=0)
    for i, row in enumerate(matrix):
        nonnull = sum(1 for c in row if c is not None and str(c).strip() != "")
        if nonnull >= 2 or (width <= 1 and nonnull >= 1):
            return i
    return 0


def _marked_cell(cell) -> object:
    """Cell value with contiguous UNDERLINED runs wrapped in <u>…</u> (and merged across
    formatting breaks, e.g. a coloured mutant residue inside the underline). Used to
    surface minimal-epitope substrings that papers mark by underlining within a longer
    peptide (Li Table S2: <u>STYTAYIV</u> inside QLASTYTAYIVGYVHYGDWLK). Non-rich cells
    pass through unchanged so numeric columns keep their type."""
    _rt = _require("openpyxl.cell.rich_text")
    CellRichText, TextBlock = _rt.CellRichText, _rt.TextBlock

    v = cell.value
    if isinstance(v, CellRichText):
        parts, open_u = [], False
        for blk in v:
            if isinstance(blk, TextBlock):
                text = blk.text
                u = bool(getattr(blk.font, "u", None)) if blk.font is not None else False
            else:
                text, u = str(blk), False
            if u and not open_u:
                parts.append("<u>"); open_u = True
            elif not u and open_u:
                parts.append("</u>"); open_u = False
            parts.append(text)
        if open_u:
            parts.append("</u>")
        return "".join(parts)
    if isinstance(v, str) and cell.font is not None and cell.font.underline:
        return f"<u>{v}</u>"
    return v


def _byte_cap(rows: list) -> tuple[list, str]:
    """Drop trailing rows until the JSON fits the per-result byte budget. Returns
    (kept_rows, note) so the agent knows it was capped and how to narrow."""
    if len(json.dumps(rows, default=str)) <= _MAX_RESULT_BYTES:
        return rows, ""
    kept = rows
    while kept and len(json.dumps(kept, default=str)) > _MAX_RESULT_BYTES:
        kept = kept[:-max(1, len(kept) // 10)]
    return kept, " (byte-capped — narrow with row_filter/columns or a smaller max_rows)"


_ZERO_HINT_MAX_SHEETS = 60          # a workbook wider than this is summarised, not enumerated
_ZERO_HINT_DOMAIN_SHOWN = 8
# Wall-clock ceiling for the whole diagnostic. It runs only on an error path and buys back a lane
# worth a 40-minute run, so seconds are cheap -- but a pathological workbook must not hang a tool
# call, so the sweep stops and says how far it got.
_ZERO_HINT_BUDGET_S = 30.0


def _operand_type_note(row_filter: dict, op: str, domain: list) -> str:
    """'' or a note that the operand and the column disagree on TYPE.

    The Hu re-extract's second probe was the same filter with the operand re-typed as a string; it
    got another bare zero. When the column stores numbers and the filter asks for '1' (or the
    reverse), that is the whole story and it should be said outright."""
    if op == "not_empty":
        return ""
    wanted = row_filter.get("equals") if op == "equals" else None
    if wanted is None or not domain:
        return ""
    # The test is not "what type dominates the column" -- a scored column is legitimately MIXED
    # (0 / 1 / 'n.d.'). It is simply: does a cell equal to the operand exist under the OTHER type?
    if isinstance(wanted, str):
        try:
            as_num = float(wanted)
        except ValueError:
            return ""
        if any(isinstance(v, (int, float)) and not isinstance(v, bool) and float(v) == as_num
               for v in domain):
            return (f"the operand {wanted!r} is a STRING but this column stores that value as a "
                    f"NUMBER -- use {as_num:g} (unquoted)")
    elif isinstance(wanted, (int, float)) and not isinstance(wanted, bool):
        if any(isinstance(v, str) and v.strip() == f"{wanted:g}" for v in domain):
            return (f"the operand {wanted!r} is a NUMBER but this column stores that value as "
                    f"TEXT -- use \"{wanted:g}\" (quoted)")
    return ""


def _column_values(rows: list, key, header_hint: int | None, resolve) -> tuple[list, int] | None:
    """(values under `key`, first data row) for one probed sheet, or None when the column is absent.

    `key` is an int index or a header NAME; a name is located by scanning the first rows for it, so
    a sibling sheet whose header sits at a different row is still searched."""
    if isinstance(key, int):
        start = (header_hint + 1) if header_hint is not None else 0
        if not any(key in r for r in rows[start:]):
            return None
        return [r.get(key) for r in rows[start:]], start
    for i, row in enumerate(rows[:12]):
        for col, val in row.items():
            if isinstance(resolve(val), str) and resolve(val).strip() == key:
                return [r.get(col) for r in rows[i + 1:]], i + 1
    return None


def _filter_zero_match_diagnostic(path: str, sheet: str | None, row_filter: dict,
                                  header_row: int | None) -> str:
    """Why a `row_filter` matched nothing, and where it WOULD match.

    Root cause this closes: the Hu 33479501 re-extract probed the memory-timepoint column on the
    WRONG SHEET -- `col_idx:15` against 'Suppl Dataset 4a' -- got a bare `0 rows after filter`,
    retried the operand as a string, got zero again, and abandoned the whole Year 3-4.5 lane. Sheet
    4b matches that same filter on 83 rows. A silent zero is indistinguishable from 'this paper has
    no such data', so it reads as a finding rather than a mistake.

    Uses the census's zip reader, never openpyxl: re-opening Hu's 257MB workbook to explain a failed
    filter would cost ~140s. Sheets past the scan cap are counted, not silently ignored."""
    try:
        from . import tcr_source_census as _census
        import zipfile as _zip
    except Exception:                              # noqa: BLE001
        return ""
    try:
        key = table_map._col_key(row_filter)
        op, err = table_map._sole_operator(row_filter)
        if err or key is None:
            return ""
    except ValueError:
        return ""
    here_domain: list = []
    elsewhere: list = []
    skipped = 0
    unchecked = 0
    deadline = time.monotonic() + _ZERO_HINT_BUDGET_S
    try:
        with _zip.ZipFile(path) as zf:
            members = _census._sheet_index(zf)
            if len(members) > _ZERO_HINT_MAX_SHEETS:
                return ""
            # The sheet the caller was on is examined FIRST, so its value domain is reported even
            # if the sweep of the others runs out of budget.
            members = sorted(members, key=lambda nm: nm[0] != sheet)
            for name, member in members:
                if time.monotonic() > deadline:
                    unchecked += 1
                    continue
                try:
                    if zf.getinfo(member).file_size > _census._MAX_SCAN_BYTES:
                        skipped += 1
                        continue
                except KeyError:
                    continue
                rows = _census._probe_rows_full(zf, member)
                if not rows:
                    continue
                needed = {v[1] for r in rows[:12] for v in r.values() if isinstance(v, tuple)}
                strings = _census._resolve_shared_strings(zf, needed) if needed else {}

                def resolve(v, _s=strings):
                    return _s.get(v[1], "") if isinstance(v, tuple) else v

                found = _column_values(rows, key, header_row if name == sheet else None, resolve)
                if found is None:
                    continue
                values, _ = found
                cell_ids = {v[1] for v in values if isinstance(v, tuple)}
                cells = _census._resolve_shared_strings(zf, cell_ids) if cell_ids else {}
                vals = [cells.get(v[1]) if isinstance(v, tuple) else v for v in values]
                n = sum(1 for v in vals
                        if table_map._matches(v, row_filter, op) and v not in (None, ""))
                if name == sheet:
                    here_domain = vals
                elif n:
                    elsewhere.append((name, n))
    except Exception:                              # noqa: BLE001 - a hint must never break the read
        return ""

    parts = []
    if here_domain:
        counts: dict = {}
        for v in here_domain:
            if v not in (None, ""):
                counts[repr(v)] = counts.get(repr(v), 0) + 1
        top = sorted(counts.items(), key=lambda kv: -kv[1])[:_ZERO_HINT_DOMAIN_SHOWN]
        if top:
            shown = ", ".join(f"{v}x{n}" for v, n in top)
            more = f" (+{len(counts) - len(top)} more distinct)" if len(counts) > len(top) else ""
            parts.append(f"that column on this sheet holds: {shown}{more}")
            type_note = _operand_type_note(row_filter, op, here_domain)
            if type_note:
                parts.append(type_note)
        else:
            parts.append("that column is empty on this sheet")
    else:
        parts.append("that column does not exist on this sheet")
    if elsewhere:
        elsewhere.sort(key=lambda t: -t[1])
        named = "; ".join(f"{n!r} ({c} rows)" for n, c in elsewhere[:6])
        more = f" (+{len(elsewhere) - 6} more)" if len(elsewhere) > 6 else ""
        parts.append(f"the SAME filter matches on other sheet(s) in this file: {named}{more} "
                     f"-- you are probably on the wrong sheet")
    else:
        parts.append("no other sheet in this file matches it either, so this is likely a genuine "
                     "zero rather than a wrong-sheet mistake")
    if skipped or unchecked:
        why = []
        if skipped:
            why.append(f"{skipped} too large to scan")
        if unchecked:
            why.append(f"{unchecked} past the {_ZERO_HINT_BUDGET_S:g}s budget")
        parts.append(f"({' and '.join(why)} -- those sheets were NOT checked)")
    return " | ".join(parts)


def _render_matrix(matrix: list, *, prefix: str, raw_unit: str, filtered_unit: str,
                   max_rows: int = 500, row_filter: dict | None = None,
                   columns: list | None = None, header_row: int | None = None,
                   observe=None, zero_hint=None) -> str:
    """Shared matrix -> capped-preview renderer for read_table (xlsx) and read_docx (docx).

    Keeps the byte-cap / trailing-empty-trim / title-row header skip / row_filter / columns
    behaviour IDENTICAL across both readers. `prefix` is the first output line (e.g.
    'sheets=[...]' or 'tables=3'); `raw_unit`/`filtered_unit` name the source in the two
    output shapes (the plain preview uses raw_unit; the filtered/projected preview uses
    filtered_unit). The xlsx caller passes the same strings it used inline before, so its
    output is unchanged.
    """
    matrix = _trim_width(matrix)
    matrix = [[_clip_cell(c) for c in row] for row in matrix]

    if row_filter is None and not columns:
        head, cap = _byte_cap(matrix[:max_rows])
        return (
            f"{prefix}\nshowing {len(head)}/{len(matrix)} rows of "
            f"{raw_unit}{cap}:\n" + json.dumps(head, default=str)
        )
    if not matrix:
        return f"{prefix}\n{filtered_unit} is empty"
    hi = header_row if header_row is not None else _header_index(matrix)
    headers = [str(h) if h is not None else "" for h in matrix[hi]]
    width = len(headers)
    # Columns and row_filter may reference a column by NAME or by 0-based POSITION
    # (int in `columns`, or col_idx/col_letter in row_filter) — positions reach the
    # blank/duplicated-header columns that names cannot.
    bad = []
    for c in (columns or []):
        if isinstance(c, int):
            if not (0 <= c < width):
                bad.append(f"idx {c} (width {width})")
        elif c not in headers:
            bad.append(c)
    fkey = None
    if row_filter:
        try:
            fkey = table_map._col_key(row_filter)
        except ValueError as e:
            return str(e)
        if isinstance(fkey, int):
            if not (0 <= fkey < width):
                bad.append(f"idx {fkey} (width {width})")
        elif fkey is not None and fkey not in headers:
            bad.append(fkey)
    if bad:
        return f"column(s) {bad} not found; headers: {headers}"
    data = []
    for r in matrix[hi + 1:]:
        if all(c is None for c in r):
            continue
        d = {i: r[i] for i in range(len(r))}                                   # positional keys
        d.update({headers[i]: r[i] for i in range(min(len(headers), len(r)))})  # name keys
        data.append(d)
    total = len(data)
    if row_filter:
        ok, kept = table_map.apply_filter(data, row_filter)
        if not ok:
            return kept  # error message from the DSL
        data = kept
        if observe is not None:
            observe(len(data), total)

    def _label(c):
        if isinstance(c, int):
            h = headers[c] if 0 <= c < width else ""
            return f"{h}[{c}]" if h else f"[{c}]"
        return c

    if columns:
        data = [{_label(c): r.get(c) for c in columns} for r in data]
    else:  # strip positional keys so the plain preview is unchanged
        data = [{k: v for k, v in r.items() if not isinstance(k, int)} for r in data]
    shown, cap = _byte_cap(data[:max_rows])
    note = f" matching {row_filter}" if row_filter else ""
    head = (f"{prefix}\n{filtered_unit}: {len(data)}/{total} data rows"
            f"{note}; showing {len(shown)}{cap}:")
    if row_filter and not data and zero_hint is not None:
        # A bare zero reads as a finding ('this paper has no such data'). Say why, and where the
        # same filter WOULD match -- see `_filter_zero_match_diagnostic`.
        hint = zero_hint()
        if hint:
            head += f"\nZERO MATCHES -- {hint}"
    return head + "\n" + json.dumps(shown, default=str)


# A workbook larger than this is streamed read-only (see read_table_rows). Module-level so tests can
# lower it to exercise the streaming branch on a small file.
_BIG_XLSX_BYTES = 50_000_000


def read_table_rows(path: str, sheet: str | None = None, max_rows: int = 500,
                    row_filter: dict | None = None, columns: list | None = None,
                    underline: bool = False, header_row: int | None = None) -> str:
    """Parse an .xlsx into a capped row preview (Layer 1, highest fidelity).

    Optional INSPECTION args (so you never need a host Grep/Read on a spilled result
    file): `row_filter` ({"col": H, "in"|"equals"|"not_empty": ...}) keeps only matching
    rows and reports the matched/total counts; `columns` (list of header names) projects
    to just those columns. With either arg the result is header-keyed dicts; without
    both it is the raw row matrix (header row first).

    `underline=True` reveals underlined sub-sequences (minimal epitopes) wrapped in
    <u>…</u>. `header_row` forces the 0-based header row (else leading title rows are
    auto-skipped). Cells are clipped and the whole preview is byte-capped so an oversized
    sheet never spills — the TRUE row total is always reported.
    """
    openpyxl = _require("openpyxl")

    # A very large workbook (e.g. Hu 2021 / 33479501's 257MB clonotype supplement) STALLS the default
    # full-load path — the whole book is parsed into memory — so the agent could never read its CDR3
    # sheets (Phase-4 miss: 0 clonotypes despite tabulated sequences). Stream big files read-only
    # (openpyxl read_only + data_only), which opens the same 257MB book in ~9s. read_only mode cannot
    # reveal underline styling, so a big file forces underline OFF with a NOTE — a fair tradeoff, since
    # huge CDR3/clonotype tables carry no epitope-underline markup. Normal files are unchanged.
    try:
        big = pathlib.Path(path).stat().st_size > _BIG_XLSX_BYTES
    except OSError:
        big = False
    show_underline = underline and not big
    wb = openpyxl.load_workbook(path, data_only=True, read_only=big, rich_text=show_underline)
    try:
        ws = wb[sheet] if sheet else wb[wb.sheetnames[0]]
        if show_underline:
            matrix = [[_marked_cell(c) for c in row] for row in ws.iter_rows()]
        else:
            matrix = [list(r) for r in ws.iter_rows(values_only=True)]
        prefix = f"sheets={wb.sheetnames}"
        if underline and big:
            prefix += " | NOTE: underline reveal off (file >50MB, streamed read-only)"
        return _render_matrix(
            matrix, prefix=prefix, raw_unit=ws.title,
            filtered_unit=f"sheet {ws.title!r}", max_rows=max_rows,
            row_filter=row_filter, columns=columns, header_row=header_row,
            observe=lambda matched, total: _observe_filtered_read(
                path, ws.title, matched=matched, total=total, row_filter=row_filter),
            zero_hint=lambda: _filter_zero_match_diagnostic(path, ws.title, row_filter, header_row))
    finally:
        if big:
            wb.close()


def _docx_cell_underline(cell) -> str:
    """Cell text with contiguous UNDERLINED runs wrapped in <u>…</u> (docx analogue of the xlsx
    _marked_cell underline reveal). Paragraph breaks within a cell join with newline."""
    parts, open_u = [], False
    for pi, p in enumerate(cell.paragraphs):
        if pi:
            if open_u:
                parts.append("</u>"); open_u = False
            parts.append("\n")
        for run in p.runs:
            u = bool(run.font.underline)
            if u and not open_u:
                parts.append("<u>"); open_u = True
            elif not u and open_u:
                parts.append("</u>"); open_u = False
            parts.append(run.text)
    if open_u:
        parts.append("</u>")
    return "".join(parts)


def _docx_table_matrix(tbl, underline: bool) -> list:
    """A docx table -> the same row/cell matrix read_table produces (so _render_matrix + add_table
    treat it identically)."""
    rows = []
    for row in tbl.rows:
        rows.append([_docx_cell_underline(c) if underline else c.text for c in row.cells])
    return rows


def _unpivot_list_cells(matrix: list) -> tuple[list, str]:
    """Un-pivot a docx table that packs an ENTIRE COLUMN of values into a single cell as a
    newline-separated list — the merged-cell artifact behind silent row loss.

    Some publishers lay a manifest out so the whole table is ONE header row + ONE data row, where
    each data cell holds every value for that column joined by newlines (e.g. PMID 27274999 Supp
    Table 1: a 2x5 grid whose 5 data cells each carry all 31 peptides). python-docx returns that as a
    2-row table, and `_clip_cell` (240-char cap) then truncates each cell — dropping the tail rows
    with only a "…(+N chars)" marker and no way to recover them (root-caused: the agent saw 23 of 31
    peptides and could not reach the rest). Detect that exact shape and TRANSPOSE it back to one row
    per value.

    STRICT guard so normal tables are never reshaped: exactly two rows (header + one data row); the DATA
    row's every NON-empty cell splits into the SAME count K>=2 on newlines (empty columns allowed and
    back-filled with ""); and the HEADER row is NOT itself such a uniform K>=2 list-row (so it stays a
    genuine header — note the header cells may still wrap across lines, e.g. "Position of\npeptide", which
    is why "header has no newline" is the WRONG test). Anything else -> returned unchanged. Returns
    (matrix, note); note='' when nothing was reshaped."""
    if len(matrix) != 2:
        return matrix, ""
    header, data = matrix[0], matrix[1]

    def _segs(c):
        return [s.strip() for s in str(c).split("\n") if s.strip()] if c is not None else []

    def _uniform_k(row):
        """The single segment-count K>=2 shared by every non-empty cell, else None."""
        counts = {len(_segs(c)) for c in row if _segs(c)}
        return counts.pop() if len(counts) == 1 and next(iter(counts)) >= 2 else None

    k = _uniform_k(data)
    if k is None:                                    # data row isn't a clean K-list -> ordinary table
        return matrix, ""
    if _uniform_k(header) is not None:               # header is ALSO a uniform list -> ambiguous, bail
        return matrix, ""
    cols = [_segs(c) for c in data]
    new = [list(header)] + [[(col[i] if i < len(col) else "") for col in cols] for i in range(k)]
    return new, f"un-pivoted list-cell table (1 data row x {len(header)} cols -> {k} rows)"


def read_docx_from(path: str, table_index: int | None = None, max_rows: int = 500,
                   row_filter: dict | None = None, columns: list | None = None,
                   underline: bool = False, text_offset: int | None = None,
                   max_chars: int = 40_000) -> str:
    """Parse a .docx supplement — its TABLES (Supplementary Table Sn…) and/or its PROSE.

    `.docx` supplements are common for BMC/Molecular-Cancer and Frontiers (Data_Sheet_*.docx) and
    hold real per-entity data the xlsx/pdf readers miss (e.g. 34903219 Table S6 'New neoantigen
    mutations in recurrent tumor'). Modes:
      - table_index omitted -> SUMMARY: table count + each table's rows×cols + its caption (the
        preceding paragraph), since docx tables have no names. Pick one with table_index=N.
      - table_index=N -> that table as a capped preview (0-based), with the SAME row_filter/columns/
        underline/byte-cap behaviour as read_table; bulk-load with add_table just like xlsx.
      - text_offset set -> paragraph PROSE in a paged window (next-offset reported), like read_pdf_text.
    """
    docx = _require("docx")
    Table = _require("docx.table").Table
    Paragraph = _require("docx.text.paragraph").Paragraph

    doc = docx.Document(path)
    # associate each table with its preceding non-empty paragraph (its caption), in document order
    tables, last_para = [], ""
    for block in doc.iter_inner_content():
        if isinstance(block, Paragraph):
            if block.text.strip():
                last_para = block.text.strip()
        elif isinstance(block, Table):
            tables.append((block, last_para)); last_para = ""

    # PROSE mode (paged, mirrors read_pdf_text_from)
    if text_offset is not None:
        txt = "\n".join(p.text for p in doc.paragraphs if p.text.strip())
        total = len(txt); off = max(0, int(text_offset)); chunk = txt[off:off + max_chars]
        end = off + len(chunk)
        more = (f"\n…(+{total - end} chars; call again with text_offset={end})"
                if end < total else "")
        return f"docx prose {off}-{end}/{total} chars:\n{chunk}{more}"

    # SUMMARY mode
    if table_index is None:
        if not tables:
            return (f"docx: 0 tables, {len(doc.paragraphs)} paragraphs. "
                    f"Call with text_offset=0 to read the prose.")
        lines = []
        for i, (t, cap) in enumerate(tables):
            _, note = _unpivot_list_cells(_docx_table_matrix(t, False))
            # report the REAL record count for a list-cell pivot (its physical 2 rows are misleading)
            dims = f"{len(t.rows)}x{len(t.columns)}" + (f" [{note}]" if note else "")
            lines.append(f"  [{i}] {dims}  {cap[:70]!r}")
        summ = "\n".join(lines)
        return (f"docx: {len(tables)} table(s), {len(doc.paragraphs)} paragraphs.\n"
                f"tables (read one with table_index=N):\n{summ}\n"
                f"(prose: call with text_offset=0)")

    # SPECIFIC TABLE mode
    if not (0 <= table_index < len(tables)):
        return f"table_index {table_index} out of range (have {len(tables)} table(s): 0..{len(tables) - 1})"
    tbl, cap = tables[table_index]
    matrix = _docx_table_matrix(tbl, underline)
    matrix, note = _unpivot_list_cells(matrix)   # recover row-per-value from merged/pivoted cells
    unit = f"table[{table_index}] {cap[:70]!r}"
    prefix = f"tables={len(tables)}" + (f"\n{note}" if note else "")
    return _render_matrix(
        matrix, prefix=prefix, raw_unit=unit, filtered_unit=unit,
        max_rows=max_rows, row_filter=row_filter, columns=columns)


def read_pdf_text_from(path: str, max_chars: int = 40_000, offset: int = 0) -> str:
    """Extract selectable text from a PDF, returned in a paging window.

    The default 40k-char window stays under the CLI's per-tool-result cap so the result
    is never spilled to a file (the agent has no host Read tool to re-open a spill). When
    text remains past the window, the result says so and gives the next `offset` to call.
    """
    PdfReader = _require("pypdf").PdfReader

    txt = "\n".join((pg.extract_text() or "") for pg in PdfReader(path).pages)
    total = len(txt)
    offset = max(0, int(offset))
    chunk = txt[offset:offset + max_chars]
    end = offset + len(chunk)
    if end < total:
        nxt = (f"\n\n[truncated: chars {offset}-{end} of {total}. Call read_pdf_text again "
               f"with offset={end} for the next chunk.]")
    else:
        nxt = f"\n\n[end of text: chars {offset}-{total} of {total}.]"
    return chunk + nxt


def survey_sources(path: str, max_chars: int = 14000) -> str:
    """One-call INVENTORY of every supplement under `path` (a paper directory, recursed, or a single
    file) so the agent can LOCATE where data lives before opening anything — the fix for papers that
    bury a manifest among many oddly-named supplements (39972124: the treated-peptide list was lost
    among figure-source-data files, so a per-file hunt gave up at 3 of ~108 peptides).

    Per file, one compact block:
      - .xlsx: each sheet -> name, rows x cols, and its header row (so a 'peptide list' is recognizable
               by its columns even when the sheet/file name is unhelpful)
      - .pdf:  page count + first extractable text (caption/title; flags image-only figure PDFs)
      - .docx: table count + first paragraph (caption) + the first table's header row
    Byte-capped; any files dropped by the cap are listed (coverage is never silently truncated)."""
    openpyxl = _require("openpyxl")
    p = pathlib.Path(path)
    if not p.exists():
        return f"path not found: {path}"
    exts = (".xlsx", ".pdf", ".docx")
    files = sorted([f for f in (p.rglob("*") if p.is_dir() else [p])
                    if f.suffix.lower() in exts and not f.name.startswith("~$")])
    if not files:
        return f"no .xlsx/.pdf/.docx under {path}"

    def _clip(v, n):
        return " ".join(str(v).split())[:n]

    def _hdr(cells):
        cells = list(cells)
        while cells and cells[-1] in (None, ""):
            cells.pop()
        return _clip([("" if c is None else c) for c in cells], 220)

    # REACTIVITY MATRICES, annotated HERE rather than at the finalize gate. Measured 2026-08-25: the
    # gate's block message (which does state these columns) arrives only AFTER every add_table call
    # is already made -- the run authored 9 mappings, then saw the columns, then overrode. The load
    # decision happens at the inventory, so the column map has to be here, exactly like the
    # per-patient sheet-family flag below. Costs one structured census (~6s on Hu's 257MB book).
    reactivity: dict = {}
    try:
        from . import tcr_source_census as _census
        _dir = p if p.is_dir() else p.parent
        for _s in _census.census_tcr_sheets(_dir, with_structure=True).reactivity:
            reactivity[(_s.file.lower(), _s.sheet)] = _s
    except Exception:                       # noqa: BLE001 - the inventory must never fail on this
        reactivity = {}

    out = [f"SOURCE INVENTORY for {path} ({len(files)} file(s)):"]
    used, skipped = len(out[0]), []
    for f in files:
        lines = [f"\n# {f.name}"]
        try:
            ext = f.suffix.lower()
            if ext == ".xlsx":
                wb = openpyxl.load_workbook(str(f), read_only=True)
                fam = {}  # header-signature -> [sheet titles], to spot per-patient sheet families
                for si, ws in enumerate(wb.worksheets):
                    if si >= 80:
                        lines.append(f"  ... (+{len(wb.worksheets) - 80} more sheets)")
                        break
                    head = []
                    for row in ws.iter_rows(min_row=1, max_row=3, values_only=True):
                        if sum(1 for c in row if c not in (None, "")) > sum(1 for c in head if c not in (None, "")):
                            head = list(row)
                    lines.append(f"  sheet '{ws.title}' {ws.max_row or 0}x{ws.max_column or 0} | hdr: {_hdr(head)}")
                    sig = tuple(str(c).strip().lower() for c in head if c not in (None, ""))
                    if sig:
                        fam.setdefault(sig, []).append(ws.title)
                    rx = reactivity.get((f.name.lower(), ws.title))
                    if rx:
                        lines.append(
                            f"  >> REACTIVITY MATRIX: {rx.n_expected} POSITIVE assay call(s) = "
                            f"{rx.n_expected} evidence rows. FIRST CHOICE: call "
                            f"build_reactivity_evidence (it authors the per_cell mapping). "
                            f"Or ONE add_table using per_cell (a plain row mapping reaches only "
                            f"ONE assay column): header_row="
                            f"{max(rx.data_start - 1, 0)}, one cell rule per assay column "
                            f"{list(rx.assay_cols)} with \"equals\": 1, each carrying its own "
                            f"timepoint_phase + assay_stimulation + assay_antigen_format; "
                            f"patient col {rx.patient_col}, immunizing-peptide col {rx.imp_col}"
                            + (f", assay-peptide col {rx.assay_pep_col}"
                               if rx.assay_pep_col is not None else "")
                            + (f", id col {rx.id_col} (recovers a blank patient cell)"
                               if rx.id_col is not None else "")
                            + ". Do NOT filter to `equals: 0` rows -- a 0 is a negative call, not "
                              "an immunogenic one.")
                wb.close()
                # FLAG per-patient sheet FAMILIES (>=5 same-schema tabs) right here at the inventory, where
                # the load decision is made: prescribe ONE bulk add_table over the whole family with the
                # ready-to-copy sheet list. The 33064988 class (~34 'IAP-<patient>' immunogenicity tabs) was
                # chronically under-swept because the agent read tabs one-by-one (single add_table only) and
                # stopped — confirmed live 2026-06-09 via the [add_table] mode log (mode=single, no multi).
                for sig, names in fam.items():
                    if len(names) >= 5:
                        shown = names[:60]
                        more = f" (+{len(names) - 60} more)" if len(names) > 60 else ""
                        lines.append(
                            f"  >> PER-PATIENT SHEET FAMILY: {len(names)} sheets share one schema -> load "
                            f"ALL in ONE add_table call with sheets={shown}{more} (do NOT read them "
                            f"one-by-one); derive the per-sheet key via "
                            "{'col':'__sheet__','extract':'<pattern>'}.")
            elif ext == ".pdf":
                PdfReader = _require("pypdf").PdfReader
                rd = PdfReader(str(f))
                t = ""
                for pg in rd.pages[:2]:
                    t += (pg.extract_text() or "")
                    if t.strip():
                        break
                lines.append(f"  pdf {len(rd.pages)}p | "
                             + (_clip(t, 130) if t.strip() else "(image-only, no extractable text)"))
            elif ext == ".docx":
                _docx = _require("docx")
                d = _docx.Document(str(f))
                cap = next((pp.text.strip() for pp in d.paragraphs if pp.text.strip()), "")
                tinfo = ""
                if d.tables and d.tables[0].rows:
                    t0 = d.tables[0]
                    tinfo = f" | T0 {len(t0.rows)}x{len(t0.columns)} hdr: {_hdr([c.text for c in t0.rows[0].cells])}"
                lines.append(f"  docx {len(d.tables)} table(s) | {_clip(cap, 90)}{tinfo}")
        except Exception as e:
            lines.append(f"  (unreadable: {e})")
        block = "\n".join(lines)
        if used + len(block) > max_chars:
            skipped.append(f.name)
            continue
        out.append(block)
        used += len(block)
    if skipped:
        out.append(f"\n[capped: {len(skipped)} file(s) not shown: {', '.join(skipped[:20])}"
                   + ("…" if len(skipped) > 20 else "")
                   + " — call survey_sources on a subdirectory or single file to see them]")
    return "\n".join(out)


def render_figure_image(pdf_path, page, region=None, *, max_side=1568, dpi=300):
    """Rasterize a figure page (or a fractional crop) to a base64 PNG.

    region = (x0, y0, x1, y1) as fractions of the page (top-left origin); None = whole
    page. The longest side is capped at max_side so a full page stays within the vision
    input budget (crops keep their detail since they are small in area). Returns
    (base64_png_str, width, height).
    """
    # validate before any file I/O so a bad region fails fast (and testably)
    if region is not None:
        x0, y0, x1, y1 = region
        if not (0.0 <= x0 < x1 <= 1.0 and 0.0 <= y0 < y1 <= 1.0):
            raise ValueError(f"region must be ascending fractions in [0,1]; got {region!r}")

    Image = _require("PIL.Image")
    _require("fitz")  # pymupdf — used by render_panel.render_page below
    from .figure_benchmark.render_panel import render_page, crop_frac

    img = render_page(pdf_path, int(page), dpi=dpi)
    if region is not None:
        img = crop_frac(img, region)
    if max(img.width, img.height) > max_side:
        # LANCZOS: sharp, version-stable downscaling — resolution matters for the
        # vision read on the other end.
        scale = max_side / max(img.width, img.height)
        img = img.resize(
            (max(1, int(img.width * scale)), max(1, int(img.height * scale))),
            Image.Resampling.LANCZOS,
        )
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return base64.b64encode(buf.getvalue()).decode("ascii"), img.width, img.height


# --- incremental record assembly (build on disk across turns) ----------------
# section name on ExtractedPaper -> its item model in schema. curator_notes is
# intentionally excluded (human-authored, never extracted).
SECTION_MODEL = {
    "patients": "ExtractedPatient",
    "immunizing_peptides": "ImmunizingPeptide",
    "epitopes": "MinimalEpitope",
    "pools": "ExtractedPeptidePool",
    "evidence": "ExtractedEvidence",
    "survival_outcomes": "SurvivalOutcome",
    # v2.9.1 P21: paper-level cohort delivery/treatment latency (e.g. surgery->vaccine).
    "cohort_latencies": "CohortLatency",
    # v2.8 P19: the prediction->selection->outcome funnel. Wired here so the agent
    # can actually write candidates via add_entities/add_table (the schema model +
    # guards existed but this dispatch did not, so candidates were unreachable).
    "candidates": "NeoantigenCandidate",
    # v2.10 P22: paper-level (cohort-aggregate) response->benefit bridge signals. The
    # PER-PATIENT clinical_benefit_signals ride the patients add-path (nested, no entry
    # here); this entry is only for a cohort-level signal not attributable to one patient.
    "clinical_benefit_signals": "ClinicalBenefitSignal",
    # v2.11 P20: gene-level neoantigen mutations (the funnel's genomic upstream) — wired here so
    # the agent can write them via add_entities/add_table from a mutation-dynamics table.
    "neoantigen_mutations": "NeoantigenMutation",
    # v2.14 #4: the SCREENING BUCKET — per-target manifest rows (bulk "No response" denominator),
    # kept OUT of `evidence`. Wired so a manifest table loads here via add_table/add_entities.
    "screening_readouts": "ScreeningReadout",
    # v2.17 TCR extraction MVP — the four TCR paper-level lists, wired so the agent can write them
    # via add_entities/add_table (from Methods prose, a data-availability statement, or a clean
    # clonotype/flag supp table). tcr_seq_status is a SCALAR paper field (set via meta), not here.
    "data_depositions": "DataDeposition",
    "tcr_seq_methods": "TcrSeqMethod",
    "neoantigen_tcr_flags": "NeoantigenTcrFlag",
    "tcr_clonotypes": "TcrClonotype",
}


def _partial_path(out_path: str) -> str:
    return str(out_path) + ".partial.json"


# --- STICKY manifest-loaded lanes (2026-07-07) --------------------------------------------------------
# A deterministic loader (add_epitope_manifest / add_reactive_tcr) owns the lane it populates: it read the
# WHOLE source table. Live run tcr-wired run2 showed the agent then calling clear_entities on the loaded
# epitope lane and hand-re-adding 2 (99 -> 2 collapse). So a loaded lane is LOCKED: clear_entities refuses
# it. The lock lives in a side file (not the record) so it never reaches schema validation; it is cleared
# on init and removed on finalize success alongside the partial.
def _lock_path(out_path: str) -> str:
    return str(out_path) + ".locked.json"


def _quarantine_path(out_path: str) -> str:
    """Sidecar holding rows finalize removed because no source-supported linkage exists for them
    (see quarantine_orphan_epitopes). Kept OUT of the record so it never reaches schema validation,
    and kept ON DISK so a quarantine is auditable rather than a silent drop."""
    return str(out_path) + ".quarantine.json"


# --- LOAD LEDGER + SOURCE POINTER (2026-08-24, source-aware under-coverage gate) --------------------
# Every gate in this file checks the record against ITSELF; none can see the source, so a record that
# captured a FRACTION of what the supplement offered still finalizes clean (Hu 33479501's keep-file:
# tcr_seq_status='both', 3 methods, ZERO clonotypes, next to a supplement holding 10 loadable track
# sheets / 269 clonotypes). Two sidecars close that: the LEDGER records which source sheets a
# deterministic loader actually read, and the SOURCE pointer records where the source lives, so
# `_tcr_source_coverage_gap` can diff the ledger against `tcr_source_census`. Both live beside the
# record (never inside it) so neither reaches schema validation, exactly like the lock/quarantine
# sidecars above.
# READ observations: (file basename, sheet) -> how many rows a FILTERED read_table matched.
#
# Deliberately IN-PROCESS rather than a tool argument or a sidecar the model can address: the point
# is to record what the agent was SHOWN, and a signal the agent can decline to emit is no signal at
# all (same reasoning as `_source_path`). Only FILTERED reads are recorded — an unfiltered preview's
# row total is the whole sheet (375 rows of Hu 4b, of which 121 are positive), so gating on it would
# false-block. Repeat reads keep the LARGEST match, so a later narrow inspection cannot lower the bar.
_READ_OBSERVATIONS: dict = {}
# Below this many matched rows a filtered read is inspection, not a bulk-load intent -> never gates.
_READ_GATE_MIN_ROWS = 20
# A sheet is under-loaded when the deterministic loaders took less than this share of what was shown.
_READ_GATE_MIN_SHARE = 0.5


def _observe_filtered_read(path: str, sheet: str, *, matched: int, total: int,
                           row_filter: dict | None) -> None:
    key = (pathlib.Path(path).name.lower(), sheet)
    prior = _READ_OBSERVATIONS.get(key)
    if prior and prior["matched"] >= matched:
        return
    _READ_OBSERVATIONS[key] = {"file": pathlib.Path(path).name, "sheet": sheet,
                               "matched": int(matched), "total": int(total),
                               "filter": row_filter}


def reset_read_observations() -> None:
    _READ_OBSERVATIONS.clear()


def _ledger_path(out_path: str) -> str:
    return str(out_path) + ".loaded.json"


def _source_path(out_path: str) -> str:
    """Sidecar naming the paper directory this record is being extracted FROM. Written by the runtime
    (`extraction_agent.extract_paper`) before the agent starts -- NOT by a tool argument, so the
    census cannot be silenced by what the model does or does not pass."""
    return str(out_path) + ".source.json"


_TCR_LOADERS = frozenset({"add_tcr_track", "add_reactive_tcr"})


def read_load_ledger(out_path: str) -> list:
    """[{file, sheet, n_loaded, n_available, loader}] -- one entry per source sheet a deterministic
    TCR loader read into this record. Missing/corrupt sidecar -> []."""
    p = pathlib.Path(_ledger_path(out_path))
    if not p.exists():
        return []
    try:
        rows = json.loads(p.read_text())
        return rows if isinstance(rows, list) else []
    except Exception:
        return []


def _append_load_ledger(out_path: str, *, path: str, sheet: str, loader: str,
                        n_loaded: int, n_available: int) -> None:
    """Record that `loader` read (path, sheet), keeping the file's BASENAME as the key: the census
    walks the paper dir while a tool call may pass any absolute/relative spelling of the same file,
    and a path-shaped mismatch must never read as 'a different sheet' (that would be a false block).
    Re-loading a sheet replaces its entry and bumps `attempts` rather than duplicating it."""
    entries = read_load_ledger(out_path)
    fname = pathlib.Path(path).name
    key = (fname.lower(), sheet)
    prior = next((e for e in entries if (str(e.get("file", "")).lower(), e.get("sheet")) == key), None)
    entry = {"file": fname, "sheet": sheet, "loader": loader,
             "n_loaded": int(n_loaded), "n_available": int(n_available),
             "attempts": int(prior.get("attempts", 1)) + 1 if prior else 1}
    entries = [e for e in entries if (str(e.get("file", "")).lower(), e.get("sheet")) != key]
    entries.append(entry)
    try:
        pathlib.Path(_ledger_path(out_path)).write_text(json.dumps(entries, indent=2))
    except OSError:                     # a sidecar write must never break a successful load
        pass


def set_source_dir(out_path: str, source_dir: str) -> tuple[bool, str]:
    """Point this record at the paper directory it is being extracted from (the source census reads
    it at finalize). Idempotent; call before/around init_record."""
    d = pathlib.Path(source_dir)
    if not d.is_dir():
        return False, f"not a directory: {source_dir}"
    payload = {"source_dir": str(d.resolve())}
    try:
        prior = json.loads(pathlib.Path(_source_path(out_path)).read_text())
        if isinstance(prior, dict) and prior.get("source_dir") == payload["source_dir"]:
            payload = {**prior, **payload}          # keep the pmid stamp across a re-run
    except Exception:
        pass
    pathlib.Path(_source_path(out_path)).parent.mkdir(parents=True, exist_ok=True)
    pathlib.Path(_source_path(out_path)).write_text(json.dumps(payload, indent=2))
    return True, f"source dir recorded: {payload['source_dir']}"


def _read_source_dir(out_path: str, pmid=None) -> str | None:
    """The recorded paper directory, or None when the census must stay SILENT: no sidecar (a record
    finalized outside the extraction runtime -- e.g. `refinalize` -- or a pre-ledger record), a
    directory that no longer exists, or a sidecar stamped with a DIFFERENT pmid than the record being
    finalized (a reused out_path pointing at another paper's supplements, which would otherwise
    manufacture required sheets this paper never had)."""
    try:
        meta = json.loads(pathlib.Path(_source_path(out_path)).read_text())
    except Exception:
        return None
    if not isinstance(meta, dict):
        return None
    d = meta.get("source_dir")
    if not d or not pathlib.Path(d).is_dir():
        return None
    stamped = meta.get("pmid")
    if stamped and pmid and str(stamped) != str(pmid):
        return None
    return str(d)


def _stamp_source_pmid(out_path: str, pmid) -> None:
    """Stamp the record's pmid onto the source sidecar at init, so `_read_source_dir` can refuse a
    sidecar left behind by a different paper at the same out_path."""
    p = pathlib.Path(_source_path(out_path))
    if not pmid or not p.exists():
        return
    try:
        meta = json.loads(p.read_text())
        if isinstance(meta, dict):
            meta["pmid"] = str(pmid)
            p.write_text(json.dumps(meta, indent=2))
    except Exception:
        pass


def _read_locks(out_path: str) -> set:
    p = pathlib.Path(_lock_path(out_path))
    if not p.exists():
        return set()
    try:
        return set(json.loads(p.read_text()))
    except Exception:
        return set()


def _lock_lanes(out_path: str, *lanes) -> None:
    locks = _read_locks(out_path) | {l for l in lanes if l}
    pathlib.Path(_lock_path(out_path)).write_text(json.dumps(sorted(locks)))


def _load_partial(out_path: str):
    """Return (record_or_None, error_or_None). error is set ONLY when the partial
    exists but is corrupt/unparseable; (None, None) means no partial yet."""
    p = pathlib.Path(_partial_path(out_path))
    if not p.exists():
        return None, None
    try:
        return json.loads(p.read_text()), None
    except Exception as e:
        return None, f"partial record at {_partial_path(out_path)} is corrupt: {e}"


def init_partial(out_path: str, paper_meta_json: str) -> tuple[bool, str]:
    """Start a partial record from paper-level fields + empty entity lists."""
    # STICKY (2026-07-07): live run3 showed the agent, blocked from clear_entities on a locked lane,
    # calling init_record AGAIN to wipe the whole record (and the lock) and rebuild by hand. Re-init is
    # NOT a bypass: refuse when an existing partial already has loader-locked lanes. (A stale lock with no
    # partial -- e.g. an interrupted prior run -- is not a live record, so it falls through and is cleared.)
    locked = _read_locks(out_path)
    if locked and pathlib.Path(_partial_path(out_path)).exists():
        return False, (
            f"a deterministic loader already populated {sorted(locked)} in this record; re-initializing "
            f"would discard them. Do NOT re-init to get around a locked lane -- keep the loaded data "
            f"(add_entities appends; clear_entities is only for unlocked lanes). If this paper truly must "
            f"be restarted, delete {_partial_path(out_path)} and its .locked.json first.")
    try:
        meta = json.loads(paper_meta_json)
    except Exception as e:
        return False, f"BAD JSON: {e}"
    if not isinstance(meta, dict):
        return False, "paper meta must be a JSON object"
    dropped = [k for k in meta if k in SECTION_MODEL and meta[k]]
    rec = {k: v for k, v in meta.items() if k not in SECTION_MODEL}
    for s in SECTION_MODEL:
        rec[s] = []
    ok, vmsg = validate_record(json.dumps(rec))
    if not ok:
        return False, f"paper metadata invalid:\n{vmsg}"
    pp = pathlib.Path(_partial_path(out_path))
    pp.write_text(json.dumps(rec))
    pathlib.Path(_lock_path(out_path)).unlink(missing_ok=True)   # fresh record -> no stale lane locks
    pathlib.Path(_quarantine_path(out_path)).unlink(missing_ok=True)   # ...and no stale quarantine sidecar
    pathlib.Path(_ledger_path(out_path)).unlink(missing_ok=True)  # ...and no stale source-load ledger
    pathlib.Path(_alias_path(out_path)).unlink(missing_ok=True)   # ...and no stale id aliases
    reset_read_observations()          # ...and no reads carried over from the previous paper
    _stamp_source_pmid(out_path, rec.get("pmid"))                 # bind the source pointer to THIS paper
    warn = (f" (note: entity-list keys {dropped} in the meta were ignored; "
            "add them with add_entities)") if dropped else ""
    return True, f"INIT ok; partial at {pp}{warn}"


# CONTENT IDENTITY for a LOADER-OWNED lane: what makes two rows the same entity, independent of the
# paper_local_id a given run happened to mint. This is the merge rule `add_epitope_manifest` already
# used inline, promoted to one registry so the freehand append path consumes it too.
#
# Root cause it closes: a deterministic loader's lane is locked against `clear_entities`, but nothing
# stopped `add_entities` from appending the same entities again under the paper's own ids. Hu
# 33479501 phase-6 run1 holds 125 immunizing peptides TWICE -- 125 loader-minted
# `IMP_<seq>_<len>_<sha>` plus 125 hand-transcribed `1-IMP04`-style. Its evidence then pointed at the
# hand-minted ids, which is why an id-keyed comparison against the keep-file reads Jaccard 0.000
# while the biology agrees ~80%. Nothing consumed the identity at WRITE time; now something does.
#
# SCOPE, and why it is exactly this: only inside a lane `add_epitope_manifest` OWNS is a sequence
# guaranteed to identify an entity, because that loader mints ids with `_imp_id` -- a pure function
# of the sequence -- and de-duplicates its own output on the same key, so the lane it produces is
# sequence-unique BY CONSTRUCTION. Any later row repeating one of those sequences is therefore a
# duplicate, provably. The rule must NOT be generalised past that: an unlocked lane may legitimately
# carry several rows per sequence. PMID 33064988 mints PER-PATIENT peptide ids ('IMP-B1-AAA',
# 'IMP-M2-AAA') that pools reference patient by patient, and the Rojas reference record holds 232
# immunizing peptides of which only 223 sequences are distinct. Collapsing those would corrupt real
# data, so lock-ownership -- not sequence-sameness alone -- is the gate.
_ENTITY_IDENTITY = {
    "immunizing_peptides": lambda d: (
        ("seq", d["sequence"].upper()) if d.get("sequence") else None),
    "epitopes": lambda d: (
        ("epi", d["sequence"].upper(), d.get("hla_allele"), d.get("mhc_class"))
        if d.get("sequence") else None),
}


def _alias_path(out_path: str) -> str:
    return str(out_path) + ".aliases.json"


def read_id_aliases(out_path: str) -> dict:
    """{dropped paper_local_id: the surviving one} for duplicates refused at append time. A sidecar,
    like the lock/quarantine/ledger files, so it never reaches schema validation."""
    p = pathlib.Path(_alias_path(out_path))
    if not p.exists():
        return {}
    try:
        d = json.loads(p.read_text())
        return d if isinstance(d, dict) else {}
    except Exception:                              # noqa: BLE001
        return {}


def _record_id_aliases(out_path: str, pairs: list) -> None:
    """Remember that `dropped` meant the same entity as `kept`, so a reference the agent writes to
    the id it THOUGHT it created still resolves. Without this the dedup below would trade a silent
    double-write for a loud dangling reference, which is better but still a failed run."""
    if not pairs:
        return
    aliases = read_id_aliases(out_path)
    for dropped, kept in pairs:
        if dropped and kept and dropped != kept:
            aliases[dropped] = kept
    try:
        pathlib.Path(_alias_path(out_path)).write_text(json.dumps(aliases, indent=2))
    except OSError:                                # a sidecar write must never break a good append
        pass


_ID_REF_FIELDS = ("immunizing_peptide_paper_id", "epitope_paper_id", "candidate_paper_id")
_ID_REF_LIST_FIELDS = ("parent_peptide_ids", "member_peptide_ids")
_ID_REF_LANES = ("evidence", "screening_readouts", "neoantigen_tcr_flags", "tcr_clonotypes",
                 "epitopes", "pools", "candidates")


def apply_id_aliases(rec: dict, aliases: dict) -> int:
    """Rewrite every reference to a de-duplicated entity onto the row that survived. Mutates `rec`;
    returns the number of references rewritten."""
    if not aliases:
        return 0
    n = 0
    for lane in _ID_REF_LANES:
        for row in rec.get(lane) or []:
            if not isinstance(row, dict):
                continue
            for f in _ID_REF_FIELDS:
                v = row.get(f)
                if v in aliases:
                    row[f] = aliases[v]
                    n += 1
            for f in _ID_REF_LIST_FIELDS:
                vals = row.get(f)
                if isinstance(vals, list):
                    new = [aliases.get(v, v) for v in vals]
                    n += sum(1 for a, b in zip(vals, new) if a != b)
                    # a remap can collide two members onto one; keep the lane a set-like list
                    seen, out = set(), []
                    for v in new:
                        if v not in seen:
                            seen.add(v)
                            out.append(v)
                    row[f] = out
    for note in rec.get("curator_notes") or []:
        if isinstance(note, dict) and isinstance(note.get("refs"), list):
            note["refs"] = [aliases.get(r, r) for r in note["refs"]]
    return n


def _dedupe_by_identity(section: str, existing: list, incoming: list,
                        enabled: bool = True) -> tuple[list, list]:
    """(fresh, [(dropped_id, kept_id)]) — incoming rows whose content identity is already present
    are not appended. Also de-duplicates WITHIN the incoming batch. `enabled=False` is a pass-through,
    for a lane no deterministic loader owns (see `_ENTITY_IDENTITY` on why the scope is that tight)."""
    identity = _ENTITY_IDENTITY.get(section) if enabled else None
    if identity is None:
        return incoming, []
    index: dict = {}
    for d in existing:
        k = identity(d)
        if k is not None:
            index.setdefault(k, d.get("paper_local_id"))
    fresh, dropped = [], []
    for nd in incoming:
        k = identity(nd)
        if k is not None and k in index:
            dropped.append((nd.get("paper_local_id"), index[k]))
            continue
        if k is not None:
            index[k] = nd.get("paper_local_id")
        fresh.append(nd)
    return fresh, dropped


def _validate_and_append(out_path: str, section: str, items: list) -> tuple[bool, str]:
    """Validate each item against its sub-model, then append the normalized batch
    (atomic: on any invalid item, append nothing). Shared by append_section and
    table_to_entities.

    Rows whose CONTENT IDENTITY is already in the lane are merged away rather than appended, and the
    id the caller minted is aliased onto the surviving row (see `_ENTITY_IDENTITY`)."""
    if section not in SECTION_MODEL:
        return False, f"unknown section {section!r}; valid: {list(SECTION_MODEL)}"
    rec, err = _load_partial(out_path)
    if err:
        return False, err
    if rec is None:
        return False, "no partial record; call init_record first"
    if not isinstance(items, list):
        return False, "items must be a JSON list"
    model = getattr(schema, SECTION_MODEL[section])
    normalized = []
    for i, it in enumerate(items):
        try:
            normalized.append(json.loads(model(**it).model_dump_json()))
        except Exception as e:
            return False, f"item {i} invalid for {section}:\n{e}\n(batch rejected; partial unchanged)"
    # Only a lane a deterministic loader OWNS can be de-duplicated on content: see `_ENTITY_IDENTITY`.
    fresh, dropped = _dedupe_by_identity(section, rec[section], normalized,
                                         enabled=section in _read_locks(out_path))
    rec[section].extend(fresh)
    pathlib.Path(_partial_path(out_path)).write_text(json.dumps(rec))
    _record_id_aliases(out_path, dropped)
    msg = f"appended {len(fresh)} to {section}; total now {len(rec[section])}"
    if dropped:
        pairs = "; ".join(f"{d!r} -> {k!r}" for d, k in dropped[:6])
        more = f" (+{len(dropped) - 6} more)" if len(dropped) > 6 else ""
        msg += (f" | MERGED AWAY {len(dropped)} row(s) already present in this lane under a "
                f"different paper_local_id: {pairs}{more}. Those ids are now ALIASES for the rows "
                f"that were already there, so references you write to them still resolve -- but "
                f"prefer the surviving id. This is what stops a lane holding the same entity twice.")
    return True, msg


def append_section(out_path: str, section: str, items_json: str) -> tuple[bool, str]:
    """Validate each item against its sub-model, then append the batch."""
    try:
        items = json.loads(items_json)
    except Exception as e:
        return False, f"BAD JSON: {e}"
    return _validate_and_append(out_path, section, items)


def set_safety_summary(out_path: str, safety_json: str) -> tuple[bool, str]:
    """Set the PAPER-LEVEL ``safety_summary`` scalar (validated against the SafetySummary sub-model).

    ``safety_summary`` is a paper-level scalar, NOT an entity list, so it has no ``add_entities``
    path -- this is its only setter (the v2.15 axis shipped without it, so the agent could never
    write it). Re-callable (overwrites). Pass ``null``/``{}`` to clear it back to None.
    """
    try:
        obj = json.loads(safety_json)
    except Exception as e:
        return False, f"BAD JSON: {e}"
    rec, err = _load_partial(out_path)
    if err:
        return False, err
    if rec is None:
        return False, "no partial record; call init_record first"
    if obj in (None, {}):
        rec["safety_summary"] = None
        pathlib.Path(_partial_path(out_path)).write_text(json.dumps(rec))
        return True, "safety_summary cleared (None)"
    try:
        normalized = json.loads(schema.SafetySummary(**obj).model_dump_json())
    except Exception as e:
        return False, f"invalid safety_summary: {e}"
    rec["safety_summary"] = normalized
    pathlib.Path(_partial_path(out_path)).write_text(json.dumps(rec))
    return True, (f"set safety_summary (max_related_grade={normalized.get('max_related_grade')}, "
                  f"any_grade3plus_related={normalized.get('any_grade3plus_related')})")


def clear_section(out_path: str, section: str) -> tuple[bool, str]:
    if section not in SECTION_MODEL:
        return False, f"unknown section {section!r}; valid: {list(SECTION_MODEL)}"
    if section in _read_locks(out_path):
        return False, (
            f"{section!r} was populated by a DETERMINISTIC loader (add_epitope_manifest / add_reactive_tcr "
            f"/ add_tcr_track) "
            f"which read the whole source table -- it is locked and clear_entities is refused (this prevents "
            f"the 99->2 epitope collapse seen live). Keep the loaded set; to ADD a row the table missed use "
            f"add_entities (it appends, dedup runs at finalize). To REPLACE the lane, re-run the loader.")
    rec, err = _load_partial(out_path)
    if err:
        return False, err
    if rec is None:
        return False, "no partial record; call init_record first"
    rec[section] = []
    pathlib.Path(_partial_path(out_path)).write_text(json.dumps(rec))
    return True, f"cleared {section}"


def partial_status(out_path: str) -> tuple[bool, str]:
    rec, err = _load_partial(out_path)
    if err:
        return False, err
    if rec is None:
        return False, "no partial record; call init_record first"
    counts = {s: len(rec.get(s, [])) for s in SECTION_MODEL}
    meta = {k: rec.get(k) for k in ("pmid", "journal", "year", "title", "cohort_size", "indication_summary")}
    msg = f"counts={counts} meta={meta}"
    # TCR-FLAG TARGET INDEX: shown exactly when the record has TCR content but an empty flag lane -- the
    # state the finalize flag-gate blocks on. The flags need REAL target ids (loader-minted ids are opaque
    # and schema `_check_target` rejects a guess), so hand them over here rather than make the agent guess.
    if (rec.get("tcr_clonotypes") or rec.get("tcr_seq_methods")) and not rec.get("neoantigen_tcr_flags"):
        idx = _target_id_index(rec.get("immunizing_peptides") or [], rec.get("epitopes") or [],
                               rec.get("candidates") or [])
        if idx:
            msg += ("\nneoantigen_tcr_flags is EMPTY while TCR content exists -- flag targets you may "
                    "reference (GENE:paper_local_id, copy verbatim): " + idx)
    return True, msg


# Quantitative immunogenicity assays where a magnitude (SFC etc.) is normally reported.
_MAGNITUDE_ASSAYS = {"elispot", "tetramer", "ics", "flow_cytometry"}


def _missing_magnitude_evidence(rec: "ExtractedPaper") -> list:
    """Immunogenic, quantitative-assay evidence rows with NO magnitude at all.

    A row passes when magnitude.value is set OR magnitude.raw is non-empty -- so a raw
    like "not reported" or a per-patient SFC set counts as a deliberate call. Only a
    completely absent magnitude is flagged. SOFT (block-once with override) because
    magnitudes are genuinely optional; see
    docs/superpowers/specs/2026-06-03-magnitude-consistency-design.md.
    """
    miss = []
    for e in rec.evidence:
        if e.outcome == "immunogenic" and e.assay in _MAGNITUDE_ASSAYS:
            m = e.magnitude
            has = m is not None and (m.value is not None or (m.raw and m.raw.strip()))
            if not has:
                miss.append(e.immunizing_peptide_paper_id or e.patient_paper_id or "?")
    return miss


def _evidence_anchor_gap(rec: "ExtractedPaper") -> list[str]:
    """Messages where the recorded evidence materially disagrees with the paper's STATED counts
    (n_immunogenic_reported / n_tested_negative_reported) — the funnel-nudge analogue for evidence.

    REDESIGNED 2026-06-09 after the iris validation batch showed this nudge fired on 5/7 papers,
    almost all FALSE (it routed every known-good paper to needs_review). Two root causes, two fixes:

    (1) IMMUNOGENIC RECALL is now POOLING-AWARE. A pooled response is ONE target_kind='pool' row that
        stands in for its member peptides, but the paper counts those members toward
        n_immunogenic_reported. Counting the pool row as 1 dragged the distinct-target total below the
        reported count (Rojas: 23 distinct + a pool = 30 once members are counted, vs reported 25; HCC
        34903219: 0 distinct -> 39 expanded vs 36). So each immunogenic pool row is expanded by its
        pool's member count before the recall comparison.

    (2) NEGATIVE side is OVER-ONLY (never under). The canonical negative-grain rule records only the
        non-responses the paper NAMES (per patient/pool/peptide) — legitimately FEWER than a cohort-total
        n_tested_negative_reported (agents routinely set it to total-minus-immunogenic, e.g. 39972124's
        83 or Rojas's 8, with 0-few named). UNDER is therefore EXPECTED and must not fire; only a
        material OVER-enumeration (a negative per blank cell — the real grain error) is flagged
        (threshold 1.5x: a 4-vs-3 +1 is noise, a 50-vs-3 blowup is not).

    Both sides are EXEMPT under companion_paper_ref (v2.11.4): a secondary paper defers its whole
    response characterization — immunogenic AND negative — to the companion paper. Empty when no anchor
    is set, the record is a companion paper, or the counts agree."""
    def _tgt(e):
        return (e.immunizing_peptide_paper_id or e.epitope_paper_id
                or e.pool_paper_id or e.candidate_paper_id or "")
    companion = bool(getattr(rec, "companion_paper_ref", None))
    msgs = []
    # --- immunogenic RECALL (pooling-aware): expand each immunogenic pool row by its member count ---
    pool_members = {p.paper_local_id: len(getattr(p, "member_peptide_ids", None) or [])
                    for p in (getattr(rec, "pools", None) or [])}
    seen = set()
    eff = 0
    for e in rec.evidence:
        # pre_existing is a paper-stated positive (31440238 Fig 5 pre-vaccine CD4). The
        # anchor is the count the paper states, not vaccine-induced-only (plan P1).
        if e.outcome not in ("immunogenic", "positive", "pre_existing"):
            continue
        if getattr(e, "target_kind", None) == "pool" and getattr(e, "pool_paper_id", None):
            eff += max(1, pool_members.get(e.pool_paper_id, 1))
        else:
            k = (e.patient_paper_id, _tgt(e), getattr(e, "t_cell_subset", None))
            if k not in seen:
                seen.add(k)
                eff += 1
    if (not companion and rec.n_immunogenic_reported is not None
            and eff < rec.n_immunogenic_reported):
        msgs.append(
            f"recorded {eff} immunogenic responses (pool rows counted by member) but the paper reports "
            f"{rec.n_immunogenic_reported} -- you likely MISSED rows")
    # --- negative GRAIN (over-enumeration only; under is grain-rule-legitimate) ---
    # 'negative' (e.g. tumour non-recognition) is a DIFFERENT finding and is excluded; only
    # outcome=='not_immunogenic' counts against the tested-negative anchor.
    neg = sum(1 for e in rec.evidence if e.outcome == "not_immunogenic")
    if (not companion and rec.n_tested_negative_reported is not None
            and neg > rec.n_tested_negative_reported * 1.5):
        msgs.append(
            f"recorded {neg} not_immunogenic rows but the paper reports "
            f"{rec.n_tested_negative_reported} tested-negative -- you OVER-enumerated; emit a negative "
            f"row only at the grain the paper NAMES (per patient / pool / peptide), never for blank "
            f"table cells")
    return msgs


# --- data-resolution grain (schema v2.16) -------------------------------------
# Ladder finest -> coarsest (mirrors vocab.DATA_RESOLUTIONS). Lower rank = finer grain.
_GRAIN_RANK = {g: i for i, g in enumerate(
    ("per_sequence", "per_mutation", "per_target_gene", "cohort_summary", "clinical_only"))}


def derive_data_resolution(rec: "ExtractedPaper") -> str | None:
    """The ACHIEVED data-resolution grain, derived deterministically from populated content (NOT the
    declared `data_resolution` field). The finest populated layer wins. Used to (a) cross-check the
    agent's declared grain and (b) gate the per-sequence completeness anchors for genuinely coarse
    papers. Defensive (getattr defaults) so it is safe on partial/duck-typed records.
    (per_target_gene is agent-declared-only — not deterministically distinguishable here.)"""
    _g = lambda n: getattr(rec, n, None) or []
    if any(getattr(p, "sequence", None) for p in _g("immunizing_peptides")) \
            or any(getattr(e, "sequence", None) for e in _g("epitopes")):
        return "per_sequence"
    if _g("neoantigen_mutations"):
        return "per_mutation"
    if _g("screening_readouts") or getattr(rec, "n_immunogenic_reported", None) is not None:
        return "cohort_summary"
    if _g("survival_outcomes") or getattr(rec, "safety_summary", None):
        return "clinical_only"
    return None


def _is_faithfully_coarse(rec: "ExtractedPaper") -> bool:
    """True when the paper genuinely reports BELOW per-sequence grain AND no per-sequence manifest was
    available in the source -- so the per-sequence completeness anchors (peptide recall, evidence
    breadth) must NOT fire (n_selected_reported is then a CITED count, like companion_paper_ref). If a
    manifest WAS present (peptide_manifest_present=True), a shortfall is real under-extraction -> the
    anchors still fire. None grain (empty record) is NOT coarse -- it should be flagged, not admitted."""
    # A per-sequence/per-patient layer already exists (peptides or per-patient evidence) -> the anchors
    # are meaningful and must apply; this is not a coarse-only paper.
    if getattr(rec, "immunizing_peptides", None) or getattr(rec, "evidence", None):
        return False
    achieved = derive_data_resolution(rec)
    return (achieved is not None
            and _GRAIN_RANK.get(achieved, 0) > _GRAIN_RANK["per_sequence"]
            and not getattr(rec, "peptide_manifest_present", None))


def _peptide_recall_gap(rec: "ExtractedPaper") -> str | None:
    """Peptide RECALL anchor (prototype): the recorded peptides must meet the paper's STATED
    selected/administered count `n_selected_reported`.

    The strict peptide reconciliation (`_peptide_count_recon_msg`: declared == actual) is GAMEABLE —
    an agent that cannot locate the peptide table can satisfy it by LOWERING every patient's
    n_peptides_synthesized to match the few peptides it scraped (root-caused on 39972124: declared 108
    across 16 patients, then lowered to 3 to match 3 text-scraped peptides → schema-valid but ~97% of
    the peptides missing). Anchoring to the paper's OWN number closes that hole: lowering per-patient
    declarations no longer escapes the check. SOFT/block-once (overridable) and fires only when the
    anchor is set and the record holds materially fewer peptides than the paper selected — so a
    shared-peptide paper (few distinct peptides, many patients) is unaffected as long as its
    n_selected_reported matches its distinct peptide count.

    EXEMPTION (v2.11.4): a SECONDARY-ANALYSIS paper that defers its manifest to a companion paper
    (companion_paper_ref set) reports n_selected_reported as a CITED cohort count, not a target it
    enumerates — the per-sequence list lives in the companion paper (e.g. 39972124's 108 neoantigens
    are in Rojas). So the recall shortfall is expected, not a miss; do not fire."""
    if getattr(rec, "companion_paper_ref", None):
        return None
    # v2.16: a genuinely coarse-grained paper (achieved grain below per_sequence, no manifest available)
    # cites n_selected_reported as a cohort count, not a per-sequence target it enumerates -- admit it.
    if _is_faithfully_coarse(rec):
        return None
    anchor = rec.n_selected_reported
    if anchor is None:
        return None
    actual = len(rec.immunizing_peptides)
    if actual >= anchor:
        return None
    return (
        f"recorded {actual} immunizing_peptides but the paper reports {anchor} selected/administered "
        f"(n_selected_reported) -- you likely MISSED the peptide table; locate it (often a 'vaccine "
        f"peptides' / per-patient supplement) and add_table it. Do NOT lower n_peptides_synthesized "
        f"to pass reconciliation."
    )


def _evidence_breadth_gap(rec: "ExtractedPaper") -> str | None:
    """Evidence COVERAGE/breadth nudge (prototype): flag when recorded responses cover only a small
    fraction of the patients who were actually vaccinated.

    Immunogenicity is usually assessed per vaccinated patient, so evidence covering few of them is
    likely incomplete (root-caused on 33064988: 823 peptides across 48 patients, but evidence for only
    6 — the per-patient immunogenicity lived in ~34 separate supplement sheets and the single run
    opened just a handful before finishing). 'Vaccinated' = patients declaring >=1 synthesized peptide.
    SOFT/block-once (overridable). Threshold is deliberately conservative — fires only at SEVERE
    sparsity (coverage < 1/3 of a >=6-patient vaccinated cohort) — because moderate non-response is real
    biology: across the audited refs the lowest legitimate coverage is Rojas 8/16 (50%) and Keskin 5/8
    (62%), while the under-extraction (33064988) sits at 6/48 (12%). So <1/3 separates skipped-sheets
    from a genuine low response rate without false-positiving on the gold corpus."""
    # v2.16: a genuinely coarse-grained paper reports cohort-level (not per-patient) immunogenicity, so
    # zero per-patient evidence rows is faithful, not incomplete -- don't apply the per-sequence breadth
    # anchor. (A manifest being present => still apply it.)
    if _is_faithfully_coarse(rec):
        return None
    # v2.13 A1: a methodological arm (model_antigen_validation / healthy_donor) is NOT a vaccinated
    # disease cohort — exclude it from the denominator so it doesn't dilute coverage. None/patient/
    # tumor_model count (None = legacy untagged).
    _real = ("patient", "tumor_model", None)
    vaccinated = {p.paper_local_id for p in rec.patients
                  if (p.n_peptides_synthesized or 0) > 0
                  and getattr(p, "cohort_kind", None) in _real}
    if len(vaccinated) < 6:
        return None
    covered = {e.patient_paper_id for e in rec.evidence if e.patient_paper_id}
    n_cov = len(vaccinated & covered)
    if n_cov * 3 >= len(vaccinated):   # coverage >= 1/3 -> plausible response rate, don't nudge
        return None
    return (
        f"immune evidence covers {n_cov} of {len(vaccinated)} vaccinated patients -- immunogenicity is "
        f"usually assessed per patient, so this looks INCOMPLETE. If responses live in per-patient "
        f"figures/sheets (e.g. one immunogenicity sheet per patient), iterate ALL of them; do not stop "
        f"after a few. If the trial truly immune-monitored only this subset, proceed with the override."
    )


# class-II coverage anchor (v2.12 P23; cue-split v2.19): fire only when a CLASS-II RESTRICTION is
# quoted (named DR/DP/DQ / I-A/I-E / 'class II-restricted') and ZERO non-inferred class-II records
# exist. A CD4/helper subset cue is NOT a restriction — Keskin-style CD4 ELISpots against a long
# peptide stay t_cell_subset='cd4' with mhc_class unset/not_determined. SOFT/block-once, overridable
# with allow_missing_class_ii (HARD -> needs_review) for the restriction-cue case only.
_CII_RESTRICTION_CUE = re.compile(
    r"class[ \-]?ii|hla-?d|\bdr\b|\bdp\b|\bdq\b|drb|dpa|dpb|dqa|dqb|h-2i|\bi-[ae][a-z]", re.I
)
_NAMED_CII_ALLELE = re.compile(r"hla-?d[rpq]|drb|dpa|dpb|dqa|dqb|h-2i|\bi-[ae][a-z]", re.I)


def _class_ii_minting_gap(rec: "ExtractedPaper") -> list[str]:
    # v2.13 (Fable review): a GUESS must not silence the anti-guess guard — count only NON-inferred
    # class-II epitopes (mhc_class_inferred=True is a heuristic length-call, not a reported class).
    n_class_ii_epi = sum(1 for e in rec.epitopes
                         if getattr(e, "mhc_class", None) == "II"
                         and not getattr(e, "mhc_class_inferred", False))
    n_class_ii_ev = sum(1 for e in rec.evidence if getattr(e, "mhc_class", None) == "class_ii")
    if n_class_ii_epi or n_class_ii_ev:
        return []  # something class-II IS reported (non-inferred) -> not the under-typing failure mode
    texts = [t for t in (
        *[e.quoted_text for e in rec.evidence],
        *[e.quoted_text for e in rec.epitopes],
    ) if t]
    has_restriction_cue = any(_CII_RESTRICTION_CUE.search(t) for t in texts) or any(
        _NAMED_CII_ALLELE.search(e.hla_allele or "") for e in rec.epitopes if getattr(e, "hla_allele", None)
    )
    has_named = any(_NAMED_CII_ALLELE.search(t) for t in texts) or any(
        _NAMED_CII_ALLELE.search(e.hla_allele or "") for e in rec.epitopes if getattr(e, "hla_allele", None)
    )
    msgs = []
    if has_restriction_cue:
        msgs.append(
            "a class-II restriction cue is quoted but ZERO class-II records were minted "
            "(no MinimalEpitope mhc_class='II', no evidence mhc_class='class_ii') -- type "
            "evidence.mhc_class='class_ii' and/or mint a class-II epitope for a defined sequence "
            "with that quoted restriction (CD4/helper alone is not a class-II cue)")
    if has_named:
        msgs.append(
            "a named DR/DP/DQ (or I-A/I-E) restriction is quoted but ZERO class-II epitopes were minted "
            "-- add the class-II MinimalEpitope for the quoted restriction")
    return msgs


# CD4 evidence grain (v2.19 follow-up): a CD4 measurement is the LONG immunizing peptide, not the
# nested class-I 8-11mer. Keskin 30568305 (archived: outputs/archive/2026-08-22_classii_keskin/classii_keskin_2026-08-22) typed subset/class honestly
# then hung CD4 ELISpots on invented EPI_RESP_* aliases of Table-S5 9-mers (HAKLETALR etc.). Gold
# already points those rows at the IMPs. Block-once; HARD override allow_cd4_on_class_i_minimal
# only if the paper quotes that exact 9-mer as the CD4 target (needs_review).
_CLASS_I_MINIMAL_MAX_AA = 11


def _cd4_class_i_minimal_target_gap(rec: "ExtractedPaper") -> list[str]:
    epis = {getattr(e, "paper_local_id", None): e for e in (rec.epitopes or [])}
    msgs: list[str] = []
    for ev in rec.evidence or []:
        if getattr(ev, "t_cell_subset", None) != "cd4":
            continue
        if getattr(ev, "target_kind", None) != "epitope":
            continue
        eid = getattr(ev, "epitope_paper_id", None)
        epi = epis.get(eid)
        seq = (getattr(epi, "sequence", None) or "") if epi is not None else ""
        mhc = getattr(epi, "mhc_class", None) if epi is not None else None
        if epi is None or mhc != "I" or not (1 <= len(seq) <= _CLASS_I_MINIMAL_MAX_AA):
            continue
        parents = [p for p in (getattr(epi, "parent_peptide_ids", None) or []) if p]
        parent_hint = (
            "parent immunizing_peptide_paper_id(s) already on that epitope: " + ", ".join(parents)
            if parents else
            "choose the long immunizing peptide whose sequence contains this fragment"
        )
        msgs.append(
            f"CD4 evidence targets epitope {eid!r} (sequence {seq}, mhc_class=I, len={len(seq)}) "
            f"-- a CD4 / long-peptide measurement is immunizing_peptide grain, not the nested "
            f"class-I 8-11mer. Set target_kind='immunizing_peptide' and immunizing_peptide_paper_id "
            f"to the vaccine peptide ({parent_hint}). Do not mint an EPI_RESP alias of this 9-mer."
        )
    return msgs


# CD8 evidence grain: block CD8-on-long-IMP only when the row's OWN quote NAMES a nested
# class-I 8-11mer (sequence token or Keskin-style EPT12A). Bare "epitopes (EPT)" is not a
# named bite -- SLX4MUT with four Table-S5 EPTs and no which-one stays on the IMP (decision A).
# HARD override allow_cd8_on_immunizing_peptide = paper named an EPT and we still used the IMP.
_NAMED_EPT_ID = re.compile(r"\bEPT\s*\d+[A-Z]\b", re.I)


def _nested_class_i_minimals(rec: "ExtractedPaper", imp_id: str | None, imp_seq: str) -> list:
    hits = []
    for e in rec.epitopes or []:
        if getattr(e, "mhc_class", None) != "I":
            continue
        seq = getattr(e, "sequence", None) or ""
        if not (1 <= len(seq) <= _CLASS_I_MINIMAL_MAX_AA):
            continue
        parents = getattr(e, "parent_peptide_ids", None) or []
        if (imp_id and imp_id in parents) or (imp_seq and seq and seq in imp_seq):
            hits.append(e)
    return hits


def _ev_get(obj, name, default=None):
    if isinstance(obj, dict):
        return obj.get(name, default)
    return getattr(obj, name, default)


def _ev_hay(ev) -> str:
    return " ".join(
        t for t in (
            _ev_get(ev, "quoted_text"),
            _ev_get(ev, "assay_detail"),
            _ev_get(ev, "section_ref"),
        ) if t
    )


def _quote_names_nested_ept(ev, nested: list) -> bool:
    """True iff the quote uniquely names one class-I bite in `nested` (tests / H1)."""
    return _unique_named_class_i_epitope(ev, nested) is not None


def _named_class_i_hits(ev, class_i: list) -> list:
    hay = _ev_hay(ev)
    if not hay:
        return []
    hay_u = hay.upper()
    hits = []
    for e in class_i:
        seq = (_ev_get(e, "sequence") or "").upper()
        if seq and re.search(r"\b" + re.escape(seq) + r"\b", hay_u):
            hits.append(e)
            continue
        pid = _ev_get(e, "paper_local_id") or ""
        if pid and re.search(r"\b" + re.escape(pid) + r"\b", hay):
            hits.append(e)
    ept_toks = [m.group(0).upper().replace(" ", "") for m in _NAMED_EPT_ID.finditer(hay)]
    if ept_toks and not hits:
        for tok in ept_toks:
            for e in class_i:
                blob = " ".join(
                    t for t in (_ev_get(e, "paper_local_id"), _ev_get(e, "quoted_text")) if t
                ).upper().replace(" ", "")
                if tok and tok in blob:
                    hits.append(e)
    return hits


def _unique_named_class_i_epitope(ev, class_i: list):
    """Exactly one resolvable class-I 8-11mer from the row's quote. Zero names
    (Sahin/SLX4) and two+ distinct named EPT ids (Ott COL22A1 three EPTs) → None.

    Two+ EPT ids win even when only one of those epitopes is still in the catalog
    (stick-outs quarantined): do not collapse the CD8 row onto the one minted 9-mer.
    A contrastive quote that names two EPT ids plus one sequence therefore also
    stays on the IMP — same rule, documented limitation, not a sequence override.
    """
    hay = _ev_hay(ev)
    ept_toks = {m.group(0).upper().replace(" ", "") for m in _NAMED_EPT_ID.finditer(hay)}
    if len(ept_toks) >= 2:
        return None
    hits = _named_class_i_hits(ev, class_i)
    uniq = {}
    for e in hits:
        key = (_ev_get(e, "sequence") or "").upper() or (_ev_get(e, "paper_local_id") or "")
        uniq[key] = e
    if len(uniq) != 1:
        return None
    return next(iter(uniq.values()))


def _cd8_long_peptide_ept_gap(rec: "ExtractedPaper") -> list[str]:
    """Fire only when a UNIQUE named class-I bite is resolvable and the row is still
    on the long IMP (H1 should have re-pointed). Several named EPTs → stay on IMP."""
    imps = {getattr(p, "paper_local_id", None): p for p in (rec.immunizing_peptides or [])}
    class_i = [
        e for e in (rec.epitopes or [])
        if getattr(e, "mhc_class", None) == "I"
        and 8 <= len(getattr(e, "sequence", None) or "") <= _CLASS_I_MINIMAL_MAX_AA
        and getattr(e, "paper_local_id", None)
    ]
    msgs: list[str] = []
    for ev in rec.evidence or []:
        if getattr(ev, "t_cell_subset", None) != "cd8":
            continue
        if getattr(ev, "target_kind", None) != "immunizing_peptide":
            continue
        iid = getattr(ev, "immunizing_peptide_paper_id", None)
        imp = imps.get(iid)
        iseq = (getattr(imp, "sequence", None) or "") if imp is not None else ""
        if len(iseq) <= _CLASS_I_MINIMAL_MAX_AA:
            continue
        chosen = _unique_named_class_i_epitope(ev, class_i)
        if chosen is None:
            continue
        eid = getattr(chosen, "paper_local_id", None) or getattr(chosen, "sequence", "")
        msgs.append(
            f"CD8 evidence targets immunizing peptide {iid!r} (len={len(iseq)}) but the quote "
            f"names a unique class-I 8-11mer {eid!r}. Use target_kind='epitope' and "
            f"epitope_paper_id={eid!r}. Override only if the paper assayed the LONG peptide "
            f"for CD8 despite naming that EPT."
        )
    return msgs


def _repoint_cd8_named_ept(rec: dict) -> int:
    """H1: if a CD8-on-IMP row's quote names exactly one class-I 8-11mer sequence that
    already exists as an epitope, re-point to that epitope (even if canonicalize_evidence
    hung the row on a sibling overlapping IMP). Unnamed (Sahin/SLX4) and ambiguous
    (two named 9-mers) are left alone. Mutates rec; returns rows re-pointed."""
    imps = {p.get("paper_local_id"): p for p in rec.get("immunizing_peptides") or []}
    class_i = [
        e for e in (rec.get("epitopes") or [])
        if e.get("mhc_class") == "I"
        and 8 <= len(e.get("sequence") or "") <= _CLASS_I_MINIMAL_MAX_AA
        and e.get("paper_local_id")
    ]
    changed = 0
    for ev in rec.get("evidence") or []:
        if ev.get("t_cell_subset") != "cd8" or ev.get("target_kind") not in (
            "immunizing_peptide", "candidate"
        ):
            continue
        iid = ev.get("immunizing_peptide_paper_id") or ev.get("candidate_paper_id")
        imp = imps.get(iid) or {}
        iseq = imp.get("sequence") or ""
        # Candidate 8-11mers still get EPT-id / sequence matching; long IMPs skip if already minimal.
        if ev.get("target_kind") == "immunizing_peptide" and len(iseq) <= _CLASS_I_MINIMAL_MAX_AA:
            continue
        chosen = _unique_named_class_i_epitope(ev, class_i)
        if chosen is None:
            continue
        eid = _ev_get(chosen, "paper_local_id")
        ev["target_kind"] = "epitope"
        ev["epitope_paper_id"] = eid
        ev["immunizing_peptide_paper_id"] = None
        ev["candidate_paper_id"] = None
        changed += 1
    return changed


# Overlapping ASP 15-mers dumped as class-II MinimalEpitopes (Ott SI5: 551 windows inside
# ~100 IMPs). Rojas named-DR epitopes are ~1 per IMP (max 4 windows).
_ASP_II_WINDOW_MIN_AA = 12
_ASP_II_WINDOW_MAX_AA = 18
_ASP_II_DUMP_MIN_ROWS = 50
_ASP_II_DUMP_MIN_MAX_PER_PARENT = 6


def _lane(rec, name):
    if isinstance(rec, dict):
        return rec.get(name) or []
    return getattr(rec, name, None) or []


def _mhc_class_is_ii(e) -> bool:
    return str(_ev_get(e, "mhc_class") or "").upper() in ("II", "2", "CLASS_II")


def _has_predicted_affinity(e) -> bool:
    """True iff a predicted_affinity Measurement carries a value or a non-empty raw token.

    A predicted HLA allele alone does NOT count — Ott SI5 stamps DP/DR on every sliding
    15-mer. Named class-II findings (Sahin/Hu/Ott tetramers, Rojas DR) store a score.
    """
    pa = _ev_get(e, "predicted_affinity")
    if not pa:
        return False
    if isinstance(pa, dict):
        if pa.get("value") is not None:
            return True
        return bool((pa.get("raw") or "").strip())
    if getattr(pa, "value", None) is not None:
        return True
    return bool((getattr(pa, "raw", None) or "").strip())


def _seq_contained_in_an_imp(seq, imps) -> bool:
    if not seq:
        return False
    for p in imps:
        iseq = _ev_get(p, "sequence") or ""
        if iseq and seq in iseq and seq != iseq:
            return True
    return False


def _evidence_epitope_ids(rec) -> set:
    return {
        _ev_get(e, "epitope_paper_id")
        for e in _lane(rec, "evidence")
        if _ev_get(e, "target_kind") == "epitope" and _ev_get(e, "epitope_paper_id")
    }


def _is_droppable_asp_ii_window(e, imps, ev_epi: set) -> bool:
    """Class-II 12-18mer inside an IMP, no evidence, no affinity score.

    Predicted HLA on the window is not a keep — that is the Ott SI5 dump shape.
    """
    if not _mhc_class_is_ii(e):
        return False
    seq = _ev_get(e, "sequence") or ""
    if not (_ASP_II_WINDOW_MIN_AA <= len(seq) <= _ASP_II_WINDOW_MAX_AA):
        return False
    pid = _ev_get(e, "paper_local_id")
    if pid and pid in ev_epi:
        return False
    if _has_predicted_affinity(e):
        return False
    return _seq_contained_in_an_imp(seq, imps)


def _drop_class_ii_asp_windows(rec: dict) -> int:
    """Remove ASP-window class-II epitopes (and refs to them). Mutates rec; returns n dropped."""
    imps = rec.get("immunizing_peptides") or []
    ev_epi = _evidence_epitope_ids(rec)
    drop_ids = {
        e.get("paper_local_id")
        for e in (rec.get("epitopes") or [])
        if e.get("paper_local_id") and _is_droppable_asp_ii_window(e, imps, ev_epi)
    }
    if not drop_ids:
        return 0
    rec["epitopes"] = [e for e in (rec.get("epitopes") or [])
                       if e.get("paper_local_id") not in drop_ids]
    for lane in _EPITOPE_REF_LANES:
        rows = rec.get(lane) or []
        rec[lane] = [r for r in rows
                     if not (isinstance(r, dict) and r.get("epitope_paper_id") in drop_ids)]
    for note in rec.get("curator_notes") or []:
        if isinstance(note, dict) and note.get("refs"):
            note["refs"] = [r for r in note["refs"] if r not in drop_ids]
    return len(drop_ids)


def _asp_ii_window_parents(e, imps) -> list[str]:
    seq = _ev_get(e, "sequence") or ""
    parents = []
    if not (_ASP_II_WINDOW_MIN_AA <= len(seq) <= _ASP_II_WINDOW_MAX_AA):
        return parents
    for p in imps:
        iseq = _ev_get(p, "sequence") or ""
        if seq and iseq and seq in iseq and seq != iseq:
            pid = _ev_get(p, "paper_local_id")
            if pid:
                parents.append(pid)
    return parents


def _class_ii_asp_dump_gap(rec: "ExtractedPaper") -> list[str]:
    """O2 leftovers: overlapping class-II 12-18mers inside IMPs that the auto-drop
    did not remove (they have an affinity score) and that no evidence targets.

    Auto-droppable windows (no affinity, no evidence) are not counted here.
    """
    imps = _lane(rec, "immunizing_peptides")
    if not imps:
        return []
    ev_epi = _evidence_epitope_ids(rec)
    per_parent: dict[str, int] = {}
    n_windows = 0
    for e in _lane(rec, "epitopes"):
        if not _mhc_class_is_ii(e):
            continue
        if _is_droppable_asp_ii_window(e, imps, ev_epi):
            continue
        if _ev_get(e, "paper_local_id") in ev_epi:
            continue
        parents = _asp_ii_window_parents(e, imps)
        if not parents:
            continue
        n_windows += 1
        for pid in parents:
            per_parent[pid] = per_parent.get(pid, 0) + 1
    if n_windows < _ASP_II_DUMP_MIN_ROWS:
        return []
    mx = max(per_parent.values()) if per_parent else 0
    if mx < _ASP_II_DUMP_MIN_MAX_PER_PARENT:
        return []
    return [
        f"{n_windows} class-II 12-18mer epitopes with predicted affinity are overlapping "
        f"windows inside immunizing peptides (max {mx} per IMP) with no evidence targeting "
        f"them. Those look like an ASP/prediction dump that survived auto-drop — CD4 stays "
        f"on the IMP; mint class-II only when the paper names a DR/DP/DQ (or I-A/I-E) "
        f"restriction for a defined minimal sequence. Drop the extra window rows. Override "
        f"only if they really are named class-II minimal determinants."
    ]


# minimal-epitope GRAIN anchor (v2.12 P23 follow-up): a MinimalEpitope is the MINIMAL binder — a class-I
# epitope is an 8-11mer (the whole gold corpus class-I epitopes are 8-12 aa). Two failure modes route the
# record to needs_review (both seen in the lite-lane Li 33879241 run, which produced 16 "class-I epitopes"
# that were actually 19-29mer LONG peptides and ZERO of the 50 real minimal epitopes):
#   (a) a class-I epitope is >14 aa -> it is almost certainly a long IMMUNIZING peptide mislabeled as an
#       epitope (move it to immunizing_peptides; mint the predicted minimal epitope instead); and
#   (b) the paper PREDICTS minimal epitopes (a NetMHC/affinity table: predicted_affinity set, or an
#       IC50/%rank/'minimal epitope' cue) yet ZERO class-I 8-11mer minimal epitopes were minted across a
#       real (>=6) immunizing-peptide set -> the entire minimal-epitope layer was dropped.
# SOFT/block-once, overridable allow_missing_minimal_epitopes (HARD -> needs_review). companion_paper_ref
# (manifest deferred to a prior paper) exempts side (b). Calibrated: 0-fire on keskin/rojas/li gold.
_PRED_CUE = re.compile(r"netmhc|ic50|%?\s?rank|minimal epitope|predicted (minimal )?epitope|mhc[- ]?i prediction", re.I)


def _minimal_epitope_grain_gap(rec: "ExtractedPaper") -> list[str]:
    # v2.13 (Fable review #1): the "class-I epitope >14 aa is a mislabeled long peptide" check is now a
    # HARD schema reject (MinimalEpitope._class_i_is_a_minimal_binder) — a valid record can't reach here
    # with one. This soft nudge keeps only side (b): the DROPPED minimal-epitope LAYER (a prediction table
    # was read but no minimal epitope minted). Companion-deferred manifests are exempt.
    epis = rec.epitopes or []
    n_minimal_ci = sum(1 for e in epis
                       if getattr(e, "mhc_class", None) == "I" and 8 <= len(e.sequence or "") <= 11)
    n_imp = len(rec.immunizing_peptides or [])
    companion = getattr(rec, "companion_paper_ref", None)
    has_pred_cue = any(getattr(e, "predicted_affinity", None) for e in epis) or any(
        _PRED_CUE.search(t) for t in (
            *[e.quoted_text for e in epis if getattr(e, "quoted_text", None)],
            *[getattr(e, "section_ref", "") or "" for e in epis],
        ))
    msgs = []
    if not companion and has_pred_cue and n_minimal_ci == 0 and n_imp >= 6:
        msgs.append(
            f"the paper predicts minimal epitopes (affinity/NetMHC cues present) but ZERO class-I 8-11mer "
            f"minimal epitopes were minted across {n_imp} immunizing peptides -- the minimal-epitope layer "
            f"was dropped; mint the predicted class-I minimal epitope per peptide")
    return msgs


# Fields that do NOT define an epitope's scientific identity — they are bookkeeping
# (the label, the many-to-many parent links) or provenance, so two records that differ
# ONLY in these are the SAME epitope recorded twice and must be merged.
_EPITOPE_ID_EXCLUDE = frozenset({
    "paper_local_id", "parent_peptide_ids", "quoted_text", "section_ref",
    "provenance", "confidence", "needs_review",
})


def _epitope_identity(e: dict) -> str:
    """Canonical key for an epitope's SCIENTIFIC identity (everything except label/parents/
    provenance). Same key => the same epitope; different HLA / affinity / gene / mutation =>
    different key (so Rojas same-sequence-different-restriction records stay separate)."""
    core = {k: v for k, v in e.items() if k not in _EPITOPE_ID_EXCLUDE}
    return json.dumps(core, sort_keys=True, default=str)


def canonicalize_epitopes(rec: dict) -> int:
    """Collapse MinimalEpitope records that are IDENTICAL on every scientific field and differ
    ONLY in paper_local_id / parent_peptide_ids into ONE record holding the UNION of their
    parent_peptide_ids — the schema's intended MANY-TO-MANY form (one minimal epitope tiling
    several long peptides is one epitope with several parents, not N duplicate rows). Records
    differing in ANY scientific field (e.g. same sequence but a different HLA restriction) are
    LEFT SEPARATE. Evidence epitope_paper_id and curator_note refs pointing at a dropped id are
    remapped to the surviving id. Deterministic (smallest paper_local_id survives; parents
    sorted) and idempotent. Mutates `rec`; returns the number of records merged away."""
    eps = rec.get("epitopes") or []
    if not eps:
        return 0
    groups: dict = {}
    order: list = []
    for e in eps:
        k = _epitope_identity(e)
        if k not in groups:
            groups[k] = []
            order.append(k)
        groups[k].append(e)
    canon: list = []
    remap: dict = {}
    merged = 0
    for k in order:
        members = groups[k]
        if len(members) == 1:
            canon.append(members[0])
            continue
        survivor = dict(min(members, key=lambda e: e.get("paper_local_id") or ""))
        survivor["parent_peptide_ids"] = sorted(
            {p for e in members for p in (e.get("parent_peptide_ids") or [])}
        )
        if any(e.get("needs_review") for e in members):
            survivor["needs_review"] = True
        for e in members:
            pid = e.get("paper_local_id")
            if pid and pid != survivor.get("paper_local_id"):
                remap[pid] = survivor["paper_local_id"]
        merged += len(members) - 1
        canon.append(survivor)
    rec["epitopes"] = canon
    if remap:
        for ev in rec.get("evidence") or []:
            if ev.get("epitope_paper_id") in remap:
                ev["epitope_paper_id"] = remap[ev["epitope_paper_id"]]
        for note in rec.get("curator_notes") or []:
            note["refs"] = [remap.get(r, r) for r in (note.get("refs") or [])]
    return merged


# --- ORPHANED-EPITOPE REPAIR (2026-08-16) -------------------------------------------------------------
# ROOT CAUSE (tcr_phase5 live runs): `add_epitope_manifest` -> epitope_manifest_adapter emits manifest
# epitopes with NO parent_peptide_ids at all, and then LOCKS the epitopes lane. So every manifest paper
# produced a 100%-orphaned epitope lane that the outer_guard rejects and that the agent CANNOT repair
# (clear_entities is refused on a locked lane; add_entities only appends). 3 of 4 phase-5 papers died
# here (33479501 203/203, 37165196 225/682 loader-loaded, 30568305 99/99); the one paper that never
# called the loader (39972124) was clean, and the whole phase-4 generation -- which predates the loader
# -- has ZERO orphans. So this is a CODE defect in the loader lane, not agent misbehaviour.
#
# THE REPAIR: a minimal epitope is by definition a SUBSTRING of the long peptide that was immunized, so
# the parent link is DERIVABLE from the record itself -- the same sequence-containment rule
# `build_crossreactivity_evidence` already uses. We never invent an id: a parent must be an
# immunizing peptide already in the record whose sequence CONTAINS the epitope's. Ambiguity is resolved
# by gene_symbol when both sides declare one; a genuinely-tiled epitope keeps ALL its containing parents
# (schema.MinimalEpitope.parent_peptide_ids is many-to-many by design).
_MAX_PARENT_PEPTIDES = 12   # schema.MinimalEpitope.parent_peptide_ids max_length


def link_orphan_epitopes(rec: dict) -> int:
    """Give every ORPHANED epitope (empty parent_peptide_ids) the immunizing peptide(s) whose sequence
    CONTAINS it. Never touches an epitope that already has parents, never mints an id, and never guesses:
    an epitope with no containing peptide stays orphaned (quarantine over fabrication).

    Disambiguation, in order: (1) candidates = loaded immunizing peptides containing the epitope
    sequence; (2) if the epitope and >=1 candidate agree on gene_symbol, keep only the gene-matched
    ones; (3) more than `_MAX_PARENT_PEPTIDES` survivors means the link is not pinned down by the
    source -- leave it orphaned rather than truncate an arbitrary 12. Deterministic (ids sorted) and
    idempotent. Mutates `rec`; returns the number of epitopes newly linked."""
    eps = rec.get("epitopes") or []
    imps = rec.get("immunizing_peptides") or []
    if not eps or not imps:
        return 0
    catalog = [(_norm_seq(p.get("sequence")), (p.get("gene_symbol") or "").strip().upper(),
                p.get("paper_local_id"))
               for p in imps if p.get("sequence") and p.get("paper_local_id")]
    linked = 0
    for e in eps:
        if e.get("parent_peptide_ids"):
            continue
        seq = _norm_seq(e.get("sequence"))
        if not seq:
            continue
        hits = [(gene, pid) for pseq, gene, pid in catalog if seq in pseq]
        if not hits:
            continue
        gene = (e.get("gene_symbol") or "").strip().upper()
        same_gene = [pid for g, pid in hits if gene and g == gene]
        ids = sorted({pid for pid in (same_gene or [pid for _, pid in hits])})
        if len(ids) > _MAX_PARENT_PEPTIDES:
            continue                     # not pinned down by the source -> stays orphaned
        e["parent_peptide_ids"] = ids
        linked += 1
    return linked


# An orphan that survives `link_orphan_epitopes` cannot be joined to any peptide, so it is not usable
# data -- but it is also not a licence to invent a parent. It is QUARANTINED: pulled out of the record
# (together with the rows that reference it, which carry an exactly-one-target discriminator and so
# cannot simply be re-pointed) and handed back for the caller to write to a sidecar. Bounded on purpose:
# quarantining a LARGE share of the lane would mean the peptide lane itself is wrong (wrong sheet, wrong
# column) and must be fixed, not swept away -- above the cap we refuse and report instead.
_MAX_QUARANTINE_FRACTION = 0.10
_EPITOPE_REF_LANES = ("evidence", "screening_readouts", "neoantigen_tcr_flags", "tcr_clonotypes")


def quarantine_orphan_epitopes(rec: dict) -> tuple[dict, str | None]:
    """Remove epitopes still orphaned after `link_orphan_epitopes`, plus every row pointing at them.

    Returns ({"epitopes": [...], <lane>: [...]}, None) with the removed rows (empty dict = nothing to
    do), or ({}, error) when the orphan share exceeds `_MAX_QUARANTINE_FRACTION` (or nothing linked at
    all) -- a systematic load error the curator must fix, NOT row-level source noise. Mutates `rec`."""
    eps = rec.get("epitopes") or []
    # An epitope that evidence already targets is not dump noise: Hu Supp11b's
    # assayed 15-mers are not nested in the immunizing peptides. Empty parents
    # are the schema's honest "paper did not resolve a long peptide" case.
    # Quarantining them deletes the measurements. Keep evidence-linked orphans.
    ev_epi = {
        r.get("epitope_paper_id")
        for r in rec.get("evidence") or []
        if isinstance(r, dict) and r.get("epitope_paper_id")
    }
    orphan_ids = {e.get("paper_local_id") for e in eps
                  if not e.get("parent_peptide_ids")
                  and e.get("paper_local_id") not in ev_epi}
    orphan_ids.discard(None)
    if not orphan_ids:
        return {}, None
    if (len(orphan_ids) == len(eps)                                     # never empty the whole lane
            or len(orphan_ids) > max(1, int(len(eps) * _MAX_QUARANTINE_FRACTION))):
        return {}, (
            f"[gate] epitope linkage broken: {len(orphan_ids)}/{len(eps)} epitopes have EMPTY "
            f"parent_peptide_ids AND their sequences are not contained in ANY loaded immunizing "
            f"peptide, so no parent can be derived. That share is too large to be source noise -- the "
            f"immunizing_peptides lane is probably incomplete or read from the wrong sheet/column. Load "
            f"the peptide manifest these epitopes came from (add_epitope_manifest / add_table on the "
            f"long-peptide column), then finalize again. Do NOT invent parent ids. "
            f"First orphans: {sorted(orphan_ids)[:3]}"
        )
    removed: dict = {"epitopes": [e for e in eps if e.get("paper_local_id") in orphan_ids]}
    rec["epitopes"] = [e for e in eps if e.get("paper_local_id") not in orphan_ids]
    for lane in _EPITOPE_REF_LANES:
        rows = rec.get(lane) or []
        hit = [r for r in rows if isinstance(r, dict) and r.get("epitope_paper_id") in orphan_ids]
        if hit:
            removed[lane] = hit
            rec[lane] = [r for r in rows if r not in hit]
    for note in rec.get("curator_notes") or []:
        if isinstance(note, dict) and note.get("refs"):
            note["refs"] = [r for r in note["refs"] if r not in orphan_ids]
    return removed, None


# ---------------------------------------------------------------------------
# Evidence-lane determinism (stability lab, 2026-07-06). Diagnosis on Keskin: the 3 live runs recorded
# the IDENTICAL biology (Pt7/Pt8 responses to ARHGAP35/GPC1/SLX4/SHANK2/SVEP1/COX18) but scored
# stability 0.0 because they encoded it three different ways: patient id "7" vs "Pt7", and the SAME
# response attached to a candidate in one run, an immunizing peptide in another. Zero hallucinations —
# pure keying noise. These three deterministic finalize passes canonicalize that encoding so the runs
# agree, at the grain the curator chose: (patient, mutated-gene, mhc_class, assay).
# ---------------------------------------------------------------------------
_PATIENT_ALIAS = re.compile(r"^(?:patients?|pts?|subjects?|cases?|p)?[\s._\-#]*(\d+)$", re.I)


def canon_patient_id(pid):
    """Canonicalize a numeric patient alias to 'Pt<N>' ('7', 'Patient 7', 'P-7', 'pt_7' -> 'Pt7');
    leave non-numeric aliases (e.g. 'MEL-21', a mouse strain) untouched. Deterministic + idempotent.

    FORMAT CHOICE ('Pt<N>'): the live agent emits the SAME patient in several formats run-to-run ('7'
    vs 'Pt7'), which splits one fact into two metric rows — so a single global form is required for
    run-to-run stability, and no per-paper form can (the entity format itself wobbles). 'Pt<N>' matches
    the Keskin gold exactly. KNOWN DIVERGENCE: the Rojas (37165196) and Li (33879241) reference_records
    use 'P<N>', so once those papers are run live their patient-keyed correctness will read low against
    the current golds purely on this cosmetic axis. Resolve by canonicalizing those golds with THIS
    function (they are blessed-LLM records, so the format is not authoritative) before comparing."""
    if not isinstance(pid, str):
        return pid
    m = _PATIENT_ALIAS.match(pid.strip())
    return f"Pt{int(m.group(1))}" if m else pid.strip()


# Lanes carrying a `patient_paper_id` reference (the patient's OWN id lives on ExtractedPatient
# .paper_local_id, handled separately). Kept explicit so a new lane is a deliberate add.
_PATIENT_REF_LANES = ("immunizing_peptides", "epitopes", "candidates", "pools", "evidence",
                      "screening_readouts", "neoantigen_tcr_flags", "tcr_clonotypes")


def canonicalize_patient_ids(rec: dict) -> int:
    """FIX 1. Rewrite every patient id to canonical 'Pt<N>' form — the patient's own paper_local_id
    AND every patient_paper_id reference — using one pure function, so definitions and references stay
    aligned (the metric keys evidence on patient id, so '7' vs 'Pt7' split the same fact into two rows).
    Mutates `rec`; returns the count of ids rewritten."""
    n = 0
    for p in rec.get("patients") or []:
        old = p.get("paper_local_id"); new = canon_patient_id(old)
        if new != old:
            p["paper_local_id"] = new; n += 1
    for lane in _PATIENT_REF_LANES:
        for x in rec.get(lane) or []:
            old = x.get("patient_paper_id")
            if old is None:
                continue
            new = canon_patient_id(old)
            if new != old:
                x["patient_paper_id"] = new; n += 1
    return n


def _gene_to_canonical_imp(rec: dict) -> dict:
    """gene_symbol -> the deterministic representative immunizing-peptide (paper_local_id, sequence):
    lexicographically smallest sequence, so every run picks the same one."""
    by: dict = {}
    for x in rec.get("immunizing_peptides") or []:
        g, s, pid = x.get("gene_symbol"), x.get("sequence"), x.get("paper_local_id")
        if g and s and pid:
            by.setdefault(g, []).append((s, pid))
    return {g: sorted(v)[0] for g, v in by.items()}


def _target_id_to_gene(rec: dict) -> dict:
    m: dict = {}
    for lane in ("immunizing_peptides", "epitopes", "candidates"):
        for x in rec.get(lane) or []:
            if x.get("paper_local_id"):
                m[x["paper_local_id"]] = x.get("gene_symbol")
    return m


def _evidence_measurement_key(e: dict) -> tuple:
    """FULL measurement identity — two evidence rows with this key are the same fact.

    `assay_antigen_format` and `assay_detail` are part of it: a reactivity matrix measures the
    same patient/epitope/phase/stimulation peptide-pulsed, by minigene, against autologous
    tumour, and again with the response antibody-blocked. Leaving them out collapses those
    distinct measurements onto one row (16 lost on Hu 33479501 'Suppl Dataset 4b' alone)."""
    return (e.get("patient_paper_id"), e.get("target_kind"), e.get("immunizing_peptide_paper_id"),
            e.get("epitope_paper_id"), e.get("pool_paper_id"), e.get("candidate_paper_id"),
            e.get("assay"), e.get("outcome"), e.get("t_cell_subset"), e.get("mhc_class"),
            e.get("assay_stimulation"), e.get("timepoint_phase"),
            e.get("assay_antigen_format"), e.get("assay_detail"))


def canonicalize_evidence(rec: dict) -> int:
    """FIX 2. Put immunogenic evidence at the canonical grain (patient, mutated-gene, mhc_class, assay).
    The response to an ADMINISTERED neoantigen belongs on its immunizing peptide (schema §NeoantigenCandidate
    — a selected candidate's outcome rides its bridged IMP, not the funnel entity). So any immunogenic
    `candidate`- or `immunizing_peptide`-target row is re-pointed to its gene's DETERMINISTIC representative
    immunizing peptide — making the runs agree even under the sequence-keyed metric (epitopes are already
    run-invariant via the manifest loader). `epitope`- and `pool`-target rows are LEFT ALONE (an
    epitope-specific assay is legitimately epitope-grained; pools are Fix 3's job). Identical measurements
    produced by the re-point are de-duplicated (full measurement identity, so CD4 vs CD8 stay separate).
    Mutates `rec`; returns rows re-pointed + rows de-duplicated."""
    ev = rec.get("evidence") or []
    if not ev:
        return 0
    g2i = _gene_to_canonical_imp(rec)
    i2g = _target_id_to_gene(rec)
    changed = 0
    for e in ev:
        if e.get("outcome") != "immunogenic" or e.get("target_kind") not in ("candidate", "immunizing_peptide"):
            continue
        tid = e.get("candidate_paper_id") or e.get("immunizing_peptide_paper_id")
        gene = i2g.get(tid)
        if not gene or gene not in g2i:
            continue                                  # can't resolve a canonical IMP -> leave honest
        _, pid = g2i[gene]
        if e.get("target_kind") == "immunizing_peptide" and e.get("immunizing_peptide_paper_id") == pid:
            continue                                  # already canonical
        e["target_kind"] = "immunizing_peptide"
        e["immunizing_peptide_paper_id"] = pid
        e["candidate_paper_id"] = e["epitope_paper_id"] = e["pool_paper_id"] = None
        changed += 1
    seen, out, dropped = set(), [], 0
    for e in ev:
        key = _evidence_measurement_key(e)
        if key in seen:
            dropped += 1
            continue
        seen.add(key)
        out.append(e)
    rec["evidence"] = out
    return changed + dropped


def derive_pool_evidence(rec: dict) -> int:
    """FIX 3. Roll up a pool-level immunogenic row per (patient, pool) that the record already ENTAILS:
    the patient has an immunogenic response to a peptide that is a MEMBER of that pool, yet no pool-target
    row records it (the recall gap — Keskin gold carries these coarse per-pool summaries; the runs miss
    them). This is a faithful roll-up of evidence already present (never invention): it fires only where a
    member response exists, references an EXISTING pool, and is flagged needs_review + magnitude.raw (so it
    doesn't trip the magnitude nudge). Idempotent. Mutates `rec`; returns rows added."""
    pools = rec.get("pools") or []
    ev = rec.get("evidence") or []
    if not pools or not ev:
        return 0
    id2seq: dict = {}
    for lane in ("immunizing_peptides", "epitopes", "candidates"):
        for x in rec.get(lane) or []:
            if x.get("paper_local_id") and x.get("sequence"):
                id2seq[x["paper_local_id"]] = x["sequence"]
    have = {(e.get("patient_paper_id"), e.get("pool_paper_id"))
            for e in ev if e.get("outcome") == "immunogenic" and e.get("target_kind") == "pool"}
    resp, assay_of = {}, {}
    for e in ev:
        if e.get("outcome") != "immunogenic":
            continue
        tid = (e.get("immunizing_peptide_paper_id") or e.get("epitope_paper_id")
               or e.get("candidate_paper_id"))
        s = id2seq.get(tid)
        if s:
            resp.setdefault(e.get("patient_paper_id"), set()).add(s)
            assay_of.setdefault(e.get("patient_paper_id"), e.get("assay"))
    new = []
    for pool in pools:
        pt, pid = pool.get("patient_paper_id"), pool.get("paper_local_id")
        if not pt or not pid or (pt, pid) in have:
            continue
        members = {id2seq.get(m) for m in (pool.get("member_peptide_ids") or [])}
        if not (resp.get(pt) and (resp[pt] & members)):
            continue                                  # patient did not respond to any member of THIS pool
        new.append({
            "patient_paper_id": pt, "target_kind": "pool", "pool_paper_id": pid,
            "assay": assay_of.get(pt) or "elispot", "outcome": "immunogenic", "vaccine_induced": True,
            "magnitude": {"unit": "unknown", "raw": "pool-level; per-value not individually tabulated (roll-up)"},
            "needs_review": True, "confidence": 2,
            "assay_detail": "pool-level response rolled up from this patient's member-peptide responses -- needs_review",
            "quoted_text": f"Patient {pt}: neoantigen-specific response(s) to peptide(s) within this pool (roll-up of member evidence).",
            "section_ref": "derived: member-evidence roll-up",
            "provenance": [{"kind": "prose", "locator": "member-evidence roll-up"}],
        })
    if new:
        rec.setdefault("evidence", []).extend(new)
    return len(new)


def _seq_to_specific_entity(rec: dict) -> dict:
    """sequence -> (target_kind, paper_local_id) for the SPECIFIC entity a re-listed candidate duplicates:
    an epitope (short/minimal) or an immunizing peptide (long/administered). Epitopes first, min id on
    ties, so the mapping is deterministic."""
    out: dict = {}
    for kind, lane in (("epitope", "epitopes"), ("immunizing_peptide", "immunizing_peptides")):
        for x in sorted(rec.get(lane) or [], key=lambda e: e.get("paper_local_id") or ""):
            s, pid = x.get("sequence"), x.get("paper_local_id")
            if s and pid:
                out.setdefault(s, (kind, pid))
    return out


def drop_administered_candidates(rec: dict) -> int:
    """FIX (candidates lane, 2026-07-06). Keskin/Rojas/Li gold record ZERO candidates: an ADMINISTERED
    neoantigen lives in the immunizing_peptides (and, minimally, epitopes) lanes — it is not ALSO a funnel
    `candidate` (a candidate is the UNSELECTED prediction denominator). The live agent nonetheless re-lists
    administered peptides as candidates, swinging 6/0/6 run-to-run (all false positives vs gold=0). This
    drops exactly those re-listed candidates: a candidate whose sequence already appears as an epitope or
    immunizing peptide. Any evidence/screening/curator ref to a dropped candidate is RE-POINTED to that
    same-sequence entity first (so nothing dangles — the paper-level validator rejects an orphaned
    candidate ref), then the candidate is removed. A candidate whose sequence is NOT in either lane (a
    genuinely unselected prediction) is KEPT. Idempotent; mutates `rec`; returns candidates dropped.

    NOTE: no gold in the corpus carries a real unselected funnel, so this is verified not to REGRESS
    (all three golds stay at 0) but can't be positively confirmed to preserve a real funnel — the
    sequence-not-in-either-lane guard is the safeguard, and unselected candidates never match an
    administered peptide/epitope by construction."""
    cands = rec.get("candidates") or []
    if not cands:
        return 0
    seq2ent = _seq_to_specific_entity(rec)
    remap = {c["paper_local_id"]: seq2ent[c["sequence"]]
             for c in cands
             if c.get("paper_local_id") and c.get("sequence") in seq2ent}
    if not remap:
        return 0
    # every lane that can carry candidate_paper_id, with the NAME of its target-kind discriminator
    # (evidence/screening use `target_kind`; a TcrClonotype uses `specific_for_kind`).
    for lane, kind_field in (("evidence", "target_kind"), ("screening_readouts", "target_kind"),
                             ("neoantigen_tcr_flags", "target_kind"), ("tcr_clonotypes", "specific_for_kind")):
        for e in rec.get(lane) or []:
            cid = e.get("candidate_paper_id")
            if cid in remap:
                kind, eid = remap[cid]
                e[kind_field] = kind
                e["candidate_paper_id"] = None
                e["epitope_paper_id"] = eid if kind == "epitope" else None
                e["immunizing_peptide_paper_id"] = eid if kind == "immunizing_peptide" else None
                if lane in ("evidence", "screening_readouts"):
                    e["pool_paper_id"] = None            # pool target exists only on these two
    for note in rec.get("curator_notes") or []:
        if isinstance(note, dict) and note.get("refs"):
            note["refs"] = [remap[r][1] if r in remap else r for r in note["refs"]]
    rec["candidates"] = [c for c in cands if c.get("paper_local_id") not in remap]
    # a re-point can duplicate an existing row (full measurement identity) -> collapse
    seen, out = set(), []
    for e in rec.get("evidence") or []:
        key = _evidence_measurement_key(e)
        if key in seen:
            continue
        seen.add(key)
        out.append(e)
    rec["evidence"] = out
    return len(remap)


# Override tiers (v2.11.5, after the 2026-06-09 iris batch routed ALL 7 papers to needs_review):
# a HARD override marks a genuine recall/correctness gap a human must confirm -> needs_review. A SOFT
# override is expected variance or completeness-METADATA that does not impugn the extracted data, so it
# stays in the clean (extracted/) lane WITH the override still recorded in finalize_overrides_used for
# audit. The scale lane (extract_one.py) routes to needs_review unless EVERY override used is SOFT —
# so an unknown/future override defaults to needs_review (conservative). HARD is everything not SOFT.
SOFT_OVERRIDES = frozenset({
    "allow_unknown_funnel_size",   # n_predicted_reported unset — funnel-completeness metadata, not a data error
    "allow_regimen_divergence",    # per-arm delivery covariate inconsistency — a covariate, not the core readout
})
HARD_OVERRIDES = frozenset({
    "allow_missing_magnitudes", "allow_missing_pools", "allow_member_level_pool_evidence",
    "allow_candidate_bridge_mismatch", "allow_evidence_count_mismatch",
    "allow_peptide_count_mismatch", "allow_sparse_evidence", "allow_missing_class_ii",
    "allow_missing_minimal_epitopes", "allow_ungrounded_safety_grade",
    "allow_missing_tcr_gateway",     # v2.17 §5.7: a declared-but-unfilled Tier-A TCR gateway gap
    "allow_tcr_status_mismatch",     # v2.17 §5.3: status sticker insisted-on despite mismatched TCR content
    "allow_missing_tcr_flags",       # the FLAGSHIP per-(neoantigen x patient) TCR flag lane left empty
    "allow_unloaded_tcr_sheets",     # loadable TCR sheets the source offered and the record never read
    "allow_underloaded_read_sheets",  # a sheet filtered-read to N rows, then barely loaded
    "allow_underloaded_reactivity_sheets",  # a reactivity matrix the source offers, barely loaded
    "allow_cd4_on_class_i_minimal",  # CD4 evidence pointed at a class-I 8-11mer instead of the long IMP
    "allow_cd8_on_immunizing_peptide",  # CD8/EPT evidence pointed at a long IMP instead of the 9-mer
    "allow_class_ii_asp_dump",  # overlapping ASP 15-mers minted as class-II epitopes (Ott SI5)
})


def overrides_are_soft_only(overrides) -> bool:
    """True iff `overrides` is non-empty and EVERY entry is a SOFT override (so the record may stay in
    the clean lane). Empty -> False (nothing to tier). Any non-SOFT (incl. unknown) -> False -> the
    scale lane sends it to needs_review."""
    ov = list(overrides or [])
    return bool(ov) and all(o in SOFT_OVERRIDES for o in ov)


def _screening_target_key(d: dict):
    tgt = (d.get("immunizing_peptide_paper_id") or d.get("epitope_paper_id")
           or d.get("pool_paper_id") or d.get("candidate_paper_id") or "")
    return (d.get("patient_paper_id"), tgt, d.get("assay"))


def _drop_screening_covered_by_evidence(rec: dict) -> int:
    """v2.14 #4 PRECEDENCE de-dup: a NAMED ExtractedEvidence row is the real claim, so a ScreeningReadout
    row at the same (patient, target, assay) is a duplicate of the screening denominator and is DROPPED
    (the named result supersedes — a target that became a real positive/negative must not also sit in the
    screening bucket). Deterministic (root-cause), like canonicalize_epitopes. Returns count dropped."""
    scr = rec.get("screening_readouts") or []
    if not scr:
        return 0
    ev_keys = {_screening_target_key(e) for e in (rec.get("evidence") or [])}
    kept = [s for s in scr if _screening_target_key(s) not in ev_keys]
    dropped = len(scr) - len(kept)
    if dropped:
        rec["screening_readouts"] = kept
    return dropped


# v2.15 safety-grounding guard (root-caused on 39041242 / MVX-ONCO-1): the agent mapped a *serious*
# adverse event (SAE = a REGULATORY seriousness category) onto a grade>=3 CTCAE claim, and let a
# "severe/life-threatening" descriptor that the paper attributed to DISEASE PROGRESSION drive
# any_grade3plus_related=true. Seriousness != grade != relatedness. A grade>=3 RELATED claim must be
# grounded in the verbatim `raw`: if `raw` carries no explicit grade>=3 token AND it instead cites an
# SAE/serious wording or an attribution disclaimer (disease progression / not treatment-related), the
# structured booleans contradict their own quote -> route to needs_review for a human grade check.
_SAFETY_GRADE3_RE = re.compile(r"grades?\s*(?:[≥>]=?\s*)?(?:3|4|5|iii|iv|v\b|three|four|five)", re.I)
_SAFETY_DISCLAIMER_RE = re.compile(
    r"disease progression|progression of disease|not[\s-]*(?:treatment|drug|vaccine|study[\s-]*drug)[\s-]*related"
    r"|not related to (?:treatment|the study|the vaccine|the drug)|unrelated to (?:treatment|the study)"
    r"|regardless of attribution|not (?:treatment|drug|vaccine)-related", re.I)
_SAFETY_SAE_TERM_RE = re.compile(r"serious adverse event|\bSAEs?\b", re.I)
_SAFETY_SEVERITY_RE = re.compile(r"severe|life[\s-]?threatening|fatal|death|\bgrade", re.I)


def _safety_grade_ungrounded(full) -> str | None:
    """Return a reason if safety_summary asserts a grade>=3 TREATMENT-RELATED event that its own `raw`
    quote does not support (no explicit grade>=3) and that the quote actively undercuts (cites an
    SAE/'serious' wording, an attribution disclaimer, or no severity word at all). None = grounded/ok."""
    ss = getattr(full, "safety_summary", None)
    if ss is None:
        return None
    mg = ss.max_related_grade
    claims_g3 = (ss.any_grade3plus_related is True) or (mg is not None and mg >= 3)
    if not claims_g3:
        return None
    raw = ss.raw or ""
    if _SAFETY_GRADE3_RE.search(raw):
        return None  # explicitly grounded in a numeric grade>=3 token -> trust it
    sae = bool(_SAFETY_SAE_TERM_RE.search(raw))
    disclaimer = bool(_SAFETY_DISCLAIMER_RE.search(raw))
    no_severity = not _SAFETY_SEVERITY_RE.search(raw)
    if not (sae or disclaimer or no_severity):
        return None  # 'severe'/'life-threatening' present, no SAE/disclaimer cue -> plausibly legit
    why = []
    if sae:
        why.append("its raw cites a 'serious adverse event'/SAE (a regulatory category, NOT a grade)")
    if disclaimer:
        why.append("its raw attributes the severe events to disease progression / not-treatment")
    if no_severity and not (sae or disclaimer):
        why.append("its raw contains no grade/severity wording at all")
    return (f"safety_summary asserts grade>=3 treatment-related (any_grade3plus_related="
            f"{ss.any_grade3plus_related}, max_related_grade={mg}) but no explicit grade>=3 in `raw`; "
            + " and ".join(why))


def _tcr_gateway_gap(rec: "ExtractedPaper") -> str:
    """v2.17 §5.7 QC gateway anchor. The paper asserts TCR-seq was done (tcr_seq_status is a
    non-'none' value) but NEITHER a TcrSeqMethod NOR a DataDeposition was extracted — so the
    Tier-A gateway fact (HOW it was done / WHERE the recoverable data live) is known-missing, not
    legitimately absent. 'mentioned_only' is the honest "done but not extractable" state and still
    fires here (it IS a declared gap): the override routes the record to needs_review. Returns a
    reason string, else '' (no gap). tcr_seq_status None/'none' -> no TCR claim -> no gap."""
    status = getattr(rec, "tcr_seq_status", None)
    if status in (None, "none"):
        return ""
    if not rec.tcr_seq_methods and not rec.data_depositions:
        return (
            f"tcr_seq_status={status!r} asserts TCR-seq was done but no TcrSeqMethod and no "
            f"DataDeposition were extracted (the Tier-A 'how / where the data live' gateway fact)"
        )
    return ""


def _tcr_flag_gap(rec: "ExtractedPaper") -> str:
    """FLAGSHIP anchor (2026-08-16), symmetric with `_tcr_gateway_gap`. `neoantigen_tcr_flags` -- the
    per-(neoantigen x patient) 'was a TCR identified, and CD4 or CD8?' cell collaborators actually query
    (schema §3d) -- is the ONE TCR lane no gate enforced. It regressed to ZERO on 3 of 4 live phase-5
    papers (Keskin 3->0, Hu 5->0, Rojas 1->0) the moment the deterministic clonotype loaders removed the
    hand-building the flags used to ride along with: the gateway gate is satisfied by tcr_seq_methods OR
    data_depositions, so a record could finalize clean with 664 clonotypes and no flag at all.

    Fires when the record asserts real TCR-seq (status bulk/single_cell/both -- 'mentioned_only' is the
    honest no-detail state and is the GATEWAY's business, not this one) AND carries substantive TCR
    content (clonotypes or a seq method) AND the flag lane is empty. Returns a reason string, else ''."""
    status = getattr(rec, "tcr_seq_status", None)
    if status in (None, "none", "mentioned_only"):
        return ""                                   # no asserted TCR-seq -> nothing to flag
    if not (rec.tcr_clonotypes or rec.tcr_seq_methods):
        return ""                                   # gateway-level emptiness; _tcr_gateway_gap owns it
    if rec.neoantigen_tcr_flags:
        return ""
    return (
        f"tcr_seq_status={status!r} with {len(rec.tcr_clonotypes)} clonotype(s) and "
        f"{len(rec.tcr_seq_methods)} TCR-seq method(s) extracted, but neoantigen_tcr_flags is EMPTY -- "
        f"the FLAGSHIP per-(neoantigen x patient) 'was a TCR identified, CD4 or CD8?' fact is missing"
    )


def source_census(rec: "ExtractedPaper", out_path: str):
    """The source census for this record, or None when there is nothing trustworthy to census.

    Shared by the TCR and reactivity gates so a finalize pays for it ONCE -- a full census of Hu's
    257MB / 50-sheet supplement is ~2s, and both gates ask the same question of the same directory.
    Every failure path returns None, because a gate must block only on a sheet the census positively
    identified, never on one it could not read."""
    src = _read_source_dir(out_path, getattr(rec, "pmid", None))
    if not src:
        return None
    try:
        from . import tcr_source_census as _census        # lazy: keeps agent_core reader-dep-free
        return _census.census_tcr_sheets(src)
    except Exception:                                     # noqa: BLE001 - an unavailable census never blocks
        return None


# A reactivity sheet is under-loaded below this share of the positive calls the census counted.
_REACTIVITY_MIN_SHARE = 0.5


def _reactivity_source_coverage_gap(rec: "ExtractedPaper", out_path: str, census=None) -> str:
    """'' or a line naming reactivity matrices the source offers that this record barely loaded.

    The counterpart to `_read_under_coverage_gap`, and strictly stronger: that one keys off what the
    agent READ, so an agent that never issues a filtered read is never gated. This one asks the
    SOURCE. On Hu 33479501 the difference is concrete -- the census finds a third reactivity sheet
    ('Supp11b NetMHCpan2.4 CD4 T cell', 11 positive calls) that NONE of the five committed runs ever
    cited, because none of them read it.

    Counts positive assay CELLS, which is the grain the maintainer chose (2026-08-24): one evidence
    row per positive call. `add_table` with a `per_cell` mapping produces exactly that number."""
    census = census if census is not None else source_census(rec, out_path)
    if census is None or not getattr(census, "reactivity", None):
        return ""
    loaded: dict = {}
    for e in read_load_ledger(out_path):
        key = (str(e.get("file", "")).lower(), e.get("sheet"))
        loaded[key] = loaded.get(key, 0) + int(e.get("n_loaded") or 0)
    gaps = []
    short = [s for s in census.reactivity
             if loaded.get(s.key, 0) < _REACTIVITY_MIN_SHARE * s.n_expected]
    if not short:
        return ""
    # The census ALREADY knows which columns these are. Telling the agent beats making it infer them:
    # a per_cell mapping is authored fresh every run, so an inferred one is the step the run-to-run
    # swing moved to once the loader existed. Observed on Hu 33479501 with identical code -- one run
    # wrote 4b exactly (279 calls), the next wrote `equals: 0` rules on 4a as well (36 rows against 18
    # calls) and landed 4 short on 4b. Structure costs ~2s per workbook, so it is derived ONLY here,
    # on the path that is already blocking; the passing path never pays for it.
    detailed = census
    try:
        if any(s.imp_col is None for s in short):
            from . import tcr_source_census as _census
            src = _read_source_dir(out_path, getattr(rec, "pmid", None))
            if src:
                detailed = _census.census_tcr_sheets(src, with_structure=True)
    except Exception:                                  # noqa: BLE001 - a richer message never blocks
        detailed = census
    by_key = {s.key: s for s in getattr(detailed, "reactivity", [])}
    for s in short:
        d = by_key.get(s.key, s)
        gaps.append(
            f"{s.file} [{s.sheet}]: {s.n_expected} positive assay call(s), loaded "
            f"{loaded.get(s.key, 0)}. Its columns are already known -- header_row="
            f"{max(d.data_start - 1, 0)}, one per_cell rule per assay column {list(d.assay_cols)} "
            f"with \"equals\": 1, patient col {d.patient_col}, immunizing-peptide col {d.imp_col}"
            + (f", assay-peptide col {d.assay_pep_col}" if d.assay_pep_col is not None else "")
            + (f", id col {d.id_col} (recovers a blank patient cell)" if d.id_col is not None else ""))
    shown = gaps[:8]
    more = f" (+{len(gaps) - len(shown)} more)" if len(gaps) > len(shown) else ""
    return ("the source offers " + str(len(census.reactivity)) + " reactivity matrix/matrices this "
            "record under-loaded: " + "; ".join(shown) + more)


def _tcr_source_coverage_gap(rec: "ExtractedPaper", out_path: str, census=None) -> str:
    """UNDER-COVERAGE anchor (2026-08-24). Every other gate here asks "is this record internally
    consistent?"; this one asks "did it take what the SOURCE offered?" -- the question that let Hu
    33479501's keep-file finalize clean asserting tcr_seq_status='both' with 3 TCR-seq methods and
    ZERO clonotypes, beside a supplement holding 10 loadable tet-spec track sheets (269 clonotypes /
    654 timepoint observations). A 3-run replicate routed 11 / 6 / 1 sheets to the loaders on
    identical input and all three finalized clean: the loaders are deterministic GIVEN a sheet, so
    the non-determinism is entirely in WHICH sheets get routed, and only the source can see it.

    Fires when `tcr_source_census` finds a REQUIRED (loadable, in-scope) TCR sheet that the record's
    load ledger never recorded. Unlike its siblings this predicate needs `out_path` -- both sidecars
    and the source pointer hang off it.

    SILENT, by construction, whenever the census cannot be trusted rather than blocking on a guess:
    no source pointer (a record finalized outside the extraction runtime, e.g. `refinalize`, or one
    written before the ledger existed), a source dir that no longer exists, a pointer stamped with
    another paper's pmid, an unreadable workbook, or a census that simply found nothing loadable
    (Rojas 37165196 has no clonotype-shaped sheet at all -- its zero is HONEST and must not be
    flagged). Only a census that positively identified an unloaded required sheet returns a reason.

    Out-of-scope sheets never appear in `required` (see `tcr_source_census.OUT_OF_SCOPE_RULES`), so
    this gate can never coerce Keskin's non-antigen-sorted 'Tumor-associated' table into a record."""
    census = census if census is not None else source_census(rec, out_path)
    if census is None or not census.required:
        return ""
    # Only a TCR loader counts as covering a required TCR sheet. `add_table` also ledgers now, and
    # a generic bulk-load of a clonotype sheet must never be mistaken for the deterministic TCR
    # load this gate exists to demand -- that would turn the gate silently OFF, its worst failure.
    loaded = {(str(e.get("file", "")).lower(), e.get("sheet")) for e in read_load_ledger(out_path)
              if e.get("loader") in _TCR_LOADERS}
    missing = [s for s in census.required if s.key not in loaded]
    if not missing:
        return ""
    shown = missing[:12]
    lines = "; ".join(f"{s.file} [{s.sheet}] -> {s.loader}" for s in shown)
    more = f" (+{len(missing) - len(shown)} more)" if len(missing) > len(shown) else ""
    return (
        f"the source offers {len(census.required)} loadable TCR sheet(s) but this record loaded "
        f"{len(census.required) - len(missing)}; NEVER LOADED: {lines}{more}"
    )


def _read_under_coverage_gap(out_path: str) -> str:
    """'' or a line naming sheets the agent FILTERED-READ and then barely loaded.

    The Hu 33479501 re-extract read 'Suppl Dataset 4b' filtered to 121 rows -- twice, the second
    time with max_rows=121, so it had counted them -- then hand-wrote about three rows per patient
    and finalized clean at 53 evidence. Every existing gate passed: the patient-coverage gate
    ('covers 1 of 8 patients') is satisfied by ONE row per patient. Nothing compared what the agent
    was shown against what it kept.

    Deliberately narrow, so it can only fire on a real gap:
      - only FILTERED reads count (an unfiltered preview's total is the whole sheet);
      - only reads matching >= _READ_GATE_MIN_ROWS (a small filter is inspection, not a load);
      - a sheet is clear as soon as ANY deterministic loader took >= half of what was shown, and
        `add_table` with per_cell fan-out routinely exceeds it (121 rows -> 279 evidence).
    """
    if not _READ_OBSERVATIONS:
        return ""
    loaded: dict = {}
    for e in read_load_ledger(out_path):
        key = (str(e.get("file", "")).lower(), e.get("sheet"))
        loaded[key] = loaded.get(key, 0) + int(e.get("n_loaded") or 0)
    gaps = []
    for (fname, sheet), obs in sorted(_READ_OBSERVATIONS.items()):
        if obs["matched"] < _READ_GATE_MIN_ROWS:
            continue
        got = loaded.get((fname, sheet), 0)
        if got >= _READ_GATE_MIN_SHARE * obs["matched"]:
            continue
        gaps.append(f"{obs['file']} [{sheet}]: filter {obs['filter']} matched "
                    f"{obs['matched']}/{obs['total']} rows, deterministic loaders took {got}")
    if not gaps:
        return ""
    shown = gaps[:6]
    more = f" (+{len(gaps) - len(shown)} more)" if len(gaps) > len(shown) else ""
    return "; ".join(shown) + more


def _tcr_scope_note(rec: "ExtractedPaper", out_path: str) -> str:
    """One line naming the TCR sheets the census found but did NOT require, and why.

    POLICY VISIBILITY (docs/POLICY_tcr_sheet_scope_2026-08-24.md). `_tcr_source_coverage_gap` can only
    be made QUIETER by an out-of-scope classification, never louder -- so required->out_of_scope is the
    dangerous misclassification, and its failure mode is silent: coverage the gate simply stops asking
    about. An exclusion nobody can see is indistinguishable from a sheet the census failed to recognise.
    This surfaces every exclusion on the successful finalize so the policy stays auditable as the corpus
    grows. Purely informational -- it never blocks and never changes the record.

    Reuses the census the gate already ran when possible; silent whenever the census is unavailable
    (same conditions as the gate: no source pointer, pmid mismatch, unreadable workbook)."""
    src = _read_source_dir(out_path, getattr(rec, "pmid", None))
    if not src:
        return ""
    try:
        from . import tcr_source_census as _census
        census = _census.census_tcr_sheets(src)
    except Exception:                                     # noqa: BLE001 - never let a note break finalize
        return ""
    if not census.out_of_scope:
        return ""
    by_reason: dict = {}
    for s in census.out_of_scope:                         # `reason` is the rule's documented text
        short = (s.reason or "not antigen-sorted").split(",")[0].split(" (")[0].strip()[:70]
        by_reason.setdefault(short, []).append(s.sheet)
    parts = "; ".join(f"{reason} -> {', '.join(sheets[:3])}"
                      + (f" (+{len(sheets) - 3} more)" if len(sheets) > 3 else "")
                      for reason, sheets in sorted(by_reason.items()))
    return (f" | TCR scope: {len(census.out_of_scope)} sheet(s) present but NOT REQUIRED "
            f"[{parts}] (policy 2026-08-24: not antigen-sorted -> no neoantigen specificity; "
            f"loading them is neither required nor forbidden)")


def _target_id_index(imps: list, epis: list, cands: list = (), limit: int = 25) -> str:
    """'<GENE>:<paper_local_id>' index of the TCR-flaggable targets a record actually holds, in record
    order, truncated per kind. This is the ANSWER to the live runs' complaint that loader-minted ids are
    'opaque, unguessable': schema `_check_target` hard-rejects an unresolvable target, so the agent needs
    the real ids in hand before it can write a flag. Rows are plain dicts (partial) or model-like objects."""
    def _get(row, field):
        return row.get(field) if isinstance(row, dict) else getattr(row, field, None)

    parts = []
    for kind, rows in (("immunizing_peptide", imps), ("epitope", epis), ("candidate", cands)):
        rows = list(rows or [])
        if not rows:
            continue
        ids = [f"{_get(r, 'gene_symbol') or '?'}:{_get(r, 'paper_local_id')}" for r in rows]
        shown = ids[:limit]
        parts.append(f"{kind} ({len(ids)}) " + ", ".join(shown)
                     + (f", +{len(ids) - len(shown)} more" if len(ids) > len(shown) else ""))
    return " | ".join(parts)


def _tcr_flag_target_index(full: "ExtractedPaper", limit: int = 25) -> str:
    return _target_id_index(full.immunizing_peptides, full.epitopes, full.candidates, limit=limit)


# Deposition data_types that by themselves imply the paper HAS TCR-seq (so a lone deposition of one
# of these, with no method/clonotype/flag, still means tcr_seq_status should be non-'none').
_TCR_DEPOSITION_TYPES = frozenset({"raw_tcr_seq", "processed_clonotypes", "paired_scvdj"})


def _tcr_status_consistency_gap(rec: "ExtractedPaper") -> str:
    """v2.17 §5.3 consistency — the COMPLEMENT of the gateway anchor. The gateway fires when
    tcr_seq_status asserts TCR-seq but nothing was extracted; THIS fires when TCR content WAS
    extracted but the status sticker doesn't reflect it — either None/'none'/'mentioned_only' despite
    real methods/clonotypes/flags, or it understates the modalities the methods declare (e.g.
    status='single_cell' when a bulk method is also present -> should be 'both'). The live Phase-4 run
    under-set the sticker on Rojas (None) and Hu (single_cell despite a bulk method). Returns a reason
    string, else ''. Purely a label check — the extracted DATA is fine; overriding routes to needs_review."""
    status = getattr(rec, "tcr_seq_status", None)
    has_tcr = bool(
        rec.tcr_seq_methods or rec.tcr_clonotypes or rec.neoantigen_tcr_flags
        or any(d.data_type in _TCR_DEPOSITION_TYPES for d in rec.data_depositions)
    )
    if not has_tcr:
        return ""  # no TCR content -> nothing to be consistent with (gateway covers the inverse)
    if status in (None, "none"):
        return (f"TCR content was extracted (methods/clonotypes/flags) but tcr_seq_status={status!r}; "
                f"set it to bulk / single_cell / both to match")
    if status == "mentioned_only":
        return ("tcr_seq_status='mentioned_only' contradicts the TCR methods/clonotypes/flags actually "
                "extracted ('mentioned_only' means no detail was extractable); set bulk / single_cell / both")
    # modality consistency: methods declare which modalities ran; 'both' is an acceptable superset.
    mods = {m.modality for m in rec.tcr_seq_methods}
    expected = ("both" if {"bulk_tcr_seq", "single_cell_tcr_seq"} <= mods
                else "bulk" if "bulk_tcr_seq" in mods
                else "single_cell" if "single_cell_tcr_seq" in mods else None)
    if expected and status != "both" and status != expected:
        return (f"tcr_seq_status={status!r} disagrees with the extracted tcr_seq_methods modalities "
                f"{sorted(mods)} (expected {expected!r} or 'both')")
    return ""


def finalize_partial(out_path: str, allow_missing_magnitudes: bool = False,
                     allow_missing_pools: bool = False,
                     allow_member_level_pool_evidence: bool = False,
                     allow_candidate_bridge_mismatch: bool = False,
                     allow_unknown_funnel_size: bool = False,
                     allow_regimen_divergence: bool = False,
                     allow_evidence_count_mismatch: bool = False,
                     allow_peptide_count_mismatch: bool = False,
                     allow_sparse_evidence: bool = False,
                     allow_missing_class_ii: bool = False,
                     allow_missing_minimal_epitopes: bool = False,
                     allow_ungrounded_safety_grade: bool = False,
                     allow_missing_tcr_gateway: bool = False,
                     allow_tcr_status_mismatch: bool = False,
                     allow_missing_tcr_flags: bool = False,
                     allow_unloaded_tcr_sheets: bool = False,
                     allow_underloaded_read_sheets: bool = False,
                     allow_underloaded_reactivity_sheets: bool = False,
                     allow_cd4_on_class_i_minimal: bool = False,
                     allow_cd8_on_immunizing_peptide: bool = False,
                     allow_class_ii_asp_dump: bool = False) -> tuple[bool, str]:
    """Validate the assembled record against the FULL schema and write it.

    After the structural guards pass, SOFT checks block ONCE (each overridable):
    a magnitude check (immunogenic quantitative-assay responses with no magnitude); a
    pool-entity check (a patient has pooled evidence but no ExtractedPeptidePool entity);
    a pool-evidence-collapse check (pooled responses encoded as per-member rows instead of
    one pool row → run-to-run evidence-count variance); a candidate-bridge check (v2.8 P19:
    a selected candidate's sequence differs from its bridged IMP → label noise); and a
    funnel-size check (v2.8 P19: candidates exist but n_predicted_reported is unrecorded).
    """
    rec, err = _load_partial(out_path)
    if err:
        return False, err
    if rec is None:
        return False, "no partial record; call init_record first"
    # Deterministic epitope canonicalization: merge same-identity duplicate epitopes into the
    # schema's many-to-many form so the count is run-invariant (not the agent's varying dedup).
    n_epi_merged = canonicalize_epitopes(rec)
    # Orphaned-epitope repair (2026-08-16): derive the missing parent links by sequence containment,
    # then quarantine whatever still cannot be linked. Runs HERE (not only inside the manifest loader)
    # so peptides loaded AFTER the epitopes still rescue their orphans.
    n_epi_linked = link_orphan_epitopes(rec)
    quarantined, qerr = quarantine_orphan_epitopes(rec)
    if qerr:
        return False, qerr
    if quarantined:
        pathlib.Path(_quarantine_path(out_path)).write_text(json.dumps(
            {"reason": "epitope sequence not contained in any loaded immunizing peptide -- "
                       "no parent linkage the source supports; removed instead of fabricated",
             "removed": quarantined}, indent=2))
    # Evidence-lane determinism (stability lab 2026-07-06): canonicalize patient ids -> 'Pt<N>',
    # re-point immunogenic candidate/IMP evidence onto its gene's canonical immunizing peptide, and roll
    # up entailed per-pool rows. Together these take Keskin evidence stability 0.0 -> 1.0 (the runs record
    # the same biology; this removes the keying/grain noise). Order matters: ids first (so evidence + pool
    # refs are already canonical), then the target re-point, then the pool roll-up over canonical refs.
    # Aliases FIRST: a reference written to an id that was merged away at append time must resolve
    # before anything else reads the target graph (the re-point and roll-up passes below both do).
    n_alias = apply_id_aliases(rec, read_id_aliases(out_path))
    n_pt_canon = canonicalize_patient_ids(rec)
    n_ev_canon = canonicalize_evidence(rec)
    # Candidate CD8 rows must become IMP/epitope BEFORE H1 (Ott fix2 hung CD8 on candidates).
    n_cand_dropped = drop_administered_candidates(rec)
    # H1: CD8-on-IMP/candidate whose quote names exactly one nested 8-11mer → that epitope.
    n_cd8_repointed = _repoint_cd8_named_ept(rec)
    # O2: drop class-II 12-18mer ASP windows (inside an IMP, no affinity, no evidence).
    n_asp_dropped = _drop_class_ii_asp_windows(rec)
    n_pool_derived = derive_pool_evidence(rec)
    n_scr_dropped = _drop_screening_covered_by_evidence(rec)  # v2.14 #4: named evidence supersedes screening
    ok, msg = save_record(json.dumps(rec), out_path)  # full validation + write
    if not ok:
        return False, msg  # cross-entity errors; fix via clear/append then finalize again
    gok, gmsg = outer_guard(out_path)
    if not gok:
        return False, gmsg  # structural guard failed; keep the partial to fix
    full = ExtractedPaper(**json.loads(pathlib.Path(out_path).read_text()))
    # structural OK -> soft magnitude nudge (block-once, overridable)
    if not allow_missing_magnitudes:
        miss = _missing_magnitude_evidence(full)
        if miss:
            return False, (
                f"[finish] {len(miss)} immunogenic {sorted(_MAGNITUDE_ASSAYS)} responses have "
                f"NO magnitude (e.g. {miss[:3]}). Check the figure source-data for SFC/numeric "
                f"values and add them (value + sfc_per_1e6, or the per-patient set as raw + "
                f"needs_review). If a response truly has no reported magnitude, record a raw "
                f"saying so. To proceed anyway, call finalize with allow_missing_magnitudes=true."
            )
    # soft pool nudge (block-once, overridable): pooled evidence with no pool entity
    if not allow_missing_pools:
        need = _patients_needing_pool(full)
        if need:
            return False, (
                f"[finish] {len(need)} patient(s) {need} have pooled ELISpot evidence "
                f"(quoted_text names a neoantigen pool) but NO ExtractedPeptidePool entity. "
                f"Add one pool per patient (member_peptide_ids = that patient's in-pool "
                f"peptides; quoted_text from the pooled-response sentence/table value). If the "
                f"paper does not resolve pool membership, call finalize with allow_missing_pools=true."
            )
    # soft pool-evidence-collapse nudge (block-once, overridable): pooled responses encoded as
    # per-member rows instead of ONE pool row — the canonical rule that makes evidence count run-
    # invariant (root-caused on Rojas P25: "De novo response in pool" on 7 member rows vs gold's
    # 1 pool row + 3 deconvoluted member rows).
    if not allow_member_level_pool_evidence:
        split = _pool_evidence_not_collapsed(full)
        if split:
            return False, (
                f"[finish] {len(split)} patient(s) {split} have a POOLED response encoded as "
                f"multiple per-member evidence rows (quoted_text says 'pool') with NO consolidating "
                f"pool-target row. Canonical rule: a response measured at the POOL level is ONE "
                f"evidence row with target_kind='pool' (pool_paper_id = that patient's pool); keep a "
                f"per-peptide row ONLY for a member the source reports INDIVIDUALLY (deconvoluted, no "
                f"'in pool' wording). Replace the 'in pool' member rows with one pool row. If the "
                f"source truly deconvolutes every member, call finalize with "
                f"allow_member_level_pool_evidence=true."
            )
    # soft candidate-bridge nudge (v2.8 P19, block-once, overridable): a selected
    # candidate's sequence disagrees with the IMP it bridges to → label noise.
    if not allow_candidate_bridge_mismatch:
        mism = _candidate_bridge_seq_mismatch(full)
        if mism:
            return False, (
                f"[finish] {len(mism)} candidate(s) {mism} have a selected_peptide_id whose "
                f"immunizing peptide has a DIFFERENT sequence than the candidate — a wrong bridge "
                f"would attach the candidate's prioritization scores to the wrong outcome. Verify "
                f"each candidate's sequence matches its bridged IMP (fix the sequence or the "
                f"selected_peptide_id), or set needs_review on the affected candidate(s). To "
                f"proceed anyway, call finalize with allow_candidate_bridge_mismatch=true."
            )
    # soft funnel-completeness nudge (v2.8 P19, block-once, overridable): candidates
    # exist but the paper-stated predicted count (the denominator) is unrecorded.
    if not allow_unknown_funnel_size:
        if _funnel_size_unknown(full):
            return False, (
                f"[finish] {len(full.candidates)} candidate(s) recorded but n_predicted_reported "
                f"is unset, so a truncated funnel (e.g. 50 of 322 predicted) can't be told from a "
                f"complete one and any selection/immunogenic-per-predicted rate has a silently wrong "
                f"denominator. Set n_predicted_reported (and n_selected_reported) to the count(s) the "
                f"paper STATES. If the paper does not report a predicted total, call finalize with "
                f"allow_unknown_funnel_size=true."
            )
    # soft regimen-consistency nudge (v2.9 P21, block-once, overridable): patients of one arm
    # carry divergent trial-constant delivery regimen → likely an inconsistent per-patient re-derivation.
    if not allow_regimen_divergence:
        div = _regimen_divergence(full)
        if div:
            return False, (
                f"[finish] patients sharing vaccine_platform {div} have DIVERGENT trial-constant "
                f"delivery regimen (adjuvant / dose / schedule). The regimen is normally identical "
                f"for every patient of an arm — extract it ONCE and apply it to all of them "
                f"(per-patient fields like surgery→dose latency and doses-received MAY differ). If "
                f"the arm genuinely uses different regimens (e.g. dose-escalation cohorts), call "
                f"finalize with allow_regimen_divergence=true."
            )
    # soft evidence-completeness nudge (v2.11.1 P20.1, block-once, overridable): recorded evidence
    # materially disagrees with the paper's STATED counts → the granularity/recall variance the
    # diagnostic root-caused. Fires only when an anchor was set.
    if not allow_evidence_count_mismatch:
        gaps = _evidence_anchor_gap(full)
        if gaps:
            return False, (
                "[finish] evidence vs the paper's stated counts: " + "; ".join(gaps) +
                ". Fix it (add the missed immunogenic rows / match the named negative grain) or, if the "
                "paper's count is at a different grain (e.g. pooled), call finalize with "
                "allow_evidence_count_mismatch=true."
            )
    # soft peptide-recall nudge (prototype, block-once, overridable): recorded peptides fall short of
    # the paper's STATED selected/administered count — closes the gameable strict reconciliation
    # (lowering n_peptides_synthesized to match a missed table). Fires only when n_selected_reported set.
    if not allow_peptide_count_mismatch:
        pgap = _peptide_recall_gap(full)
        if pgap:
            return False, (
                "[finish] peptide recall vs the paper's stated count: " + pgap +
                " To proceed anyway (e.g. the paper's count includes peptides not individually listed), "
                "call finalize with allow_peptide_count_mismatch=true."
            )
    # soft evidence-breadth nudge (prototype, block-once, overridable): responses cover few of the
    # vaccinated patients → likely an incomplete per-patient sweep (e.g. many per-patient sheets).
    if not allow_sparse_evidence:
        bgap = _evidence_breadth_gap(full)
        if bgap:
            return False, (
                "[finish] " + bgap + " To proceed anyway, call finalize with allow_sparse_evidence=true."
            )
    # soft class-II minting nudge (v2.12 P23; cue-split v2.19): a quoted class-II restriction
    # (DR/DP/DQ / I-A/I-E / 'class II-restricted') is present but no class-II record was minted.
    # CD4/helper alone does not fire this gate.
    if not allow_missing_class_ii:
        c2 = _class_ii_minting_gap(full)
        if c2:
            return False, (
                "[finish] class-II coverage: " + "; ".join(c2) +
                ". Mint the class-II epitope(s)/type the evidence, or if the paper truly reports only "
                "class-I responses, call finalize with allow_missing_class_ii=true."
            )
    # soft minimal-epitope GRAIN nudge (v2.12 P23 follow-up, block-once, overridable): a class-I epitope
    # is a long peptide mislabeled, or the predicted minimal-epitope layer was dropped (the lite-Li miss).
    if not allow_missing_minimal_epitopes:
        meg = _minimal_epitope_grain_gap(full)
        if meg:
            return False, (
                "[finish] minimal-epitope grain: " + "; ".join(meg) +
                ". Fix the epitope grain (minimal binder in epitopes, long peptide in immunizing_peptides), "
                "or if the paper truly reports no minimal epitopes, call finalize with "
                "allow_missing_minimal_epitopes=true."
            )
    # soft safety-grounding nudge (v2.15, block-once, overridable): a grade>=3 treatment-related claim
    # whose own `raw` quote does not support it (SAE/'serious' wording, disease-progression disclaimer,
    # or no grade at all) -- the seriousness/grade/relatedness conflation root-caused on 39041242.
    if not allow_ungrounded_safety_grade:
        sg = _safety_grade_ungrounded(full)
        if sg:
            return False, (
                "[finish] " + sg + ". 'Serious'/SAE is a regulatory category, NOT a CTCAE grade; only a "
                "grade>=3 (or 'severe'/'life-threatening') AE attributed to the STUDY TREATMENT sets "
                "any_grade3plus_related=true. Re-call set_safety_summary with the correct grade (and a "
                "`raw` that states it), or if the paper truly reports a grade>=3 related AE whose grade "
                "you cannot quote verbatim, call finalize with allow_ungrounded_safety_grade=true."
            )
    # SOURCE UNDER-COVERAGE gate (2026-08-24, block-once, overridable) -- runs FIRST among the TCR
    # gates because loading a missed sheet changes the content the others judge (status sticker,
    # flags), so fixing coverage first is the sequence that converges. Computed ONCE and reused in
    # the override-audit tuple below: unlike its record-only siblings this predicate does source I/O
    # (a full census of Hu's 257MB / 50-sheet supplement is ~0.6s), so re-deriving it there would
    # double that for no new information.
    # READ UNDER-COVERAGE gate (2026-08-24): the agent was SHOWN N filtered rows of a sheet and kept
    # a fraction of them. Runs before the TCR gates for the same reason they run early -- loading a
    # sheet properly changes the content every later gate judges.
    read_gap = _read_under_coverage_gap(out_path)
    if not allow_underloaded_read_sheets and read_gap:
        return False, (
            "[finish] source under-coverage: " + read_gap + ". You read these sheets with a filter, "
            "so the row count you saw is the count the sheet offers -- keeping a fraction of it is "
            "recall loss, not selection. Load them with add_table instead of transcribing by hand; "
            "for a reactivity/immunogenicity sheet whose assay columns sit side by side on one row, "
            "use the per_cell mapping so EVERY positive cell becomes its own evidence row (a plain "
            "row mapping reaches only one assay column). If the filter was genuinely just an "
            "inspection, or the sheet's rows legitimately collapse to far fewer records, call "
            "finalize with allow_underloaded_read_sheets=true (the record then routes to needs_review)."
        )
    census = source_census(full, out_path)          # ~2s on Hu's supplement -- paid once, shared
    reactivity_gap = _reactivity_source_coverage_gap(full, out_path, census)
    if not allow_underloaded_reactivity_sheets and reactivity_gap:
        return False, (
            "[finish] evidence source coverage: " + reactivity_gap + ". A reactivity matrix scores "
            "several assays per row (timepoints x ex-vivo / pre-stimulated / minigene / autologous "
            "tumour) and EVERY positive cell is one evidence row. Load each sheet with add_table "
            "using a `per_cell` mapping — a plain row mapping reaches only one assay column and "
            "silently drops the rest — giving each cell rule its own timepoint_phase, "
            "assay_stimulation and assay_antigen_format. If a sheet's calls are genuinely already "
            "recorded from another table, or it reports something this paper does not attribute to "
            "a patient, call finalize with allow_underloaded_reactivity_sheets=true (the record "
            "then routes to needs_review)."
        )
    coverage_gap = _tcr_source_coverage_gap(full, out_path, census)
    if not allow_unloaded_tcr_sheets and coverage_gap:
        return False, (
            "[finish] TCR source coverage: " + coverage_gap + ". These sheets ARE loadable by a "
            "deterministic loader (the census probed their headers, not their names), and a loader "
            "reads the WHOLE sheet -- including every timepoint column -- so this is recall you "
            "cannot match by hand. Call the named tool once per sheet (add_tcr_track: pass `path`, "
            "`sheet`, `patient`; add_reactive_tcr: also `antigen`), then finalize again. If a sheet "
            "genuinely belongs to a patient/antigen this paper does not report, or the loader comes "
            "back empty on it, call finalize with allow_unloaded_tcr_sheets=true (the record then "
            "routes to needs_review)."
        )
    # soft TCR-gateway nudge (v2.17 §5.7, block-once, overridable): the paper asserts TCR-seq was done
    # (tcr_seq_status non-'none') but no method AND no deposition were extracted — the Tier-A gateway fact
    # is known-missing. HARD override (a real recall gap): overriding routes the record to needs_review.
    if not allow_missing_tcr_gateway:
        tg = _tcr_gateway_gap(full)
        if tg:
            return False, (
                "[finish] TCR gateway: " + tg + ". Add the TCR-seq method (platform / software / "
                "clonotype_definition, from Methods) and/or the data-deposition accession (from Data "
                "availability). If the paper only MENTIONS TCR-seq with no extractable recipe or "
                "accession, set tcr_seq_status='mentioned_only' and call finalize with "
                "allow_missing_tcr_gateway=true (the record then routes to needs_review)."
            )
    # soft TCR-status-consistency nudge (v2.17 §5.3, block-once, overridable): TCR content was extracted
    # but tcr_seq_status is missing/understated (the Phase-4 Rojas/Hu miss). Trivially fixable (set the
    # one field); HARD override so an insisted-on mismatch routes to needs_review.
    if not allow_tcr_status_mismatch:
        ts = _tcr_status_consistency_gap(full)
        if ts:
            return False, (
                "[finish] TCR status: " + ts + ". Set tcr_seq_status on the paper to match the TCR "
                "content you extracted (use 'both' when you recorded both a bulk and a single-cell "
                "method). To proceed anyway, call finalize with allow_tcr_status_mismatch=true."
            )
    # FLAGSHIP TCR-flag nudge (2026-08-16, block-once, overridable), SYMMETRIC with the gateway above:
    # real TCR-seq content but an EMPTY neoantigen_tcr_flags lane. HARD override -> needs_review. The
    # message hands over the ACTUAL target ids, because the deterministic loaders mint opaque ids
    # (IMP_<12aa>_<len>_<sha1_6>) and schema._check_target HARD-REJECTS an unresolvable one, so a blind
    # hand-write cannot pass validation -- the three failing live runs all called the ids "unguessable".
    if not allow_missing_tcr_flags:
        fg = _tcr_flag_gap(full)
        if fg:
            idx = _tcr_flag_target_index(full)
            return False, (
                "[finish] TCR flags: " + fg + ". Add one neoantigen_tcr_flags row per (neoantigen x "
                "patient) that has TCR evidence: target_kind 'immunizing_peptide'|'epitope'|'candidate' "
                "+ that entity's OWN paper_local_id (exactly one of immunizing_peptide_paper_id / "
                "epitope_paper_id / candidate_paper_id), patient_paper_id, tcr_identified in "
                "cd4/cd8/both/none/not_assessed, plus quoted_text + section_ref naming the figure/table "
                "you read it from. The target id MUST already exist in this record -- a guessed id is "
                "REJECTED by validation -- so call `partial_status` FIRST and copy the real ids it lists "
                "(the loaders mint ids like IMP_<12aa>_<len>_<sha1_6>; never invent one). "
                + (f"Targets here: {idx}. " if idx else "")
                + "If the paper truly makes no per-target TCR call, call finalize with "
                "allow_missing_tcr_flags=true (the record then routes to needs_review)."
            )
    # CD4 vs class-I 9-mer grain (v2.19 follow-up): block CD4 evidence hung on a nested class-I
    # 8-11mer (the Keskin EPI_RESP_* shape). Point at the parent IMP instead.
    if not allow_cd4_on_class_i_minimal:
        cg = _cd4_class_i_minimal_target_gap(full)
        if cg:
            return False, (
                "[finish] CD4 evidence grain: " + "; ".join(cg) +
                " To proceed anyway because the paper quotes that exact 9-mer as the CD4 target, "
                "call finalize with allow_cd4_on_class_i_minimal=true (routes to needs_review)."
            )
    # CD8 vs long-IMP grain: block CD8 hung on a long immunizing peptide only when the
    # row's own quote NAMES a nested class-I 8-11mer (EPT12A / sequence). Bare "(EPT)"
    # with several Table-S5 bites (SLX4MUT) stays on the IMP (decision A).
    if not allow_cd8_on_immunizing_peptide:
        sg = _cd8_long_peptide_ept_gap(full)
        if sg:
            return False, (
                "[finish] CD8 evidence grain: " + "; ".join(sg) +
                " To proceed anyway because the paper assayed the LONG peptide for CD8, "
                "call finalize with allow_cd8_on_immunizing_peptide=true (routes to needs_review)."
            )
    # Overlapping ASP 15-mers minted as class-II epitopes (Ott SI5 dump). Named-DR Rojas
    # epitopes are ~1 per IMP and do not fire.
    if not allow_class_ii_asp_dump:
        ag = _class_ii_asp_dump_gap(full)
        if ag:
            return False, (
                "[finish] class-II ASP dump: " + "; ".join(ag) +
                " To proceed anyway, call finalize with allow_class_ii_asp_dump=true "
                "(routes to needs_review)."
            )
    # OVERRIDE AUDIT TRAIL (v2.11.3): reaching here means every soft guard either did not fire or was
    # OVERRIDDEN. Re-derive which fired AND were allowed (= overridden) and persist them on the record,
    # so a record that finalized DESPITE a fired guard is not silently clean — the scale lane routes a
    # non-empty `finalize_overrides_used` to needs_review (the live test showed agents override expensive
    # recall nudges rather than recover; this makes the override an auditable QC signal).
    overrides_used = [name for name, allowed, fired in (
        ("allow_missing_magnitudes", allow_missing_magnitudes, bool(_missing_magnitude_evidence(full))),
        ("allow_missing_pools", allow_missing_pools, bool(_patients_needing_pool(full))),
        ("allow_member_level_pool_evidence", allow_member_level_pool_evidence,
         bool(_pool_evidence_not_collapsed(full))),
        ("allow_candidate_bridge_mismatch", allow_candidate_bridge_mismatch,
         bool(_candidate_bridge_seq_mismatch(full))),
        ("allow_unknown_funnel_size", allow_unknown_funnel_size, _funnel_size_unknown(full)),
        ("allow_regimen_divergence", allow_regimen_divergence, bool(_regimen_divergence(full))),
        ("allow_evidence_count_mismatch", allow_evidence_count_mismatch, bool(_evidence_anchor_gap(full))),
        ("allow_peptide_count_mismatch", allow_peptide_count_mismatch, bool(_peptide_recall_gap(full))),
        ("allow_sparse_evidence", allow_sparse_evidence, bool(_evidence_breadth_gap(full))),
        ("allow_missing_class_ii", allow_missing_class_ii, bool(_class_ii_minting_gap(full))),
        ("allow_missing_minimal_epitopes", allow_missing_minimal_epitopes,
         bool(_minimal_epitope_grain_gap(full))),
        ("allow_ungrounded_safety_grade", allow_ungrounded_safety_grade,
         bool(_safety_grade_ungrounded(full))),
        ("allow_missing_tcr_gateway", allow_missing_tcr_gateway, bool(_tcr_gateway_gap(full))),
        ("allow_tcr_status_mismatch", allow_tcr_status_mismatch, bool(_tcr_status_consistency_gap(full))),
        ("allow_missing_tcr_flags", allow_missing_tcr_flags, bool(_tcr_flag_gap(full))),
        ("allow_unloaded_tcr_sheets", allow_unloaded_tcr_sheets, bool(coverage_gap)),
        ("allow_underloaded_read_sheets", allow_underloaded_read_sheets, bool(read_gap)),
        ("allow_underloaded_reactivity_sheets", allow_underloaded_reactivity_sheets,
         bool(reactivity_gap)),
        ("allow_cd4_on_class_i_minimal", allow_cd4_on_class_i_minimal,
         bool(_cd4_class_i_minimal_target_gap(full))),
        ("allow_cd8_on_immunizing_peptide", allow_cd8_on_immunizing_peptide,
         bool(_cd8_long_peptide_ept_gap(full))),
        ("allow_class_ii_asp_dump", allow_class_ii_asp_dump,
         bool(_class_ii_asp_dump_gap(full))),
    ) if allowed and fired]
    # Always rewrite the audit trail, including []. A seed/partial that still carries a
    # previous HARD stamp must not keep it after the gap no longer fires (replay / refinalize).
    full = full.model_copy(update={"finalize_overrides_used": overrides_used})
    pathlib.Path(out_path).write_text(full.model_dump_json(indent=2, exclude_none=True))
    pathlib.Path(_partial_path(out_path)).unlink(missing_ok=True)  # success: clean up partial
    pathlib.Path(_lock_path(out_path)).unlink(missing_ok=True)     # and its lane locks
    if n_epi_merged:
        gmsg += f" | canonicalized epitopes: merged {n_epi_merged} duplicate record(s)"
    if n_epi_linked:
        gmsg += (f" | epitope linkage: derived parent_peptide_ids for {n_epi_linked} orphaned "
                 f"epitope(s) by sequence containment")
    if quarantined:
        gmsg += (f" | QUARANTINED {len(quarantined.get('epitopes', []))} unlinkable epitope(s) "
                 f"(+{sum(len(v) for k, v in quarantined.items() if k != 'epitopes')} referencing row(s)) "
                 f"-> {_quarantine_path(out_path)} -- needs_review")
    if n_scr_dropped:
        gmsg += f" | screening: dropped {n_scr_dropped} row(s) superseded by named evidence"
    if n_asp_dropped:
        gmsg += (f" | dropped {n_asp_dropped} ASP-window class-II epitope(s) "
                 f"(12-18mer inside an IMP, no affinity, no evidence)")
    if overrides_used:
        gmsg += f" | OVERRIDES USED -> needs_review: {overrides_used}"
    gmsg += _tcr_scope_note(full, out_path)
    return True, gmsg


def refinalize(in_path: str, out_path: str, **allow_kwargs) -> tuple[bool, str]:
    """Replay finalize's deterministic transforms on a finished record into a NEW path.

    Does not overwrite `in_path`. Seeds `out_path.partial.json` from the finished JSON
    (drops `finalize_overrides_used` so the audit trail is recomputed), then
    `finalize_partial`. Default: no overrides — if a HARD gap still fires, this fails
    rather than re-stamping. Pass allow_* kwargs only for gaps that should still be
    overridable on this replay.
    """
    src = pathlib.Path(in_path)
    dest = pathlib.Path(out_path)
    if not src.is_file():
        return False, f"no record at {in_path}"
    try:
        if src.resolve() == dest.resolve():
            return False, "refinalize refuses to overwrite the source; pass a new out_path"
    except OSError:
        pass
    try:
        rec = json.loads(src.read_text())
    except Exception as e:
        return False, f"could not read {in_path}: {e}"
    if not isinstance(rec, dict):
        return False, "source is not a JSON object"
    rec.pop("finalize_overrides_used", None)
    dest.parent.mkdir(parents=True, exist_ok=True)
    pathlib.Path(_lock_path(str(dest))).unlink(missing_ok=True)
    pathlib.Path(_quarantine_path(str(dest))).unlink(missing_ok=True)
    # A replay carries no load ledger and no source pointer of its OWN; anything left at `dest` by a
    # previous record would be read as this record's coverage. Clear both -> the coverage census is
    # unavailable here and stays SILENT, which is the designed degradation for an out-of-runtime
    # finalize (it must never block a replay, and must never bless one either).
    pathlib.Path(_ledger_path(str(dest))).unlink(missing_ok=True)
    pathlib.Path(_source_path(str(dest))).unlink(missing_ok=True)
    pathlib.Path(_partial_path(str(dest))).write_text(json.dumps(rec))
    try:
        return finalize_partial(str(dest), **allow_kwargs)
    except TypeError as e:
        pathlib.Path(_partial_path(str(dest))).unlink(missing_ok=True)
        return False, str(e)


def _read_xlsx_dicts(path: str, sheet: str | None = None,
                     header_row: int | None = None) -> tuple[list, list]:
    """Read an .xlsx into (headers, rows) where each row is a dict keyed BOTH by header
    name AND by 0-based integer position. The positional keys let the mapping DSL address
    columns that have no usable name — merged / blank / duplicated headers, ubiquitous in
    supplementary tables (Keskin S5's peptide-ID and affinity columns are nameless).
    Leading TITLE rows are auto-skipped (or set header_row) and trailing all-empty (phantom)
    columns trimmed. Cell VALUES are kept full (no clipping) — this is the write path, not a
    preview. Fully-empty rows are skipped."""
    openpyxl = _require("openpyxl")

    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb[sheet] if sheet else wb[wb.sheetnames[0]]
    matrix = _trim_width([list(r) for r in ws.iter_rows(values_only=True)])
    if not matrix:
        raise ValueError(f"sheet {ws.title!r} has no rows (not even a header)")
    hi = header_row if header_row is not None else _header_index(matrix)
    headers = [str(h) if h is not None else "" for h in matrix[hi]]
    rows = []
    for r in matrix[hi + 1:]:
        if all(c is None for c in r):
            continue
        d = {i: r[i] for i in range(len(r))}                                   # positional keys
        d.update({headers[i]: r[i] for i in range(min(len(headers), len(r)))})  # name keys (last-wins, as before)
        rows.append(d)
    return headers, rows


def _norm_seq(s) -> str:
    """Normalize a peptide sequence for matching: drop underline markup + whitespace, uppercase."""
    if s is None:
        return ""
    return "".join(str(s).replace("<u>", "").replace("</u>", "").split()).upper()


def build_patient_pools(out_path: str, path: str, patient_col, peptide_col,
                        sheet: str | None = None, section_ref: str = "",
                        quoted_text_template: str | None = None,
                        header_row: int | None = None,
                        pool_id_prefix: str = "POOL-") -> tuple[bool, str]:
    """Deterministically build ONE ExtractedPeptidePool per patient from a per-patient
    peptide-ASSIGNMENT sheet (one row per patient x peptide — e.g. 33064988 'Vaccine peptides':
    'Patient Alias' | 'Peptide Sequence'). Groups the sheet by `patient_col`, matches each
    `peptide_col` sequence to an already-loaded `immunizing_peptides` row, and appends one pool
    per patient with `member_peptide_ids` = that patient's matched peptide ids.

    WHY: hand-building these pools via add_entities is per-patient reasoning that swung 0<->34
    pools run-to-run (live 2026-06-09). The group-by here is a pure function of the sheet, so the
    pool set is IDENTICAL every run. Sequence match is patient-disambiguated: when a sequence maps
    to several loaded peptides (one per patient, ids like 'IMP-<patient>-<seq>'), the candidate
    whose id carries the patient token wins; a lone candidate is used as-is (many-to-many is legal).
    """
    rec, err = _load_partial(out_path)
    if err:
        return False, err
    if rec is None:
        return False, "no partial record; call init_record first"
    peps = rec.get("immunizing_peptides") or []
    if not peps:
        return False, "no immunizing_peptides loaded yet; add the peptides (add_table) BEFORE build_patient_pools"
    seqmap: dict[str, list[str]] = {}
    for p in peps:
        s = _norm_seq(p.get("sequence"))
        if s:
            seqmap.setdefault(s, []).append(p.get("paper_local_id"))

    try:
        headers, rows = _read_xlsx_dicts(path, sheet, header_row)
    except Exception as e:
        return False, f"could not read sheet: {e}"

    def _resolve(col):
        if isinstance(col, int):
            return col if 0 <= col < len(headers) else None
        return col if col in headers else None
    pcol, scol = _resolve(patient_col), _resolve(peptide_col)
    if pcol is None:
        return False, f"patient_col {patient_col!r} not found; headers: {headers}"
    if scol is None:
        return False, f"peptide_col {peptide_col!r} not found; headers: {headers}"

    groups: dict[str, list[str]] = {}          # patient -> ordered-unique member ids
    n_rows = unmatched = 0
    for r in rows:
        patient = r.get(pcol)
        patient = str(patient).strip() if patient not in (None, "") else ""
        seq = _norm_seq(r.get(scol))
        if not patient or not seq:
            continue
        n_rows += 1
        cands = seqmap.get(seq) or []
        if not cands:
            unmatched += 1
            continue
        if len(cands) == 1:
            pid = cands[0]
        else:                                  # disambiguate by the patient token in the id
            tok = re.compile(rf"(?<![A-Za-z0-9]){re.escape(patient)}(?![A-Za-z0-9])")
            hit = [c for c in cands if c and tok.search(c)]
            pid = hit[0] if hit else cands[0]
        members = groups.setdefault(patient, [])
        if pid not in members:
            members.append(pid)

    pools = []
    for patient, members in groups.items():
        n = len(members)
        qt = (quoted_text_template.format(patient=patient, n=n) if quoted_text_template
              else f"Patient {patient} immunizing pool: {n} peptide(s), grouped from the "
                   f"per-patient peptide-assignment table.")
        pools.append({
            "paper_local_id": f"{pool_id_prefix}{patient}",
            "patient_paper_id": patient,
            "member_peptide_ids": members,
            "quoted_text": qt,
            "section_ref": section_ref or (sheet or "per-patient peptide-assignment table"),
        })
    if not pools:
        return False, (f"0 pools built ({n_rows} assignment rows, {unmatched} unmatched to loaded "
                       "peptides). Check patient_col/peptide_col and that peptides are loaded.")
    ok, msg = _validate_and_append(out_path, "pools", pools)
    skipped = ""  # patients are only created when they have >=1 matched member, so none are silently lost
    return ok, (f"{msg} | build_patient_pools: {len(pools)} pools over {len(pools)} patients, "
                f"{sum(len(m) for m in groups.values())} members from {n_rows} rows"
                + (f"; {unmatched} sequences unmatched to loaded peptides" if unmatched else ""))


def build_pool_evidence(out_path: str, patients: list | None = None,
                        sheets_path: str | None = None, sheet_pattern: str | None = None,
                        assay: str = "elispot", outcome: str = "immunogenic",
                        section_ref: str = "Figure 4A; Supplemental Figure 3",
                        provenance_locator: str = "Figure 4A; Supplemental Figure 3",
                        quoted_text_template: str | None = None,
                        assay_detail: str | None = None) -> tuple[bool, str]:
    """Deterministically emit ONE pool-target evidence row per MONITORED patient (A' design,
    2026-06-09). For papers whose per-patient pool immunogenicity is stated uniformly in text
    ("responses in all patients") but quantified only in per-patient FIGURES (33064988 Fig 4A /
    Supp Fig 3 — no backing table), the agent built these rows by hand via add_entities and it
    swung 0<->34 run-to-run. This generates them deterministically.

    PROVENANCE/FAITHFULNESS (per the A' review): the monitored set is the patients the paper
    INDIVIDUALLY shows (here the per-patient `IAP-<patient>` tabs / the pools build_pools made),
    NOT every vaccinated patient — the aggregate "all patients" sentence must not manufacture
    per-patient rows the figures never display. Each row is anchored to the per-patient FIGURE
    (provenance kind='figure', so a consumer never confuses it with a readable-table row),
    carries magnitude=None + needs_review=True (a clean target for an out-of-band figure-vision
    backfill), and references the patient's existing pool (looked up in the record, not assumed).

    Monitored set = explicit `patients`, OR derived from a workbook's sheet names matching
    `sheet_pattern` (e.g. 'IAP-(.+)' over mmc2.xlsx -> B1, L2, M13, ...).
    """
    rec, err = _load_partial(out_path)
    if err:
        return False, err
    if rec is None:
        return False, "no partial record; call init_record first"
    pool_by_patient = {p.get("patient_paper_id"): p.get("paper_local_id")
                       for p in (rec.get("pools") or []) if p.get("patient_paper_id")}
    if not pool_by_patient:
        return False, "no pools in the record; call build_pools BEFORE build_pool_evidence"

    if patients:
        monitored = [str(p).strip() for p in patients if str(p).strip()]
    elif sheets_path and sheet_pattern:
        openpyxl = _require("openpyxl")
        try:
            names = openpyxl.load_workbook(sheets_path, read_only=True).sheetnames
        except Exception as e:
            return False, f"could not open {sheets_path}: {e}"
        rx = re.compile(sheet_pattern)
        monitored = []
        for nm in names:
            m = rx.search(nm)
            if m:
                monitored.append((m.group(1) if m.groups() else m.group(0)).strip())
    else:
        return False, "give either patients=[...] or sheets_path+sheet_pattern (e.g. 'IAP-(.+)')"

    seen, rows, skipped = set(), [], []
    for patient in monitored:
        if patient in seen:
            continue
        seen.add(patient)
        pool_id = pool_by_patient.get(patient)
        if not pool_id:
            skipped.append(patient)               # monitored but no pool entity -> do NOT invent one
            continue
        qt = (quoted_text_template.format(patient=patient) if quoted_text_template
              else (f"Patient {patient}: de novo neoantigen-specific T cell response detected "
                    f"post-vaccination by IFN-gamma ELISpot (per-patient figure)."))
        rows.append({
            "patient_paper_id": patient,
            "target_kind": "pool",
            "pool_paper_id": pool_id,
            "assay": assay,
            "outcome": outcome,
            "vaccine_induced": True,
            "magnitude": None,
            "needs_review": True,
            "confidence": 2,
            "assay_detail": assay_detail or ("per-patient de novo ELISpot response; magnitude not "
                                             "individually tabulated (figure-derived) -- needs_review"),
            "quoted_text": qt,
            "section_ref": section_ref,
            "provenance": [{"kind": "figure", "locator": provenance_locator}],
        })
    if not rows:
        return False, (f"0 evidence rows built ({len(seen)} monitored, {len(skipped)} had no pool). "
                       "Check the monitored set / that build_pools ran first.")
    ok, msg = _validate_and_append(out_path, "evidence", rows)
    note = f"; {len(skipped)} monitored patients had no pool (skipped: {skipped[:8]})" if skipped else ""
    return ok, (f"{msg} | build_pool_evidence: {len(rows)} pool/immunogenic rows over {len(rows)} "
                f"monitored patients (magnitude=null, needs_review, figure-provenance){note}")


def build_reactivity_evidence(out_path: str, paper_dir: str | None = None,
                              sheet: str | None = None) -> tuple[bool, str]:
    """Deterministically load every (or one named) reactivity matrix the census finds.

    One evidence row per positive assay CELL. The mapping is authored here, not by
    the agent — Hu 33479501's 4b/11b swing was the agent authoring a fresh per_cell
    list each run. See vaxtract.reactivity_adapter.
    """
    rec, err = _load_partial(out_path)
    if err:
        return False, err
    if rec is None:
        return False, "no partial record; call init_record first"
    src = paper_dir or _read_source_dir(out_path, rec.get("pmid"))
    if not src:
        return False, ("no paper directory: pass paper_dir or call this inside a live extract "
                       "so the source sidecar is set")
    from . import reactivity_adapter as _rx
    items, minted, report = _rx.build_items(rec, src, sheet=sheet)
    if report.get("error"):
        return False, report["error"]
    if minted:
        ok, msg = _validate_and_append(out_path, "epitopes", minted)
        if not ok:
            return False, f"minted epitopes rejected: {msg}"
    if not items:
        return False, (
            f"0 evidence rows from {len(report.get('sheets') or [])} reactivity sheet(s) "
            f"(cells={report.get('n_cells')}, skipped_no_patient={report.get('skipped_no_patient')}, "
            f"skipped_no_target={report.get('skipped_no_target')}, skipped_dup={report.get('skipped_dup')}). "
            "Load patients and immunizing peptides first."
        )
    ok, msg = _validate_and_append(out_path, "evidence", items)
    if not ok:
        return False, msg
    for s in report.get("sheets") or []:
        # n_available is the census positive-cell count (the grain the gate measures),
        # not the sheet's raw row count.
        _append_load_ledger(out_path, path=s.get("path") or s.get("file") or src,
                            sheet=s["sheet"], loader="build_reactivity_evidence",
                            n_loaded=s["n_loaded"], n_available=s["n_expected"])
    bits = [f"{s['sheet']}:{s['n_loaded']}/{s['n_expected']}" for s in report.get("sheets") or []]
    extra = []
    if report.get("n_minted_epitopes"):
        extra.append(f"minted {report['n_minted_epitopes']} assay-peptide epitope(s)")
    if report.get("skipped_no_patient"):
        extra.append(f"{report['skipped_no_patient']} cell(s) skipped (patient not in record)")
    if report.get("skipped_no_target"):
        extra.append(f"{report['skipped_no_target']} cell(s) skipped (no IMP/epitope match)")
    if report.get("skipped_dup"):
        extra.append(f"{report['skipped_dup']} already present")
    return True, (f"{msg} | build_reactivity_evidence: {report['n_evidence']} evidence row(s) "
                  f"from {', '.join(bits) or 'no sheets'}"
                  + ("; " + "; ".join(extra) if extra else ""))


def derive_mutation_specific(mutant_reactive, wt_reactive, *, measured: bool):
    """Faithfully derive `ExtractedEvidence.mutation_specific` (mutant-preferential T-cell response)
    from an observed mutant-vs-WT comparison. SOURCE-AGNOSTIC: any adapter (the 33064988 Supp T5
    PDF, a future xlsx/docx) feeds the two observed reactivities through this ONE primitive, so the
    prediction-vs-evidence boundary is enforced in a single tested place rather than per table.

    Returns True (mutant-preferential), False (cross-reactive: mutant == WT), or None (not derivable).
    It REFUSES (None) unless there is a MEASURED two-arm comparison — specificity is a comparative
    claim, so a single arm or an in-silico prediction can never establish it (asserting from those
    would be fabrication). Contract (reviewed w/ Fable 5, 2026-06-09):
      - measured=False (in-silico / not a measured assay)  -> None  (predictions are never evidence)
      - mutant_reactive is None                            -> None  (mutant arm not observed)
      - mutant_reactive is False                           -> None  (no positive response -> specificity N/A;
                                                                      the negative finding lives on the
                                                                      immunogenicity axis, not here)
      - wt_reactive is None                                -> None  (single-arm: WT not tested)
      - else -> (mutant_reactive and not wt_reactive)

    LIMITATION (documented): the bool collapses a quantitative preference — "mutant >> WT but WT still
    weakly positive" becomes False. Fine for a binary cross-reactivity flag (33064988); a magnitude-aware
    refinement is future work IF a paper reports graded mutant-vs-WT magnitudes."""
    if not measured:
        return None
    if mutant_reactive is None or mutant_reactive is False:
        return None
    if wt_reactive is None:
        return None
    return bool(mutant_reactive) and not bool(wt_reactive)


_CROSSREACT_RX = re.compile(
    r'^(\S+)\s+([ACDEFGHIKLMNPQRSTVWY]{8,})\s+([ACDEFGHIKLMNPQRSTVWY]{8,})\s+(Yes|No)\b', re.I)


def _parse_crossreactivity_pdf(pdf_path: str) -> list:
    """Parse a mutant-vs-WT cross-reactivity table out of a PDF (33064988 Supp Table 5,
    mmc1.pdf p5: 'Peptide ID | Mutant seq | Wt seq | Cross reactive to Wt'). Scans every page
    line-by-line for the 4-field row shape (id + two AA sequences + Yes/No), so it finds the
    table without page config. Returns [(peptide_id, mutant, wt, cross_reactive_bool), ...]."""
    PdfReader = _require("pypdf").PdfReader
    rows, seen = [], set()
    for pg in PdfReader(pdf_path).pages:
        for line in (pg.extract_text() or "").splitlines():
            m = _CROSSREACT_RX.match(line.strip())
            if not m:
                continue
            pid, mut, wt, cross = m.group(1), m.group(2).upper(), m.group(3).upper(), m.group(4)
            if pid in seen:
                continue
            seen.add(pid)
            rows.append((pid, mut, wt, cross.lower() == "yes"))   # cross=True iff reactive to WT
    return rows


def load_epitope_manifest_into(out_path: str, path: str, sheet: str | None = None,
                               section_ref: str | None = None) -> tuple[bool, str]:
    """DETERMINISTIC loader for a class-I/II prediction / epitope MANIFEST (one row per predicted minimal
    epitope: an HLA/MHC-allele column + a mutant-epitope column, usually with the long immunizing-peptide
    sequence). Reads EVERY row and BOTH MHC-I and MHC-II epitope columns via `epitope_manifest_adapter`
    (a pure function), so the epitope + immunizing_peptide lanes are IDENTICAL every run.

    WHY: hand-transcribing the minimal-epitope layer collapses under budget — the same Keskin manifest
    yielded 2 vs 99 epitopes run-to-run (stability lab, 2026-07-05). This replaces that transcription with
    one call whose output is byte-identical and gold-exact on the calibrated manifests (Keskin 99/103,
    Rojas 451/223). MERGE-dedupes into whatever is already loaded — epitopes by (sequence, hla_allele,
    mhc_class), immunizing peptides by sequence — so it is safe to combine with peptides loaded elsewhere
    and is idempotent. Call ONCE for the manifest sheet INSTEAD of add_table/add_entities for those lanes.

    BOUNDARY: a manifest whose minimal epitope is NOT an explicit column (e.g. Li 33879241 Table S2 lists
    only the long 21-mer; the minimal epitope is a derived NetMHC substring) yields peptides only and NO
    epitopes — the tool loads what the columns state and never invents a predicted binder.

    LINKAGE (2026-08-16): the adapter emits epitopes with no parent_peptide_ids, so this wrapper runs
    `link_orphan_epitopes` after the merge — each loaded epitope gets the immunizing peptide(s) whose
    sequence CONTAINS it. That must happen here because the epitopes lane is locked below, which would
    otherwise leave a 100%-orphaned lane the agent cannot repair (the tcr_phase5 failure)."""
    try:
        from . import epitope_manifest_adapter as _ema
    except Exception as e:  # pragma: no cover - import guard
        return False, f"epitope_manifest_adapter import failed: {e}"
    try:
        loaded = _ema.load_epitope_manifest(path, sheet=sheet)
    except Exception as e:
        return False, f"could not read epitope manifest {path}" + (f" sheet={sheet!r}" if sheet else "") + f": {e}"
    eps = loaded.get("epitopes") or []
    imps = loaded.get("immunizing_peptides") or []
    if not eps and not imps:
        return False, (f"no epitope/immunizing-peptide manifest recognized in {path}"
                       + (f" sheet={sheet!r}" if sheet else "")
                       + " — need a header with an HLA/MHC-allele column plus a mutant-epitope or "
                         "immunizing-peptide column (this is not that kind of sheet)")
    rec, err = _load_partial(out_path)
    if err:
        return False, err
    if rec is None:
        return False, "no partial record; call init_record first"
    if section_ref:
        for x in (*eps, *imps):
            x["section_ref"] = section_ref

    report = []
    for section, items in (("immunizing_peptides", imps), ("epitopes", eps)):
        model = getattr(schema, SECTION_MODEL[section])
        normalized = []
        for i, it in enumerate(items):
            try:
                normalized.append(json.loads(model(**it).model_dump_json()))
            except Exception as e:
                return False, (f"manifest {section} row {i} invalid:\n{e}\n"
                               "(nothing added; partial unchanged)")
        # One identity rule for the whole codebase (`_ENTITY_IDENTITY`) -- this merge is where it
        # was first written, and the freehand append path now shares it.
        rec.setdefault(section, [])
        fresh, dropped = _dedupe_by_identity(section, rec[section], normalized)
        rec[section].extend(fresh)
        _record_id_aliases(out_path, dropped)
        report.append(f"{section} +{len(fresh)} of {len(items)} (total {len(rec[section])})")
    # The adapter emits manifest epitopes with NO parent_peptide_ids, and this lane is about to be
    # LOCKED -- so an unlinked epitope here is unfixable by the agent later. Derive the parents now,
    # from the peptides the record holds (this loader's own long-peptide column included). Finalize
    # repeats this, so peptides loaded after the manifest still rescue any orphan left here.
    n_linked = link_orphan_epitopes(rec)
    still = sum(1 for e in rec.get("epitopes") or [] if not e.get("parent_peptide_ids"))
    report.append(f"parent linkage: derived for {n_linked} epitope(s)"
                  + (f"; {still} STILL ORPHANED (no loaded immunizing peptide contains their sequence "
                     f"-- load that peptide lane, or they are quarantined at finalize)" if still else ""))
    pathlib.Path(_partial_path(out_path)).write_text(json.dumps(rec))
    # sticky: lock only the lanes the loader actually populated (Li-shape loads peptides-only, so the
    # agent must still be free to MINT epitopes -- don't lock an empty lane).
    _lock_lanes(out_path, *[lane for lane in ("epitopes", "immunizing_peptides") if rec.get(lane)])
    return True, "loaded epitope manifest: " + "; ".join(report)


def _clonotype_chain_key(c: dict):
    """The stable identity the ruler uses: the frozenset of (locus, CDR3) over a clonotype's chains."""
    return frozenset((ch.get("locus"), ch.get("junction_aa") or ch.get("cdr3_aa") or ch.get("cdr3_nt"))
                     for ch in (c.get("chains") or []) if (ch.get("junction_aa") or ch.get("cdr3_aa") or ch.get("cdr3_nt")))


# ---------------------------------------------------------------------------------------------
# FLAGSHIP flag DERIVATION from a reactive-TCR sheet (2026-08-16). `add_reactive_tcr` is already called
# with (patient, antigen) and its sheet often names the sorted subset -- that tuple IS a
# NeoantigenTcrFlag. Deriving it makes the flag gate cheap to satisfy WITHOUT hand-work.
#
# FABRICATION BOUNDARY (deliberately narrow -- a flag on a guessed target is worse than no flag, and
# schema `_check_target` rejects it anyway). A flag is minted ONLY when BOTH are unambiguous:
#   (1) the sorted SUBSET is an explicit CD4/CD8 token in the loader's own inputs (antigen / sheet name /
#       section_ref). `TcrIdentified` has no "found, subset unknown" value, so an unlabelled sheet cannot
#       be encoded without guessing -- we skip instead. (Keskin's Pt8 sheets are unlabelled; phase-4's
#       'cd4' there came from Fig. 3b prose, NOT the sheet.)
#   (2) the `antigen` string resolves to EXACTLY ONE target entity already in the record, by
#       paper_local_id / gene_symbol / sequence, searched immunizing_peptides -> epitopes -> candidates
#       (the reactive sort is against the STIMULATING peptide). Zero or several matches -> skip and TELL
#       the agent the candidate ids, so gate `_tcr_flag_gap` + `partial_status` can finish the job.
# The lane is NEVER locked: Rojas's real flag came from Fig. 2d, which no sheet can supply.
_FLAG_SUBSET_CD4 = re.compile(r"\bCD\s*-?\s*4\b", re.I)
_FLAG_SUBSET_CD8 = re.compile(r"\bCD\s*-?\s*8\b", re.I)
# Words that describe the ASSAY, not the antigen -- dropped before target matching so 'mutant SHANK2'
# and 'Pt7 Neoantigen-reactive CD4+' reduce to their identifying tokens (or to nothing).
_FLAG_ANTIGEN_STOPWORDS = frozenset({
    "mut", "mutant", "mutated", "mutation", "wt", "wildtype", "neoantigen", "neoantigens", "neoepitope",
    "neoepitopes", "antigen", "antigens", "epitope", "epitopes", "peptide", "peptides", "long", "minimal",
    "reactive", "reactivity", "specific", "specificity", "sorted", "sort", "pool", "pools", "pooled",
    "cell", "cells", "tcell", "tcells", "t", "cd4", "cd8", "positive", "neg", "negative", "tumor",
    "tumour", "associated", "relapse", "recurrent", "vaccine", "vaccinated", "and", "the", "of",
})


def _flag_antigen_tokens(antigen: str) -> list[str]:
    """Identifying tokens of a loader `antigen` string, assay words dropped. 'mutant SHANK2' -> ['SHANK2'];
    'mutSHANK2' -> ['mutSHANK2','SHANK2']; 'pool C neoantigen (CD4)' -> ['C'] (which resolves to nothing)."""
    out, seen = [], set()
    for tok in re.split(r"[^A-Za-z0-9\-\*:_]+", antigen or ""):
        for cand in (tok, re.sub(r"^mut[-_]?", "", tok, flags=re.I), re.sub(r"[-_]?mut$", "", tok, flags=re.I)):
            cand = cand.strip("-_")
            if not cand or cand.lower() in _FLAG_ANTIGEN_STOPWORDS or cand.lower() in seen:
                continue
            seen.add(cand.lower())
            out.append(cand)
    return out


def _resolve_flag_target(rec: dict, antigen: str) -> tuple[tuple[str, str] | None, str]:
    """(kind, paper_local_id) for `antigen`, or (None, why-not). EXACTLY-ONE-match rule, per kind, in
    tier order -- a tie is a refusal, never a pick."""
    tokens = _flag_antigen_tokens(antigen)
    if not tokens:
        return None, f"antigen {antigen!r} names no identifying token (only assay words)"
    upper = {t.upper() for t in tokens}
    for kind, lane in (("immunizing_peptide", "immunizing_peptides"), ("epitope", "epitopes"),
                       ("candidate", "candidates")):
        hits = []
        for row in rec.get(lane) or []:
            pid = row.get("paper_local_id")
            gene = (row.get("gene_symbol") or "").upper()
            seq = (row.get("sequence") or "").upper()
            if pid in tokens or (gene and gene in upper) or (seq and seq in upper):
                hits.append(pid)
        if len(hits) == 1:
            return (kind, hits[0]), ""
        if hits:
            return None, (f"antigen {antigen!r} matches {len(hits)} {kind}s ({', '.join(hits[:6])}"
                          f"{', ...' if len(hits) > 6 else ''}) -- ambiguous, not guessing which")
    return None,(f"antigen {antigen!r} matches no immunizing_peptide / epitope / candidate in the record "
                  f"(searched by paper_local_id, gene_symbol and sequence)")


def _flag_subset(*texts: str) -> tuple[str | None, str]:
    """('cd4'|'cd8', the source token) from the loader's own strings, or (None, why-not). BOTH tokens
    present -> refuse (a mixed sheet is not one flag); NEITHER -> refuse (no honest TcrIdentified value)."""
    blob = " ".join(t for t in texts if t)
    cd4, cd8 = bool(_FLAG_SUBSET_CD4.search(blob)), bool(_FLAG_SUBSET_CD8.search(blob))
    if cd4 and cd8:
        return None, "the sheet names BOTH CD4 and CD8 -- one flag cannot represent a mixed sort"
    if cd4:
        return "cd4", "CD4"
    if cd8:
        return "cd8", "CD8"
    return None, ("neither CD4 nor CD8 appears in the antigen / sheet name / section_ref, and "
                  "`tcr_identified` has no 'found, subset unknown' value")


def derive_reactive_tcr_flags(rec: dict, *, patient_paper_id: str, antigen: str, sheet: str,
                              section_ref: str, n_clonotypes: int) -> tuple[list[dict], str]:
    """One `add_reactive_tcr` call -> ([flag dict] or []), plus a human reason when nothing was derived.
    PURE w.r.t. the record (reads it, returns rows; the caller appends). Never mints an unresolvable
    target, never guesses a subset, never invents provenance (quoted_text/section_ref come from the sheet
    this loader actually read)."""
    if patient_paper_id not in {p.get("paper_local_id") for p in rec.get("patients") or []}:
        return [], (f"patient {patient_paper_id!r} is not in the record's patients lane "
                    f"(add the patient, then the flag can be derived on a re-run)")
    subset, subset_why = _flag_subset(antigen, sheet, section_ref)
    target, target_why = _resolve_flag_target(rec, antigen)
    if subset is None or target is None:
        why = "; ".join(w for w, bad in ((subset_why, subset is None), (target_why, target is None)) if bad)
        return [], why
    kind, tid = target
    existing = {(f.get("patient_paper_id"), f.get("target_kind"),
                 f.get("immunizing_peptide_paper_id") or f.get("epitope_paper_id")
                 or f.get("candidate_paper_id")): f.get("tcr_identified")
                for f in rec.get("neoantigen_tcr_flags") or []}
    prior = existing.get((patient_paper_id, kind, tid))
    if prior is not None:
        # never overwrite a row already in the lane; if this sheet sorted the OTHER subset, say so --
        # 'both' is the agent's call to make (it can see whether the two sorts are the same claim).
        return [], (f"a flag for ({patient_paper_id}, {tid}) already exists (tcr_identified={prior!r}) "
                    f"-- left as it is" + (f"; this sheet sorted {subset!r}, so the truth may be 'both' "
                                           f"-- update that row if the paper supports it"
                                           if prior != subset else ""))
    quoted = (f"{patient_paper_id} {antigen}-reactive single-cell TCR sheet [{sheet}]: {n_clonotypes} "
              f"clonotype(s) recovered from cells sorted as {subset_why}+")
    flag = schema.NeoantigenTcrFlag(
        target_kind=kind,
        immunizing_peptide_paper_id=tid if kind == "immunizing_peptide" else None,
        epitope_paper_id=tid if kind == "epitope" else None,
        candidate_paper_id=tid if kind == "candidate" else None,
        patient_paper_id=patient_paper_id,
        tcr_identified=subset,
        method_label="add_reactive_tcr (single-cell)",
        quoted_text=quoted[:2000],
        section_ref=section_ref,
        needs_review=False,      # derived from an antigen-sorted sheet; H3 is inferred_association, not this
        confidence=2,
    )
    return [json.loads(flag.model_dump_json())], ""


def load_reactive_tcr_into(out_path: str, path: str, sheet: str, patient_paper_id: str, antigen: str,
                           section_ref: str | None = None, min_cells: int = 1) -> tuple[bool, str]:
    """DETERMINISTIC loader for a paired single-cell 'neoantigen-reactive T cell' sheet (Keskin MOESM7/10
    shape: one row per cell, columns TRAV/TRAJ/Alpha CDR3 aa/TRBV/TRBJ/Beta CDR3 aa/Clone/Clonotype).
    One clonotype per unique (CDR3a, CDR3b) pair with cell_count >= `min_cells`; reads via
    `tcr_reactive_adapter` (a pure function), so the clonotype set is IDENTICAL every run.

    WHY: the live agent extracts the TCR lane non-deterministically (4 clonotypes one run, 0 the next) and
    captures only the dominant clones. This loads the full sorted-reactive set instead. Merge-dedupes into
    the record by chain identity (locus, CDR3) so it is safe to call once per reactive sheet and idempotent.
    tcr_seq_status is set to 'single_cell' when clonotypes land (unless a broader value is already set)."""
    try:
        from . import tcr_reactive_adapter as _tra
    except Exception as e:  # pragma: no cover
        return False, f"tcr_reactive_adapter import failed: {e}"
    try:
        clonos = _tra.load_reactive_clonotypes(
            path, sheet, patient_paper_id=patient_paper_id, antigen=antigen,
            section_ref=section_ref or f"{patient_paper_id} {antigen}-reactive single-cell TCR",
            min_cells=min_cells)
    except Exception as e:
        return False, f"could not read reactive-TCR sheet {path} [{sheet}]: {e}"
    if not clonos:
        # LEDGER an attempt that reached the adapter and came back empty: the sheet WAS routed, so a
        # coverage gate must not demand it again (that block would be unfixable), and the 0/0 row is
        # the auditable record of the empty answer.
        _append_load_ledger(out_path, path=path, sheet=sheet, loader="add_reactive_tcr",
                            n_loaded=0, n_available=0)
        return False, (f"no clonotypes with >= {min_cells} cells found in {path} [{sheet}] "
                       "(not a paired single-cell reactive-clone sheet, or all singletons below the cutoff)")
    rec, err = _load_partial(out_path)
    if err:
        return False, err
    if rec is None:
        return False, "no partial record; call init_record first"
    seen = {_clonotype_chain_key(c) for c in rec.get("tcr_clonotypes") or []}
    fresh = []
    for c in clonos:
        cd = json.loads(c.model_dump_json())
        k = _clonotype_chain_key(cd)
        if k in seen:
            continue
        seen.add(k)
        fresh.append(cd)
    rec.setdefault("tcr_clonotypes", []).extend(fresh)
    if rec.get("tcr_seq_status") in (None, "none", "mentioned_only"):
        rec["tcr_seq_status"] = "single_cell"
    # FLAGSHIP: derive the per-(neoantigen x patient) flag when this sheet unambiguously supplies BOTH the
    # sorted subset and a target that resolves in the record. Never locks the lane -- the agent must stay
    # free to add text-derived flags (Rojas's only flag came from Fig. 2d).
    sref = section_ref or f"{patient_paper_id} {antigen}-reactive single-cell TCR"
    try:
        from . import tcr_target_bind as _ttb
        if _ttb.antigen_is_mixed(rec, antigen) and _ttb.sheet_row_antigen_col(path, sheet) is None:
            bind_note = f"; mixed antigen {antigen!r} unbound (no per-row antigen column)"
        elif not _ttb.antigen_is_mixed(rec, antigen):
            for cd in fresh:
                _ttb.bind_clonotype(
                    rec, cd, antigen=antigen, section_ref=sref,
                    patient_paper_id=patient_paper_id, sheet=sheet, xlsx_path=path,
                )
            bind_note = ""
        else:
            bind_note = f"; mixed antigen {antigen!r} has a per-row column but clones are grouped by CDR3 — left unbound"
        _ttb.apply_reporter_refs(rec, clonotypes=fresh)
    except Exception:  # never let bind break a load
        bind_note = ""
    try:
        flags, why = derive_reactive_tcr_flags(rec, patient_paper_id=patient_paper_id, antigen=antigen,
                                               sheet=sheet, section_ref=sref, n_clonotypes=len(clonos))
    except Exception as e:                                      # never let derivation break the load
        flags, why = [], f"flag derivation failed ({e})"
    if flags:
        rec.setdefault("neoantigen_tcr_flags", []).extend(flags)
    pathlib.Path(_partial_path(out_path)).write_text(json.dumps(rec))
    _lock_lanes(out_path, "tcr_clonotypes")                     # sticky: the loader owns this lane
    _append_load_ledger(out_path, path=path, sheet=sheet, loader="add_reactive_tcr",
                        n_loaded=len(fresh), n_available=len(clonos))
    note = (f"; DERIVED {len(flags)} neoantigen_tcr_flags row(s) (lane NOT locked -- add any the text "
            f"supports)" if flags else
            f"; NO neoantigen_tcr_flags derived ({why}) -- that FLAGSHIP lane is yours: call "
            f"partial_status for the real target ids and add_entities one row per (neoantigen x patient)")
    return True, (f"loaded {len(fresh)} reactive clonotypes (of {len(clonos)}) from [{sheet}] "
                  f"-> tcr_clonotypes total {len(rec['tcr_clonotypes'])}; "
                  f"tcr_seq_status={rec['tcr_seq_status']}" + note + bind_note)


def load_tcr_track_into(out_path: str, path: str, sheet: str, patient_paper_id: str,
                        section_ref: str | None = None, locus: str | None = None,
                        tissue: str = "blood", sorted_compartment: str = "whole_pbmc",
                        header_row_index: int = 1) -> tuple[bool, str]:
    """DETERMINISTIC loader for a WIDE tet-spec TRACKING sheet (Hu 33479501 'Suppl 8 Pt<N> Tet-spec
    TCR<a|b> track' shape: col0 = the clone identity, col1 = a tetramer-specific flag, col2.. = that
    clone's frequency at each timepoint). Keeps ONLY the flagged rows and UNPIVOTS every timepoint column
    into the clonotype's nested `observations`, via `tcr_track_adapter` (a pure function) -- so both the
    clonotype set AND the trajectories are IDENTICAL every run. Handles both identity-cell forms that
    occur in ONE supplement (a bare CDR3, or the paper's comma-joined 'V,J,CDR3' token -> v_call/j_call
    plus the verbatim token in chain.raw) and both loci (TCRa and TCRb sheets exist per patient; `locus`
    is inferred from the identity column's header, and only needs passing when that header names neither
    chain -- it is never guessed).

    WHY: the live agent is non-deterministic on this table class (Hu swung 0 -> 39 -> 20 clonotypes
    run-to-run, 39972124 0 -> 12) and budget-shortcuts dense sheets, capturing a clone's origin-frequency
    snapshot but NOT its per-timepoint trajectory. Both failures are structural to LLM transcription of a
    large clean columnar supplement. Merge-dedupes by (patient, chain identity) so it is safe to call once
    per patient sheet and idempotent. tcr_seq_status is raised to 'bulk' ('both' when single-cell clonotypes
    were already loaded) -- a bulk repertoire track is what this sheet is.

    BOUNDARY: this is the WIDE one-row-per-CLONE-by-timepoint shape. A paired single-cell one-row-per-CELL
    reactive sheet is `add_reactive_tcr` instead; `add_table`'s flat column->field DSL can serve NEITHER
    (it cannot unpivot). Rows without the tetramer-specific flag are dropped, not loaded unflagged -- the
    tool records the paper's own specificity call and never infers one."""
    try:
        from . import tcr_track_adapter as _tta
    except Exception as e:  # pragma: no cover - import guard
        return False, f"tcr_track_adapter import failed: {e}"
    try:
        clonos = _tta.load_track_sheet(path, sheet, patient_paper_id=patient_paper_id,
                                       header_row_index=header_row_index, locus=locus, tissue=tissue,
                                       sorted_compartment=sorted_compartment)
    except Exception as e:
        return False, f"could not read tet-spec track sheet {path} [{sheet}]: {e}"
    if not clonos:
        # LEDGER the empty answer (see the twin comment in load_reactive_tcr_into).
        _append_load_ledger(out_path, path=path, sheet=sheet, loader="add_tcr_track",
                            n_loaded=0, n_available=0)
        return False, (f"no tetramer-specific tracked clonotypes found in {path} [{sheet}] — need a wide "
                       "CDR3-by-timepoint sheet with a tetramer-specific flag column set to 1 (this is not "
                       f"that kind of sheet, or header_row_index={header_row_index} is not the label row)")
    rec, err = _load_partial(out_path)
    if err:
        return False, err
    if rec is None:
        return False, "no partial record; call init_record first"
    # Dedupe key includes the PATIENT: track sheets are one per patient and a public CDR3b legitimately
    # recurs across patients, so chain-identity alone would silently drop another patient's real clone.
    # (add_reactive_tcr's paired CDR3a+CDR3b keys cannot collide with these TRB-only ones.)
    seen = {(c.get("patient_paper_id"), _clonotype_chain_key(c)) for c in rec.get("tcr_clonotypes") or []}
    fresh = []
    for c in clonos:
        cd = json.loads(c.model_dump_json())
        if section_ref:
            cd["section_ref"] = section_ref
        k = (cd.get("patient_paper_id"), _clonotype_chain_key(cd))
        if k in seen:
            continue
        seen.add(k)
        fresh.append(cd)
    rec.setdefault("tcr_clonotypes", []).extend(fresh)
    cur = rec.get("tcr_seq_status")
    if cur in (None, "none", "mentioned_only"):
        rec["tcr_seq_status"] = "bulk"
    elif cur == "single_cell":
        rec["tcr_seq_status"] = "both"          # never downgrade a sticker that is already broader
    pathlib.Path(_partial_path(out_path)).write_text(json.dumps(rec))
    _lock_lanes(out_path, "tcr_clonotypes")                     # sticky: the loader owns this lane
    _append_load_ledger(out_path, path=path, sheet=sheet, loader="add_tcr_track",
                        n_loaded=len(fresh), n_available=len(clonos))
    n_obs = sum(len(c.get("observations") or []) for c in fresh)
    # A track sheet names NO antigen and NO CD4/CD8 subset (col0 is a CDR3, col1 a tetramer flag), so a
    # flag cannot be derived from it without inventing a target -- say so instead of silently leaving the
    # flagship lane empty. The finalize flag-gate will insist on it.
    return True, (f"loaded {len(fresh)} tetramer-specific tracked clonotypes (of {len(clonos)}) carrying "
                  f"{n_obs} timepoint observations from [{sheet}] -> tcr_clonotypes total "
                  f"{len(rec['tcr_clonotypes'])}; tcr_seq_status={rec['tcr_seq_status']}"
                  f"; NO neoantigen_tcr_flags derived (a track sheet names neither the antigen nor the "
                  f"CD4/CD8 subset) -- that FLAGSHIP lane is yours: which neoantigen the tetramer carried "
                  f"is in the text/figures; call partial_status for the real target ids")


def build_crossreactivity_evidence(out_path: str, pdf_path: str,
                                   section_ref: str = "Supplemental Table 5; Figure 4B",
                                   assay: str = "elispot",
                                   provenance_locator: str = "Supplemental Table 5") -> tuple[bool, str]:
    """PAPER-SPECIFIC ADAPTER (33064988 Supplemental Table 5). Deterministically load that mutant-vs-WT
    cross-reactivity TABLE into one MinimalEpitope + one epitope-target evidence row per listed peptide.
    The table is READABLE (text PDF), but the agent hand-built these via add_entities and dropped them
    run-to-run; this parses the 13 rows once, identically every run.

    SCOPE (reviewed w/ Fable 5): this is an ADAPTER for ONE paper's table layout, NOT a general parser —
    the immunogenicity cross-reactivity *outcome*-table pattern occurs exactly once in the corpus (the
    other papers' mutant-vs-WT tables are prediction MANIFESTS, handled by add_table -> candidates/epitopes;
    routing those here would be a prediction->evidence category error). The reusable, source-agnostic part
    is `derive_mutation_specific`, which this calls. A SECOND such paper is a triggered roadmap item — add a
    sibling `build_crossreactivity_evidence_<pmid>` adapter then (do NOT parameterize this one speculatively).

    Per row: epitope sequence = the MUTANT peptide, wild_type_sequence = the WT, parent =
    the already-loaded immunizing peptide whose sequence CONTAINS the mutant (patient-token
    matched), mhc_class inferred from length (<=11 -> I else II; epitope -> needs_review since the
    class is conventional, not stated). Evidence: outcome=immunogenic, mutation_specific = the row
    is NOT cross-reactive to WT, provenance kind='table'. A row whose parent peptide isn't loaded
    is SKIPPED (no orphan epitope invented)."""
    rec, err = _load_partial(out_path)
    if err:
        return False, err
    if rec is None:
        return False, "no partial record; call init_record first"
    peps = rec.get("immunizing_peptides") or []
    if not peps:
        return False, "no immunizing_peptides loaded yet; load peptides BEFORE build_crossreactivity_evidence"
    try:
        rows = _parse_crossreactivity_pdf(pdf_path)
    except Exception as e:
        return False, f"could not parse {pdf_path}: {e}"
    if not rows:
        return False, f"no mutant-vs-WT cross-reactivity rows found in {pdf_path}"

    def _patient(pid):
        m = re.match(r'^([A-Za-z]+\d+)', pid)
        return m.group(1) if m else None

    def _parent_for(patient, mutant):
        tok = re.compile(rf"(?<![A-Za-z0-9]){re.escape(patient)}(?![A-Za-z0-9])")
        cands = [p for p in peps
                 if p.get("paper_local_id") and tok.search(p["paper_local_id"])
                 and mutant in _norm_seq(p.get("sequence"))]
        if not cands:
            return None
        exact = [p for p in cands if _norm_seq(p.get("sequence")) == mutant]
        pick = exact[0] if exact else min(cands, key=lambda p: len(_norm_seq(p.get("sequence"))))
        return pick["paper_local_id"]

    epitopes, evidence, skipped = [], [], []
    for pid, mut, wt, cross in rows:
        patient = _patient(pid)
        if not patient:
            skipped.append(pid)
            continue
        parent = _parent_for(patient, mut)
        if not parent:
            skipped.append(pid)                       # no loaded parent -> do NOT invent an orphan epitope
            continue
        ep_id = f"EP-{pid}"
        mhc_class = "I" if len(mut) <= 11 else "II"   # length convention; declared inferred + needs_review
        ep = {
            "paper_local_id": ep_id, "sequence": mut, "wild_type_sequence": wt,
            "is_neoantigen": True,                        # mutant neoantigen epitope (mutant != WT)
            "parent_peptide_ids": [parent], "mhc_class": mhc_class,
            # the class is a LENGTH HEURISTIC, not a source-stated restriction: declare it explicitly
            # (mhc_class_inferred) so the schema's class-II anchor exempts this AUDITED call instead of us
            # laundering a synthetic class cue into quoted_text. needs_review routes it to a human.
            "mhc_class_inferred": True, "needs_review": True,
            "quoted_text": (f"Mutant epitope {mut} (WT {wt}); "
                            f"{'mutant-preferential' if not cross else 'cross-reactive to WT'} "
                            "by IFN-gamma ELISpot across a titration. MHC class inferred from length "
                            "(heuristic, NOT a source-stated restriction; flagged for review)."),
            "section_ref": section_ref,
        }
        if mhc_class == "I":   # the schema requires an affinity on class-I; the table omits it -> lossless
            ep["predicted_affinity"] = {"unit": "unknown", "tier": "reported",
                                        "raw": "affinity not reported (Supplemental Table 5)"}
        epitopes.append(ep)
        evidence.append({
            "patient_paper_id": patient, "target_kind": "epitope", "epitope_paper_id": ep_id,
            "assay": assay, "outcome": "immunogenic", "vaccine_induced": True,
            # mutant_reactive=True (the epitope IS immunogenic here); wt_reactive = the cross-reactivity
            # flag. Routed through the shared primitive so the measured two-arm rule is enforced once.
            "mutation_specific": derive_mutation_specific(True, cross, measured=True),
            "magnitude": {"unit": "unknown", "tier": "reported",
                          "raw": (f"mutant {mut} "
                                  + ("preferential over WT " + wt + " (mutation-specific)" if not cross
                                     else "cross-reactive to WT " + wt + " (not mutant-preferential)")
                                  + "; absolute SFC not tabulated (Supplemental Table 5)")},
            "quoted_text": f"{pid}: mutant {mut} vs WT {wt} — cross-reactive to WT: {'Yes' if cross else 'No'}.",
            "section_ref": section_ref,
            "provenance": [{"kind": "table", "locator": provenance_locator}],
        })
    if not epitopes:
        return False, (f"0 cross-reactivity rows mapped ({len(rows)} parsed, {len(skipped)} without a "
                       "loaded parent peptide). Load the immunizing peptides first.")
    ok, emsg = _validate_and_append(out_path, "epitopes", epitopes)
    if not ok:
        return False, f"epitope append failed: {emsg}"
    ok, vmsg = _validate_and_append(out_path, "evidence", evidence)
    if not ok:
        return False, f"epitopes added but evidence append failed: {vmsg}"
    note = f"; {len(skipped)} rows had no loaded parent (skipped)" if skipped else ""
    return ok, (f"build_crossreactivity_evidence: {len(epitopes)} epitopes + {len(evidence)} "
                f"epitope/immunogenic evidence rows from {len(rows)} parsed rows{note}")


# The fields that tell two measurements of the same patient x target apart. A bulk evidence load
# whose rows collide on `_evidence_measurement_key` will be silently thinned by finalize's dedup --
# how Hu 33479501 run2 loaded 139 rows and kept 99, every one of them with an empty timepoint.
_MEASUREMENT_DISCRIMINATORS = ("timepoint_phase", "assay_stimulation", "assay_antigen_format",
                               "assay_detail", "assay", "t_cell_subset", "mhc_class")


def _colliding_measurements(items: list) -> str:
    """'' if every row is a distinct measurement, else an error naming a colliding pair.

    Rejecting here (rather than letting dedup eat them at finalize) is the difference between a
    loud, fixable load and a record that quietly holds a fraction of its source sheet."""
    seen: dict = {}
    for i, it in enumerate(items):
        if not isinstance(it, dict):
            continue
        key = _evidence_measurement_key(it)
        if key in seen:
            first = items[seen[key]]
            unset = [f for f in _MEASUREMENT_DISCRIMINATORS if not it.get(f)]
            n_dup = len(items) - len({_evidence_measurement_key(x) for x in items
                                      if isinstance(x, dict)})
            hint = (f" Every one of these is unset on the row: {unset}. Give each per_cell rule its "
                    "own timepoint_phase / assay_stimulation / assay_antigen_format so the "
                    "measurements stay distinct." if unset else
                    " Their discriminators are already set and identical, so the source really does "
                    "list this measurement twice — re-run with "
                    '"allow_duplicate_measurements": true in the mapping to load it anyway.')
            target = (it.get("epitope_paper_id") or it.get("immunizing_peptide_paper_id")
                      or it.get("candidate_paper_id") or it.get("pool_paper_id"))
            shown = {k: v for k, v in first.items() if k in _MEASUREMENT_DISCRIMINATORS}
            return (f"{n_dup} of {len(items)} generated evidence row(s) are DUPLICATE measurements "
                    f"and finalize's dedup would silently drop them; nothing was added. Rows "
                    f"{seen[key]} and {i} share one measurement identity "
                    f"(patient={it.get('patient_paper_id')!r}, target={target!r}, "
                    f"outcome={it.get('outcome')!r}).{hint} [first row: {shown}]")
        seen[key] = i
    return ""


def table_to_entities(out_path: str, section: str, path: str,
                      mapping_json: str, sheet: str | None = None,
                      header_row: int | None = None,
                      sheets: list | None = None) -> tuple[bool, str]:
    """Bulk-add entities to `section` by mapping xlsx columns to schema fields.
    Reads ALL rows, applies an optional row filter, builds each entity via the closed
    mapping DSL, and appends them atomically (any invalid row -> nothing appended).

    MULTI-SHEET (`sheets`): pass a list of sheet names to apply the SAME mapping across many
    per-entity sheets in ONE call — the cheap-compliance path for papers that split data into
    one sheet per patient (e.g. 33064988's ~34 'IAP-<patient>' immunogenicity tabs, which a
    per-sheet sweep left under-extracted). Every row gets a reserved `__sheet__` column holding
    its sheet name, so a field can be derived from it, e.g.
    `{"patient_paper_id": {"col": "__sheet__", "extract": "IAP-(.+)"}}` -> 'IAP-M13' becomes 'M13'.
    Append stays atomic across ALL sheets (any invalid row anywhere -> nothing appended)."""
    if section not in SECTION_MODEL:
        return False, f"unknown section {section!r}; valid: {list(SECTION_MODEL)}"
    rec, err = _load_partial(out_path)
    if err:
        return False, err
    if rec is None:
        return False, "no partial record; call init_record first"
    try:
        mapping = json.loads(mapping_json)
    except Exception as e:
        return False, f"BAD mapping JSON: {e}"
    if not isinstance(mapping, dict) or not isinstance(mapping.get("fields"), dict):
        return False, "mapping must be a JSON object with a 'fields' object"
    bad_rules = [k for k, v in mapping["fields"].items() if not isinstance(v, dict)]
    if bad_rules:
        return False, f"field rule(s) {bad_rules} must be objects like {{'col': 'ColName'}}"
    if mapping.get("filter") is not None and not isinstance(mapping["filter"], dict):
        return False, "filter must be a JSON object like {'col': 'H', 'in': [...]}"
    flt = mapping.get("filter")
    per_cell = mapping.get("per_cell")
    if per_cell is not None:
        err = table_map.validate_per_cell(per_cell)
        if err:
            return False, f"BAD per_cell: {err}"
    targets = list(sheets) if sheets else [sheet]   # multi-sheet mode iff `sheets` given
    multi = bool(sheets)
    items: list = []
    per_sheet: list[str] = []
    loaded_per_sheet: list[tuple] = []
    for sh in targets:
        try:
            headers, rows = _read_xlsx_dicts(path, sh, header_row)
        except Exception as e:
            return False, f"ERROR reading sheet {sh!r} in {path!r}: {e}"
        n_scanned = len(rows)
        if multi:                                   # inject the reserved sheet-name column
            for r in rows:
                r["__sheet__"] = sh
            headers = list(headers) + ["__sheet__"]
        miss = table_map.missing_columns(mapping["fields"], flt, headers, per_cell)
        if miss:
            return False, f"column(s) {miss} not found in sheet {sh!r}; headers: {headers}"
        if flt:
            ok, kept = table_map.apply_filter(rows, flt)
            if not ok:
                return False, kept  # error message
            rows = kept
        try:
            if per_cell:
                # One row -> one entity PER MATCHING CELL (reactivity matrices; see table_map).
                sheet_items = [it for r in rows
                               for it in table_map.apply_mapping_cells(r, mapping["fields"], per_cell)]
            else:
                sheet_items = [table_map.apply_mapping(r, mapping["fields"]) for r in rows]
        except ValueError as e:
            return False, str(e)
        items.extend(sheet_items)
        loaded_per_sheet.append((sh or "", len(sheet_items), n_scanned))
        per_sheet.append(f"{sh}:{len(sheet_items)}")
    if not items:
        if per_cell:
            return False, (f"0 matching CELLS across {len(rows)} row(s); nothing added to {section}. "
                           "Each per_cell rule matched no cell — check the column refs and the "
                           "operator value (a '1' may be stored as the number 1, not the string).")
        hint = (_filter_zero_match_diagnostic(path, targets[-1], flt, header_row) if flt else "")
        return False, (f"0 rows after filter; nothing added to {section}"
                       + (f". {hint}" if hint else ""))
    if section == "evidence" and not mapping.get("allow_duplicate_measurements"):
        dup = _colliding_measurements(items)
        if dup:
            return False, dup
    ok, msg = _validate_and_append(out_path, section, items)
    if ok:
        # Ledger the load so the read-coverage gate can see that a deterministic loader — not a
        # hand-written batch — is what put these rows in the record.
        for sh, n_loaded, n_scanned in loaded_per_sheet:
            _append_load_ledger(out_path, path=path, sheet=sh, loader="add_table",
                                n_loaded=n_loaded, n_available=n_scanned)
        if multi:
            msg += f" | from {len(targets)} sheets [{', '.join(per_sheet)}]"
    return ok, msg
