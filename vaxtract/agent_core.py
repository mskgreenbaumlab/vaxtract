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
    orphans = [e.paper_local_id for e in rec.epitopes if not e.parent_peptide_ids]
    if not orphans:
        return None
    return (
        f"[gate] epitope linkage broken: {len(orphans)}/{len(rec.epitopes)} epitopes have "
        f"EMPTY parent_peptide_ids (orphaned -> sequences won't join to peptides). Set "
        f"parent_peptide_ids on every epitope to its immunizing_peptide paper_local_id "
        f"(in add_table use a parent_peptide_ids field with template_list matching the "
        f"immunizing_peptides id scheme). First orphans: {orphans[:3]}"
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
    return (
        f"[gate] peptide count mismatch: patients declare {declared} synthesized "
        f"(sum of n_peptides_synthesized) but the record holds {actual} immunizing_peptides "
        f"(gap {declared - actual}). Peptides were dropped or duplicated -- add the missing "
        f"ones (INCLUDING indel/frameshift neoantigens; the sequence is the short mutant "
        f"neoantigen sequence, not the long mRNA/indel context) or correct n_peptides_synthesized."
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


# v2.9 P21: the trial-constant regimen fields of VaccineDelivery (per-patient fields excluded —
# weeks_surgery_to_first_dose and n_doses_received are MEANT to vary; source is provenance).
_REGIMEN_FIELDS = (
    "adjuvant", "adjuvant_detail", "formulation_detail", "dose_amount_raw",
    "dose_per_peptide_ug", "dose_basis", "n_priming_doses", "n_boost_doses", "schedule_detail",
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


def _render_matrix(matrix: list, *, prefix: str, raw_unit: str, filtered_unit: str,
                   max_rows: int = 500, row_filter: dict | None = None,
                   columns: list | None = None, header_row: int | None = None) -> str:
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
    return (
        f"{prefix}\n{filtered_unit}: {len(data)}/{total} data rows"
        f"{note}; showing {len(shown)}{cap}:\n" + json.dumps(shown, default=str)
    )


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

    wb = openpyxl.load_workbook(path, data_only=True, rich_text=underline)
    ws = wb[sheet] if sheet else wb[wb.sheetnames[0]]
    if underline:
        matrix = [[_marked_cell(c) for c in row] for row in ws.iter_rows()]
    else:
        matrix = [list(r) for r in ws.iter_rows(values_only=True)]
    return _render_matrix(
        matrix, prefix=f"sheets={wb.sheetnames}", raw_unit=ws.title,
        filtered_unit=f"sheet {ws.title!r}", max_rows=max_rows,
        row_filter=row_filter, columns=columns, header_row=header_row)


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
        summ = "\n".join(f"  [{i}] {len(t.rows)}x{len(t.columns)}  {cap[:70]!r}"
                         for i, (t, cap) in enumerate(tables))
        return (f"docx: {len(tables)} table(s), {len(doc.paragraphs)} paragraphs.\n"
                f"tables (read one with table_index=N):\n{summ}\n"
                f"(prose: call with text_offset=0)")

    # SPECIFIC TABLE mode
    if not (0 <= table_index < len(tables)):
        return f"table_index {table_index} out of range (have {len(tables)} table(s): 0..{len(tables) - 1})"
    tbl, cap = tables[table_index]
    matrix = _docx_table_matrix(tbl, underline)
    unit = f"table[{table_index}] {cap[:70]!r}"
    return _render_matrix(
        matrix, prefix=f"tables={len(tables)}", raw_unit=unit, filtered_unit=unit,
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
}


def _partial_path(out_path: str) -> str:
    return str(out_path) + ".partial.json"


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
    warn = (f" (note: entity-list keys {dropped} in the meta were ignored; "
            "add them with add_entities)") if dropped else ""
    return True, f"INIT ok; partial at {pp}{warn}"


def _validate_and_append(out_path: str, section: str, items: list) -> tuple[bool, str]:
    """Validate each item against its sub-model, then append the normalized batch
    (atomic: on any invalid item, append nothing). Shared by append_section and
    table_to_entities."""
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
    rec[section].extend(normalized)
    pathlib.Path(_partial_path(out_path)).write_text(json.dumps(rec))
    return True, f"appended {len(normalized)} to {section}; total now {len(rec[section])}"


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
    return True, f"counts={counts} meta={meta}"


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
        if e.outcome not in ("immunogenic", "positive"):
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


# class-II coverage anchor (v2.12 P23): the corpus repeatedly states a class-II / CD4 restriction in
# prose or names a DR/DP/DQ allele yet mints ZERO class-II epitopes (the Keskin 30568305 regression:
# GPC1/SHANK2/SVEP1 are CD4/class-II hits that landed as 0 class-II epitopes). TWO sides, both -> the
# record is suspect and routes to needs_review:
#   (a) a CD4 cue appears anywhere (evidence/epitope quoted_text or a cd4 subset) but NO class-II-typed
#       record exists (no MinimalEpitope mhc_class=='II' AND no evidence mhc_class=='class_ii'); and
#   (b) a named DR/DP/DQ (or I-A/I-E) restriction / 'class II-restricted' is quoted but 0 class-II
#       epitopes were minted.
# SOFT/block-once, overridable with allow_missing_class_ii (a HARD override -> needs_review, conservative).
_CD4_CUE = re.compile(r"\bcd4\b|helper|class[ \-]?ii|hla-?d|\bdr\b|\bdp\b|\bdq\b|h-2i|i-[ae]", re.I)
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
    cd4_subset = any(getattr(e, "t_cell_subset", None) == "cd4" for e in rec.evidence)
    has_cd4_cue = cd4_subset or any(_CD4_CUE.search(t) for t in texts)
    has_named = any(_NAMED_CII_ALLELE.search(t) for t in texts) or any(
        _NAMED_CII_ALLELE.search(e.hla_allele or "") for e in rec.epitopes if getattr(e, "hla_allele", None)
    )
    msgs = []
    if has_cd4_cue:
        msgs.append(
            "a CD4 / class-II cue appears in the text but ZERO class-II records were minted "
            "(no MinimalEpitope mhc_class='II', no evidence mhc_class='class_ii') -- a CD4 response is "
            "class-II-restricted; mint the class-II epitope(s) and/or type the evidence row(s)")
    if has_named:
        msgs.append(
            "a named DR/DP/DQ (or I-A/I-E) restriction is quoted but ZERO class-II epitopes were minted "
            "-- add the class-II MinimalEpitope for the quoted restriction")
    return msgs


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
                     allow_ungrounded_safety_grade: bool = False) -> tuple[bool, str]:
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
    # soft class-II minting nudge (v2.12 P23, block-once, overridable): a CD4/class-II cue or a named
    # DR/DP/DQ restriction is present but no class-II record was minted (the Keskin regression).
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
    ) if allowed and fired]
    if overrides_used:
        full = full.model_copy(update={"finalize_overrides_used": overrides_used})
        pathlib.Path(out_path).write_text(full.model_dump_json(indent=2, exclude_none=True))
    pathlib.Path(_partial_path(out_path)).unlink(missing_ok=True)  # success: clean up partial
    if n_epi_merged:
        gmsg += f" | canonicalized epitopes: merged {n_epi_merged} duplicate record(s)"
    if n_scr_dropped:
        gmsg += f" | screening: dropped {n_scr_dropped} row(s) superseded by named evidence"
    if overrides_used:
        gmsg += f" | OVERRIDES USED -> needs_review: {overrides_used}"
    return True, gmsg


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
    targets = list(sheets) if sheets else [sheet]   # multi-sheet mode iff `sheets` given
    multi = bool(sheets)
    items: list = []
    per_sheet: list[str] = []
    for sh in targets:
        try:
            headers, rows = _read_xlsx_dicts(path, sh, header_row)
        except Exception as e:
            return False, f"ERROR reading sheet {sh!r} in {path!r}: {e}"
        if multi:                                   # inject the reserved sheet-name column
            for r in rows:
                r["__sheet__"] = sh
            headers = list(headers) + ["__sheet__"]
        miss = table_map.missing_columns(mapping["fields"], flt, headers)
        if miss:
            return False, f"column(s) {miss} not found in sheet {sh!r}; headers: {headers}"
        if flt:
            ok, kept = table_map.apply_filter(rows, flt)
            if not ok:
                return False, kept  # error message
            rows = kept
        try:
            sheet_items = [table_map.apply_mapping(r, mapping["fields"]) for r in rows]
        except ValueError as e:
            return False, str(e)
        items.extend(sheet_items)
        per_sheet.append(f"{sh}:{len(sheet_items)}")
    if not items:
        return False, f"0 rows after filter; nothing added to {section}"
    ok, msg = _validate_and_append(out_path, section, items)
    if ok and multi:
        msg += f" | from {len(targets)} sheets [{', '.join(per_sheet)}]"
    return ok, msg
