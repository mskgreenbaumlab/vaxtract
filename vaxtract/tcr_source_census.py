"""SOURCE CENSUS — what TCR tables did this paper's supplements actually OFFER?

Motivation (2026-08-24): every finalize gate in `agent_core` checks a record against ITSELF
(provenance present, ids resolve, counts reconcile). None of them can see the source, so a record
that captured a FRACTION of what the supplement held still finalizes clean. Observed live on Hu
33479501: the designated keep-file asserts `tcr_seq_status='both'` with 3 TCR-seq methods and ZERO
clonotypes, while its 257MB supplement holds 10 loadable tetramer-track sheets (269 clonotypes /
654 timepoint observations). A 3-run replicate routed 11 / 6 / 1 sheets to the loaders on identical
input; all three finalized clean. The loaders are deterministic GIVEN a sheet -- WHICH sheets get
routed to them is the non-deterministic step, and that is what this module measures.

This module is a PURE FUNCTION of the paper directory: dir -> the loadable TCR sheets in it. It is
what `agent_core._tcr_source_coverage_gap` compares the record's load-ledger against.

THREE BUCKETS (a sheet is exactly one):
  * REQUIRED       a loadable TCR shape that curation policy says belongs in the record.
  * OUT_OF_SCOPE   a loadable TCR shape that curation policy EXCLUDES (see `OUT_OF_SCOPE_RULES`).
                   Never required -> never blocks. Loading one anyway is not a coverage error.
  * NOT_LOADABLE   everything else (no deterministic loader can read it).

SHAPE, NOT NAME. Sheet names are paper-specific and typo-ridden ('Suppl 8 Pt 5Tet-spec TCRa track'
has no space after the patient); a name-pattern census would be brittle in both directions. So a
sheet is classified by probing its header rows and reusing the SAME shape knowledge the adapters
encode -- `tcr_track_adapter._infer_locus` / `._split_identity` for the wide tet-spec track shape,
`tcr_reactive_adapter._col` / `._aa` for the paired single-cell reactive shape. If the adapters
learn a new shape, this census inherits it.

SPEED. Hu's supplement is 257MB / 50 sheets with a 159MB shared-string table; `openpyxl` (even
read_only) parses that table on open (~13s) before a single cell is read, and the census must run at
EVERY finalize. So this module never touches openpyxl: sheet names come from `xl/workbook.xml` (a
1KB zip member), the first ~3 rows come from the head of each sheet's XML, and only the handful of
shared-string indices those rows reference are resolved, by one forward scan that stops at the
largest index needed. Full 50-sheet census of Hu's 257MB workbook: ~0.6s.

The header row is index 1 (0-based) for BOTH shapes -- that is the default the wired tools use
(`add_tcr_track` header_row_index=1; `add_reactive_tcr` has no header knob and its adapter defaults
to header_row=1). A sheet whose header sits elsewhere is deliberately NOT counted as required: the
census must not demand a load the agent's tools cannot perform with their defaults.
"""
from __future__ import annotations

import dataclasses
import pathlib
import re
import zipfile

from . import tcr_reactive_adapter as _tra
from . import tcr_track_adapter as _tta

# ---------------------------------------------------------------------------------------------
# curation policy: which loadable TCR shapes are OUT OF SCOPE
# ---------------------------------------------------------------------------------------------
# POLICY, ENCODED AS DATA -- signed 2026-08-25
# (docs/POLICY_tcr_sheet_scope_2026-08-24.md; docs/RECOMMENDATIONS_2026-08-25_hu_gatecols_and_keskin.md).
# `tcr_reactive_adapter.KESKIN_SHEETS` carries the comment "neoantigen-reactive only; MOESM10 (a)
# 'Tumor-associated' is excluded": Keskin's '(a)Pt7 Tumor-associated Tcells ' is a perfectly valid
# reactive SHAPE (same TRAV/TRAJ/CDR3 columns, 228 clonotypes) but is NOT antigen-sorted -- it is the
# tumour-infiltrating repertoire, so its clonotypes carry no neoantigen specificity. If the census
# called it required, this gate would COERCE that 228-row contamination into every Keskin record.
# The rule is therefore "a loadable TCR sheet whose sheet name or title cell says the cells were not
# antigen-sorted is out of scope", expressed as patterns over (sheet name + row-0 title cell) rather
# than as a hardcoded (file, sheet) pair, so a second paper's tumour-associated table is excluded for
# the same stated reason instead of silently becoming required.
# DIRECTION OF ERROR: an out-of-scope classification can only make the gate QUIETER (a mis-excluded
# sheet is a missed block, never a false block), so this list is the safe place for policy to live.
OUT_OF_SCOPE_RULES: list[tuple[str, "re.Pattern[str]", str]] = [
    ("tumour_associated",
     re.compile(r"tumou?r[\s_\-]*(associated|infiltrat\w*)|\bTILs?\b", re.I),
     "tumour-associated / tumour-infiltrating T cells are not antigen-sorted, so their clonotypes "
     "carry no neoantigen specificity (tcr_reactive_adapter.KESKIN_SHEETS excludes this sheet class "
     "by the same rule)"),
    ("unsorted",
     re.compile(r"\bunsorted\b|\bnot\s+sorted\b|\bbulk\s+(RNA|repertoire|sequencing|T\s*cells?)\b",
                re.I),
     "an unsorted / bulk repertoire table is not an antigen-specific clonotype set"),
]

REQUIRED = "required"
OUT_OF_SCOPE = "out_of_scope"
NOT_LOADABLE = "not_loadable"

_LOADER_TOOL = {"track": "add_tcr_track", "reactive": "add_reactive_tcr",
                "reactivity": "add_table (per_cell)"}
HEADER_ROW_INDEX = 1          # 0-based; what both wired loader tools use by default
# Rows probed per sheet. The TCR shapes need 3 (title, header, first data row); a reactivity matrix
# carries a two- or three-row merged header, so its assay labels sit at index 3. Probing more rows
# is nearly free -- `_probe_rows` already reads a 1MB head of the member either way.
_PROBE_ROWS = 6


@dataclasses.dataclass(frozen=True)
class CensusSheet:
    """One supplementary sheet, classified. `file` is the BASENAME (the ledger matches on it, so a
    relative-vs-absolute path in a tool call never reads as a different sheet)."""
    file: str
    sheet: str
    shape: str                 # 'track' | 'reactive' | 'reactivity'
    scope: str                 # REQUIRED | OUT_OF_SCOPE
    loader: str                # the MCP tool that reads this shape
    reason: str
    path: str = ""             # full path, for the block message / a direct load
    n_expected: int = 0        # reactivity only: positive assay cells the sheet holds
    # Reactivity STRUCTURE. The gate needs only `n_expected`, but the evidence ruler
    # (`stability_lab/evidence_source_gold.py`) needs to read the same sheet the same way, and a
    # second hand-written column map would be free to drift from this one. One description, both
    # consumers -- so a gate that says '279 positive calls' and a gold that scores against them
    # cannot disagree about which columns those calls are in.
    data_start: int = 0                          # first data row (headers span 2-3 merged rows)
    assay_cols: tuple = ()                       # the scored 1/0/n.d. columns
    imp_col: int | None = None                   # immunizing-peptide SEQUENCE column
    assay_pep_col: int | None = None             # assay-peptide SEQUENCE column, when distinct
    id_col: int | None = None                    # '<patient>-ASP4'-style id, for a blank patient cell
    patient_col: int | None = None

    @property
    def key(self) -> tuple[str, str]:
        return (self.file.lower(), self.sheet)


@dataclasses.dataclass(frozen=True)
class SourceCensus:
    required: list[CensusSheet]
    out_of_scope: list[CensusSheet]
    not_loadable: list[tuple[str, str]]
    files_scanned: int = 0
    files_failed: list[str] = dataclasses.field(default_factory=list)
    # Reactivity matrices are NOT a TCR shape and deliberately do NOT join `required`: the TCR gate
    # demands a TCR loader, and coercing these into it would produce a wrong block message. They get
    # their own list and their own gate (`agent_core._reactivity_source_coverage_gap`).
    reactivity: list[CensusSheet] = dataclasses.field(default_factory=list)


# ---------------------------------------------------------------------------------------------
# minimal xlsx reader (no openpyxl: see the module docstring's SPEED note)
# ---------------------------------------------------------------------------------------------
# Attribute ORDER varies by writer (Excel emits Id before Target, openpyxl the reverse), so the
# element is matched first and each attribute read out of it independently.
_SHEET_EL_RE = re.compile(r"<sheet\b[^>]*/?>")
_REL_EL_RE = re.compile(r"<Relationship\b[^>]*/?>")
_ATTR_RE = re.compile(r'([\w:]+)\s*=\s*"([^"]*)"')
_ROW_OPEN_RE = re.compile(rb"<row\b[^>]*?/?>")
_ROWNUM_RE = re.compile(rb'r="(\d+)"')
_CELL_RE = re.compile(rb"<c\b([^>]*)/>|<c\b([^>]*)>(.*?)</c>", re.S)
_REF_RE = re.compile(rb'r="([A-Z]+)\d+"')
_TYPE_RE = re.compile(rb't="([a-zA-Z]+)"')
_V_RE = re.compile(rb"<v[^>]*>(.*?)</v>", re.S)
_T_RE = re.compile(rb"<t[^>]*>(.*?)</t>", re.S)
_ENTITIES = (("&lt;", "<"), ("&gt;", ">"), ("&quot;", '"'), ("&apos;", "'"), ("&amp;", "&"))
_NUMREF_RE = re.compile(r"&#(x)?([0-9A-Fa-f]+);")

_SHEET_HEAD_BYTES = 1 << 20   # enough for the first rows of any real header block
_SS_CHUNK = 1 << 22


def _unescape(s: str) -> str:
    if "&" not in s:
        return s
    for ent, ch in _ENTITIES:
        s = s.replace(ent, ch)
    # numeric character references: openpyxl writes 'CDR3&#946;' where Excel writes the literal 'β',
    # and the chain-naming header this census keys on is exactly such a character.
    return _NUMREF_RE.sub(lambda m: chr(int(m.group(2), 16 if m.group(1) else 10)), s)


def _col_index(ref: bytes) -> int:
    n = 0
    for ch in ref.decode("ascii", "replace"):
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def _sheet_index(zf: zipfile.ZipFile) -> list[tuple[str, str]]:
    """[(sheet name, zip member)] in workbook order, straight out of xl/workbook.xml."""
    wb = zf.read("xl/workbook.xml").decode("utf8", "replace")
    rels = {}
    for el in _REL_EL_RE.findall(zf.read("xl/_rels/workbook.xml.rels").decode("utf8", "replace")):
        attrs = dict(_ATTR_RE.findall(el))
        if attrs.get("Id") and attrs.get("Target"):
            rels[attrs["Id"]] = attrs["Target"]
    out = []
    for el in _SHEET_EL_RE.findall(wb):
        attrs = dict(_ATTR_RE.findall(el))
        name, rid = attrs.get("name"), attrs.get("r:id")
        if not name or not rid:
            continue
        tgt = rels.get(rid)
        if not tgt:
            continue
        tgt = tgt.lstrip("/")
        if not tgt.startswith("xl/"):
            tgt = "xl/" + tgt
        out.append((_unescape(name), tgt))
    return out


def _probe_rows(zf: zipfile.ZipFile, member: str, nrows: int) -> list[dict[int, object]]:
    """First `nrows` rows of one sheet as {column index: value}. Values are either a float, a str,
    or ('S', idx) for a not-yet-resolved shared string. Reads only the HEAD of the member.

    Rows are placed by their `r=` NUMBER, and a skipped row leaves an empty slot, because that is
    what openpyxl's `iter_rows` hands the adapters -- a sheet whose row 1 is absent from the XML must
    still put its header at index 1, or the census would probe the wrong row. Row starts are found
    with a single linear scan (a `<row>...</row>` regex backtracks across the whole buffer for every
    self-closing empty row, which is quadratic on exactly the sparse sheets this reads)."""
    try:
        with zf.open(member) as fh:
            head = fh.read(_SHEET_HEAD_BYTES)
    except (KeyError, OSError, zipfile.BadZipFile):
        return []
    starts = [(m.start(), m.end(), m.group()) for m in _ROW_OPEN_RE.finditer(head)]
    by_index: dict[int, dict[int, object]] = {}
    for pos, (start, end, tag) in enumerate(starts):
        rn = _ROWNUM_RE.search(tag)
        r_index = (int(rn.group(1)) - 1) if rn else pos
        if r_index >= nrows:
            break
        if tag.endswith(b"/>"):                       # <row r="3"/> -- a row that exists but is empty
            by_index.setdefault(r_index, {})
            continue
        stop = starts[pos + 1][0] if pos + 1 < len(starts) else len(head)
        body_bytes = head[end:stop]
        row: dict[int, object] = {}
        for i, cm in enumerate(_CELL_RE.finditer(body_bytes)):
            attrs = cm.group(1) or cm.group(2) or b""
            body = cm.group(3) or b""
            ref = _REF_RE.search(attrs)
            col = _col_index(ref.group(1)) if ref else i
            tm = _TYPE_RE.search(attrs)
            ctype = tm.group(1).decode() if tm else "n"
            vm = _V_RE.search(body)
            if ctype == "s":
                if vm:
                    try:
                        row[col] = ("S", int(vm.group(1)))
                    except ValueError:
                        pass
            elif ctype == "inlineStr":
                tt = _T_RE.search(body)
                if tt:
                    row[col] = _unescape(tt.group(1).decode("utf8", "replace"))
            elif vm:
                raw = vm.group(1).decode("utf8", "replace")
                if ctype == "str":
                    row[col] = _unescape(raw)
                else:
                    try:
                        row[col] = float(raw)
                    except ValueError:
                        row[col] = raw
        by_index[r_index] = row
    if not by_index:
        return []
    return [by_index.get(i, {}) for i in range(max(by_index) + 1)]


def _parse_row_cells(body_bytes: bytes) -> dict[int, object]:
    """One row's `<c>` cells -> {column index: float | str | ('S', shared-string index)}."""
    row: dict[int, object] = {}
    for i, cm in enumerate(_CELL_RE.finditer(body_bytes)):
        attrs = cm.group(1) or cm.group(2) or b""
        body = cm.group(3) or b""
        ref = _REF_RE.search(attrs)
        col = _col_index(ref.group(1)) if ref else i
        tm = _TYPE_RE.search(attrs)
        ctype = tm.group(1).decode() if tm else "n"
        vm = _V_RE.search(body)
        if ctype == "s":
            if vm:
                try:
                    row[col] = ("S", int(vm.group(1)))
                except ValueError:
                    pass
        elif ctype == "inlineStr":
            tt = _T_RE.search(body)
            if tt:
                row[col] = _unescape(tt.group(1).decode("utf8", "replace"))
        elif vm:
            raw = vm.group(1).decode("utf8", "replace")
            if ctype == "str":
                row[col] = _unescape(raw)
            else:
                try:
                    row[col] = float(raw)
                except ValueError:
                    row[col] = raw
    return row


def _probe_rows_full(zf: zipfile.ZipFile, member: str) -> list[dict[int, object]]:
    """EVERY row of one sheet, same value encoding as `_probe_rows`. Only ever called on a member
    already checked against `_MAX_SCAN_BYTES`, so reading it whole is bounded."""
    try:
        with zf.open(member) as fh:
            blob = fh.read()
    except (KeyError, OSError, zipfile.BadZipFile):
        return []
    starts = [(m.start(), m.end(), m.group()) for m in _ROW_OPEN_RE.finditer(blob)]
    by_index: dict[int, dict[int, object]] = {}
    for pos, (start, end, tag) in enumerate(starts):
        rn = _ROWNUM_RE.search(tag)
        r_index = (int(rn.group(1)) - 1) if rn else pos
        if tag.endswith(b"/>"):
            by_index.setdefault(r_index, {})
            continue
        stop = starts[pos + 1][0] if pos + 1 < len(starts) else len(blob)
        by_index[r_index] = _parse_row_cells(blob[end:stop])
    if not by_index:
        return []
    return [by_index.get(i, {}) for i in range(max(by_index) + 1)]


def _resolve_shared_strings(zf: zipfile.ZipFile, needed: set[int]) -> dict[int, str]:
    """Resolve just the shared-string indices the probed rows referenced.

    One forward scan over xl/sharedStrings.xml that STOPS at the largest index needed and only
    materialises the entries asked for -- Hu's table is 159MB / 4.2M entries and openpyxl parses all
    of it on open; this scans the bytes (~0.3s) and builds a dict of ~40 strings."""
    if not needed:
        return {}
    want = sorted(needed)
    out: dict[int, str] = {}
    ptr = 0
    idx = 0            # global index of the next '<si' occurrence not yet consumed
    carry = b""
    try:
        fh = zf.open("xl/sharedStrings.xml")
    except (KeyError, OSError, zipfile.BadZipFile):
        return {}
    with fh:
        while ptr < len(want):
            chunk = fh.read(_SS_CHUNK)
            if not chunk:
                break
            buf = carry + chunk
            pos, start, stalled = 0, -1, False
            while ptr < len(want):
                start = buf.find(b"<si", pos)
                if start < 0:
                    break
                if idx == want[ptr]:
                    end = buf.find(b"</si>", start)
                    if end < 0:                 # entry straddles the chunk boundary -> read more
                        stalled = True
                        break
                    out[idx] = _unescape("".join(
                        t.decode("utf8", "replace") for t in _T_RE.findall(buf[start:end])))
                    ptr += 1
                idx += 1
                pos = start + 3
            if stalled:
                carry = buf[start:]             # keep the incomplete entry, do NOT count it
            else:
                carry = buf[max(pos, len(buf) - 2):]   # keep at most a split '<s' tag
    return out


def _probe_workbook(path: pathlib.Path, nrows: int = 3) -> list[tuple[str, list[list]]]:
    """[(sheet name, rows)] with every shared string resolved. rows are dense lists (index = column),
    matching what the adapters expect from openpyxl's `values_only` iteration."""
    with zipfile.ZipFile(path) as zf:
        sheets = _sheet_index(zf)
        probed = [(name, _probe_rows(zf, member, nrows)) for name, member in sheets]
        needed = {v[1] for _, rows in probed for row in rows for v in row.values()
                  if isinstance(v, tuple)}
        strings = _resolve_shared_strings(zf, needed)
    out = []
    for name, rows in probed:
        dense = []
        for row in rows:
            width = (max(row) + 1) if row else 0
            cells: list = [None] * width
            for col, val in row.items():
                cells[col] = strings.get(val[1]) if isinstance(val, tuple) else val
            dense.append(cells)
        out.append((name, dense))
    return out


# ---------------------------------------------------------------------------------------------
# shape classification (reuses the adapters' own shape knowledge)
# ---------------------------------------------------------------------------------------------
# The track loader keeps a row only when the paper's own specificity flag (col 1) is 1. That column
# is what separates a tet-spec TRACK sheet from the paper's UNSORTED bulk repertoire sheet, which
# has the same identity column ('CDR3β') and would otherwise probe identically.
_FLAG_HEADER = re.compile(r"tet(ramer)?[\s\-]*spec|1\s*=\s*yes|specific\s*\?", re.I)


def _text(row: list, col: int) -> str:
    if col < len(row) and isinstance(row[col], str):
        return row[col].strip()
    return ""


def _is_track_sheet(rows: list[list]) -> str:
    """'' when this is not the wide tet-spec track shape, else why it is.

    Mirrors `tcr_track_adapter.parse_tetraspec_track`'s own preconditions: an identity column whose
    header names the chain (`_infer_locus`, which REFUSES to guess), the paper's specificity flag in
    col 1, at least one timepoint label from col 2, and a first data row whose identity cell parses
    to a CDR3 (`_split_identity`)."""
    if len(rows) <= HEADER_ROW_INDEX:
        return ""
    header = rows[HEADER_ROW_INDEX]
    above = rows[HEADER_ROW_INDEX - 1] if HEADER_ROW_INDEX else []
    locus = _tta._infer_locus(header[0] if header else None)
    if locus not in ("TRA", "TRB"):
        return ""
    if not (_FLAG_HEADER.search(_text(header, 1)) or _FLAG_HEADER.search(_text(above, 1))):
        return ""
    labels = [_text(header, c) for c in range(2, len(header))]
    if not any(labels):
        return ""
    data = rows[HEADER_ROW_INDEX + 1] if len(rows) > HEADER_ROW_INDEX + 1 else []
    if not data or not isinstance(data[0], str) or not _tta._split_identity(data[0])[0]:
        return ""
    return (f"wide tet-spec track sheet: {locus} identity column {header[0]!r}, specificity flag in "
            f"col 1, {sum(1 for l in labels if l)} timepoint column(s)")


def _is_reactive_sheet(rows: list[list]) -> str:
    """'' when this is not the paired single-cell reactive shape, else why it is.

    Mirrors `tcr_reactive_adapter.load_reactive_clonotypes`: it locates its columns with `_col` on
    header row 1 and keeps a row when `_aa` reads a CDR3 out of the alpha or beta column.

    BOTH chain columns are required, because that adapter's clonotype IS the unique (CDR3a, CDR3b)
    PAIR and its `cell_count` counts one row per CELL. A single-chain table (Keskin MOESM8 'Patient 8
    CDR3 sequences from analysis of bulk RNA', TRAV/TRAJ/'Alpha CDR3 amino acid seq'/UMI_count) has
    the same column names but one row per CLONOTYPE, so feeding it to that adapter would mint
    thousands of unsorted pseudo-clonotypes whose `cell_count` means nothing. It is a different
    shape, and the census must not demand it."""
    if len(rows) <= HEADER_ROW_INDEX:
        return ""
    header = rows[HEADER_ROW_INDEX]
    ca, cb = _tra._col(header, "alpha cdr3 amino"), _tra._col(header, "beta cdr3 amino")
    if ca is None or cb is None:
        return ""
    if _tra._col(header, "trav") is None or _tra._col(header, "trbv") is None:
        return ""
    data = rows[HEADER_ROW_INDEX + 1] if len(rows) > HEADER_ROW_INDEX + 1 else []
    if not data or not (_tra._aa(_tra._cell(data, ca)) or _tra._aa(_tra._cell(data, cb))):
        return ""
    return ("paired single-cell reactive sheet: alpha+beta CDR3 amino-acid columns with V/J calls, "
            "one row per cell")


# ---------------------------------------------------------------------------------------------
# reactivity matrices (the evidence lane's denominator)
# ---------------------------------------------------------------------------------------------
# A reactivity matrix is one row per assayed peptide and several ASSAY columns across, each cell a
# scored call. Hu 33479501 'Suppl Dataset 4b' has seven (two timepoints x ex-vivo / pre-stimulated /
# minigene / minigene-blocked / autologous-tumour). Counting its '1' cells is the evidence lane's
# denominator -- the thing `CHECKPOINT_2026-08-24_evidence_swing.md` believed did not exist.
#
# Detection is TWO-PHASE, because the shape's proof lives in the DATA, not the header:
#   1. a cheap header prefilter over the probed rows (this must be permissive -- it only decides
#      which sheets are worth opening);
#   2. a bounded full scan that finds columns whose entire value domain is the scored vocabulary.
#      A sheet with no such column is not a reactivity matrix, so phase 2 is also the confirmation.
_REACTIVITY_HEADER = re.compile(r"reactivit|immunogenic|\bresponse\b|\breactive\b", re.I)
# Scored-call vocabulary. A column qualifies only if EVERY non-blank value is one of these, at least
# one is a positive, and at least one is NOT -- three conditions that each kill a real false
# positive seen in the corpus:
#   * whole domain scored   -> keeps out an id/rank column (values 1,2,3,...);
#   * has a positive        -> an all-negative column is not evidence of anything;
#   * has a non-positive    -> keeps out a NORMALISATION CONSTANT. PMID 40790272 'Ext. Fig 6a' has a
#     'Baseline FC' column that is 1.0 on every row next to a 'Max response FC' header; without this
#     condition the census called it a 10-cell reactivity matrix and the gate would have demanded a
#     load of it. A column that never varies discriminates nothing.
_POSITIVE_CALL = {"1", "1.0", "+", "pos", "positive", "yes"}
_NEGATIVE_CALL = {"0", "0.0", "-", "neg", "negative", "no"}
_UNKNOWN_CALL = {"n.d.", "nd", "n.d", "not determined", "n/a", "na", "nt", "n.t."}
_SCORED = _POSITIVE_CALL | _NEGATIVE_CALL | _UNKNOWN_CALL
# Below this a sheet is not worth gating on (a stray 0/1 column is not an assay matrix).
MIN_REACTIVITY_CELLS = 10
# SCOPE, DECIDED 2026-08-24 (maintainer): Hu 33479501's 'Supp11b NetMHCpan2.4 CD4 T cell'
# (Supplementary Dataset 11, 11 positive calls) IS in scope and the gate is right to require it.
# Recorded because the sheet is easy to mistake for a prediction table that merely duplicates
# Dataset 4b -- it is not. Its scored columns are IFN-g ELISPOT readouts against NetMHCpan-predicted
# class-II peptides, i.e. measurements the record should carry, and NO run has ever cited it. Do not
# add an OUT_OF_SCOPE rule for it; a Hu record that omits those 11 calls is genuinely short.
# Full-scanning a 185MB sheet XML at every finalize is not acceptable; real reactivity matrices are
# tiny (Hu's 4a/4b are 149KB and 258KB) while the sheets that would hurt are the huge TCR ones.
_MAX_SCAN_BYTES = 8 << 20


def _combined_header_labels(rows: list, data_start: int, resolve) -> dict:
    """{column index: its full header label} for a header block of ANY depth.

    Reactivity sheets merge their headers over two or three rows -- Hu 4b is group / field / assay
    ('Assay peptide' | 'Sequence'), and its assay names sit on the THIRD row. Each header row is
    forward-filled across the blanks a merged cell leaves (starting only at its first non-empty
    cell, so a row that begins blank does not smear), then the rows are joined per column. This is
    the N-row generalisation of `epitope_manifest_adapter._combined_labels`, which folds at most two
    and would therefore read 4b's assay-label row as data."""
    width = max((max(r) + 1) for r in rows[:data_start] if r) if any(rows[:data_start]) else 0
    parts: dict = {}
    for row in rows[:data_start]:
        # Skip the TITLE row. A one-cell row spanning the sheet ('Supplementary Dataset 4a. Summary
        # of immunizing peptides and CD8+ T cell responses...') would forward-fill its whole
        # sentence into every column, and any word in it would then appear to name every column --
        # 'immunizing' in that title made column 6 look like the immunizing-peptide column.
        if sum(1 for v in row.values() if isinstance(resolve(v), str) and resolve(v).strip()) < 2:
            continue
        current = ""
        for col in range(width):
            val = resolve(row.get(col))
            text = val.strip() if isinstance(val, str) else ""
            if text:
                current = text
            elif not current:
                continue                       # nothing to fill forward from yet
            if current:
                parts.setdefault(col, []).append(current)
    out = {}
    for col, seq in parts.items():
        seen, joined = set(), []
        for t in seq:                          # a group label repeated by the fill is not new text
            if t not in seen:
                seen.add(t)
                joined.append(t)
        out[col] = " ".join(joined)
    return out


_WT_LABEL = re.compile(r"\bwild[\s\-]?type\b|\bwt\b", re.I)
_ID_VALUE = re.compile(r"^\d+\s*-\s*\S")       # '1-ASP4', '11-IMP02' -- the patient number prefixes it


def _describe_reactivity(rows: list, data_start: int, assay_cols: list, resolve) -> dict:
    """Which column is the immunizing peptide, the assay peptide, the id, the patient.

    Peptide columns are found by DATA (amino-acid strings, via the manifest adapter's own peptide
    test) and then named by their HEADER, so a wild-type control column is never mistaken for a
    target: crediting a record for naming the WT peptide would be a false positive the gold could
    not detect."""
    try:
        from . import epitope_manifest_adapter as _ema
        is_peptide = lambda s: bool(_ema._PEPTIDE.match(s))          # noqa: E731
    except Exception:                                                # noqa: BLE001
        is_peptide = lambda s: bool(re.match(r"^[ACDEFGHIKLMNPQRSTVWY]{7,}$", s))  # noqa: E731
    labels = _combined_header_labels(rows, data_start, resolve)
    data = rows[data_start:]
    width = max((max(r) + 1) for r in data if r) if any(data) else 0

    def column(col):
        return [resolve(r.get(col)) for r in data]

    peptide_cols, id_cols = [], []
    for col in range(width):
        if col in assay_cols:
            continue
        vals = [v.strip() for v in column(col) if isinstance(v, str) and v.strip()]
        if not vals:
            continue
        if sum(1 for v in vals if is_peptide(v.upper())) >= max(2, len(vals) // 2):
            peptide_cols.append(col)
        elif sum(1 for v in vals if _ID_VALUE.match(v)) >= max(2, len(vals) // 2):
            id_cols.append((-len(vals), col))       # most-populated first
    # The id column exists to recover a BLANK patient cell, so the one that covers the most rows is
    # the only useful choice. On Hu 4b the immunizing-peptide id column is blank on exactly the 11
    # parentless Pt11/Pt12 rows -- the rows that need recovering -- while the assay-peptide id
    # column is populated on all 375.
    id_cols = [c for _, c in sorted(id_cols)]

    def labelled(cols, *words):
        for c in cols:
            lab = labels.get(c, "").lower()
            if all(w in lab for w in words):
                return c
        return None

    targets = [c for c in peptide_cols if not _WT_LABEL.search(labels.get(c, ""))]
    imp = labelled(targets, "immunizing")
    assay_pep = labelled(targets, "assay") or labelled(targets, "mutated")

    def typical_len(c):
        vals = [v.strip() for v in column(c) if isinstance(v, str) and v.strip()]
        return max((len(v) for v in vals), default=0)

    if imp is None and targets:
        # No sheet-specific wording: the immunizing peptide is the LONG one -- epitopes and assay
        # peptides are nested inside it, so the longest target column is it by construction.
        imp = max(targets, key=typical_len)
    # Hu 4a labels the 9-mer 'Mutated peptide Sequence', not 'assay'. If a second
    # target column is strictly shorter than the immunizing peptide, it is the
    # assayed epitope. A single-column sheet (Hu Supp11b) stays assay_pep=None.
    if assay_pep is None and imp is not None:
        others = [c for c in targets if c != imp]
        if others:
            shorter = min(others, key=typical_len)
            if typical_len(shorter) < typical_len(imp):
                assay_pep = shorter
    patient = labelled(range(width), "patient")
    return {"imp_col": imp,
            "assay_pep_col": assay_pep if assay_pep != imp else None,
            "id_col": id_cols[0] if id_cols else None,
            "patient_col": patient if patient is not None else (0 if width else None)}


def _norm_call(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        return str(int(v)) if v.is_integer() else str(v)
    return str(v).strip().lower()


def _looks_like_reactivity(rows: list[list]) -> str:
    """Cheap header prefilter: '' or why this sheet is worth a full scan. Permissive by design."""
    for row in rows[:HEADER_ROW_INDEX + 3]:
        for cell in row:
            if isinstance(cell, str) and _REACTIVITY_HEADER.search(cell):
                return cell.strip()[:80]
    return ""


def _scan_scored_columns(zf: zipfile.ZipFile, member: str, with_structure: bool) -> dict:
    """The sheet's reactivity structure, or {} when it is not a reactivity matrix.

    Returns positives / n_cols / n_rows / data_start / assay_cols plus the peptide, id and patient
    columns -- everything both consumers need, derived once.

    A column counts only when every non-blank value it holds is a scored call and at least one is
    positive. The data block starts at the first row where two or more columns already hold scored
    calls -- these sheets carry two- and three-row merged headers, and including a header row would
    poison every column's domain with its label."""
    try:
        if zf.getinfo(member).file_size > _MAX_SCAN_BYTES:
            return {}
    except KeyError:
        return {}
    rows = _probe_rows_full(zf, member)
    if not rows:
        return {}
    start = None
    for i, row in enumerate(rows):
        if sum(1 for v in row.values() if _norm_call(v) in _SCORED) >= 2:
            start = i
            break
    if start is None:
        return {}
    # A column's unknown marker ('n.d.') is a SHARED string. Resolve only the indices sitting in
    # columns that already hold a numeric call -- typically one or two distinct strings for the
    # whole sheet -- rather than the workbook's string table, which is 159MB in Hu's supplement.
    raw: dict[int, set] = {}
    for row in rows[start:]:
        for col, val in row.items():
            if val is not None and val != "":
                raw.setdefault(col, set()).add(val if not isinstance(val, tuple) else val)
    candidates = {c for c, dom in raw.items()
                  if any(_norm_call(v) in _POSITIVE_CALL for v in dom if not isinstance(v, tuple))}
    needed = {v[1] for c in candidates for v in raw[c] if isinstance(v, tuple)}
    strings = _resolve_shared_strings(zf, needed) if needed else {}

    def norm(v):
        return _norm_call(strings.get(v[1], "\x00unresolved") if isinstance(v, tuple) else v)

    good = []
    for c in candidates:
        dom = {norm(v) for v in raw[c]}
        if dom <= _SCORED and (dom & _POSITIVE_CALL) and (dom - _POSITIVE_CALL):
            good.append(c)
    if not good:
        return {}
    good.sort()
    positives = sum(1 for row in rows[start:] for c in good
                    if norm(row.get(c)) in _POSITIVE_CALL)
    # Structure is what the evidence RULER needs; the finalize gate needs only the count. Deriving
    # it means resolving the sheet's non-assay shared strings (header labels and every peptide
    # sequence), which forces a scan deep into a 159MB string table -- ~2s per workbook on Hu, paid
    # at every finalize for nothing. So it is opt-in, and skipping it cannot change `positives`.
    structure = {"imp_col": None, "assay_pep_col": None, "id_col": None, "patient_col": None}
    if with_structure:
        text_ids = {v[1] for row in rows for col, v in row.items()
                    if isinstance(v, tuple) and col not in set(good)}
        text = _resolve_shared_strings(zf, text_ids) if text_ids else {}

        def resolve(v, _t=text, _s=strings):
            if isinstance(v, tuple):
                return _t.get(v[1], _s.get(v[1], ""))
            return v

        structure = _describe_reactivity(rows, start, good, resolve)
    return {"positives": positives, "n_cols": len(good), "n_rows": len(rows) - start,
            "data_start": start, "assay_cols": tuple(good), **structure}


def _out_of_scope_reason(sheet_name: str, rows: list[list]) -> str:
    """The curation exclusion, matched against the sheet NAME and its title cell (row 0, col 0) --
    the two places a supplement says what a table's cells were sorted on."""
    title = _text(rows[0], 0) if rows else ""
    hay = f"{sheet_name} {title}"
    for _, pattern, why in OUT_OF_SCOPE_RULES:
        if pattern.search(hay):
            return why
    return ""


def classify_sheet(sheet_name: str, rows: list[list]) -> tuple[str, str, str]:
    """(shape, scope, reason) for one probed sheet. shape is '' when nothing can load it."""
    for shape, why in (("track", _is_track_sheet(rows)), ("reactive", _is_reactive_sheet(rows))):
        if not why:
            continue
        excluded = _out_of_scope_reason(sheet_name, rows)
        if excluded:
            return shape, OUT_OF_SCOPE, excluded
        return shape, REQUIRED, why
    return "", NOT_LOADABLE, ""


# ---------------------------------------------------------------------------------------------
# the census
# ---------------------------------------------------------------------------------------------
def census_tcr_sheets(paper_dir: str | pathlib.Path,
                      with_structure: bool = False) -> SourceCensus:
    """PURE FUNCTION: paper directory -> every supplementary sheet that IS a loadable TCR shape.

    Never raises for an unreadable/foreign workbook -- an un-probeable file is recorded in
    `files_failed` and contributes NOTHING to `required`, because the gate above this must block
    only on a sheet the census positively identified, never on one it failed to read."""
    root = pathlib.Path(paper_dir)
    required: list[CensusSheet] = []
    out_of_scope: list[CensusSheet] = []
    not_loadable: list[tuple[str, str]] = []
    reactivity: list[CensusSheet] = []
    failed: list[str] = []
    files = sorted(p for p in root.rglob("*")
                   if p.suffix.lower() in (".xlsx", ".xlsm") and not p.name.startswith("~$"))
    for path in files:
        try:
            probed = _probe_workbook(path, nrows=_PROBE_ROWS)
        except Exception:                      # noqa: BLE001 - a bad workbook must not block a record
            failed.append(str(path))
            continue
        maybe_reactivity: list[tuple[str, str]] = []
        for name, rows in probed:
            shape, scope, reason = classify_sheet(name, rows)
            if not shape:
                not_loadable.append((path.name, name))
                # Only a sheet NO TCR loader claims can be a reactivity matrix, so this never
                # competes with the classification above.
                hit = _looks_like_reactivity(rows)
                if hit:
                    maybe_reactivity.append((name, hit))
                continue
            entry = CensusSheet(file=path.name, sheet=name, shape=shape, scope=scope,
                                loader=_LOADER_TOOL[shape], reason=reason, path=str(path))
            (required if scope == REQUIRED else out_of_scope).append(entry)
        if maybe_reactivity:
            reactivity.extend(_confirm_reactivity(path, maybe_reactivity, with_structure))
    return SourceCensus(required=required, out_of_scope=out_of_scope, not_loadable=not_loadable,
                        files_scanned=len(files), files_failed=failed, reactivity=reactivity)


def _confirm_reactivity(path: pathlib.Path, candidates: list[tuple[str, str]],
                        with_structure: bool = False) -> list[CensusSheet]:
    """Phase 2: full-scan each header-prefiltered sheet and keep only the ones that really carry a
    block of scored assay columns. A sheet the scan cannot confirm contributes nothing -- like the
    TCR census, this gate must fire only on a shape it positively identified."""
    out: list[CensusSheet] = []
    try:
        with zipfile.ZipFile(path) as zf:
            members = dict(_sheet_index(zf))
            for name, hit in candidates:
                member = members.get(name)
                if not member:
                    continue
                info = _scan_scored_columns(zf, member, with_structure)
                if info.get("positives", 0) < MIN_REACTIVITY_CELLS:
                    continue
                out.append(CensusSheet(
                    file=path.name, sheet=name, shape="reactivity", scope=REQUIRED,
                    loader=_LOADER_TOOL["reactivity"], path=str(path),
                    n_expected=info["positives"], data_start=info["data_start"],
                    assay_cols=info["assay_cols"], imp_col=info["imp_col"],
                    assay_pep_col=info["assay_pep_col"], id_col=info["id_col"],
                    patient_col=info["patient_col"],
                    reason=(f"reactivity matrix: {info['n_cols']} scored assay column(s) over "
                            f"{info['n_rows']} data row(s), {info['positives']} positive call(s); "
                            f"header says {hit!r}")))
    except Exception:                          # noqa: BLE001 - never block on an unreadable workbook
        return out
    return out
