"""table_map.py — a pure, closed mapping DSL: turn a spreadsheet row (dict keyed by
header) into an entity dict. No eval, no I/O, no schema dependency — deterministic and
unit-testable. Used by agent_core.table_to_entities for bulk transcription.

A column may be addressed by NAME or by POSITION. Position addressing exists for
supplementary tables with merged / blank / duplicated headers (common in Nature-style
supplements) where a needed column has no usable name — e.g. Keskin Table S5, whose
two-row header leaves the peptide-ID and affinity columns nameless. Position addressing
requires the row dict to also be keyed by 0-based integer index (agent_core._read_xlsx_dicts
and read_table_rows provide both name and index keys).

Field rule (exactly one key):
    {"col": H} | {"col_idx": N} | {"col_letter": "L"}   # a column by name / 0-based idx / Excel letter
    {"const": v} | {"template": s} | {"template_list": s}
Filter (one column ref + one operator):
    {"col"|"col_idx"|"col_letter": ..., "equals": v | "in": [...] | "not_empty": true}
Template tokens: {Header} (by name) | {#N} (0-based index) | {@L} (Excel letter)

PER-CELL FAN-OUT (`per_cell`): one sheet row -> N entities, one per MATCHING CELL.
A reactivity matrix puts several measurements side by side on one row -- Hu 33479501's
'Suppl Dataset 4b' carries seven assay columns per peptide (two timepoints x ex-vivo /
pre-stimulated / minigene / autologous-tumour), each holding 1 / 0 / 'n.d.'. Row-to-entity
mapping can only reach ONE of them, which is why that lane was hand-transcribed and swung
21-155 run-to-run. Each `per_cell` entry is a column ref + one operator + its own `fields`,
merged over (and winning against) the shared row-level `fields`:

    {"fields": {...shared...},
     "per_cell": [{"col_idx": 9,  "equals": 1,
                   "fields": {"timepoint_phase":     {"const": "post_vaccine"},
                              "assay_stimulation":   {"const": "ex_vivo"},
                              "assay_antigen_format":{"const": "peptide_pulsed"}}},
                  {"col_idx": 15, "equals": 1, "fields": {...}}]}

Without `per_cell` the behaviour is unchanged: exactly one entity per surviving row.
"""
from __future__ import annotations

import re

_TPL = re.compile(r"\{([^}]+)\}")
_COLUMN_KEYS = ("col", "col_idx", "col_letter")
_OPERATORS = ("equals", "in", "not_empty")


def col_letter_to_index(letter) -> int:
    """Excel column letter -> 0-based index ('A'->0, 'L'->11, 'AA'->26). Case-insensitive."""
    s = str(letter).strip().upper()
    if not s or not s.isalpha():
        raise ValueError(f"col_letter {letter!r} is not an Excel column letter (A, B, …, AA)")
    idx = 0
    for ch in s:
        idx = idx * 26 + (ord(ch) - ord("A") + 1)
    return idx - 1


def _col_key(spec: dict):
    """The row-dict key (str header name or int index) a col/col_idx/col_letter ref resolves to."""
    if "col" in spec:
        return spec["col"]
    if "col_idx" in spec:
        try:
            return int(spec["col_idx"])
        except (TypeError, ValueError):
            raise ValueError(f"col_idx {spec['col_idx']!r} must be a 0-based integer")
    if "col_letter" in spec:
        return col_letter_to_index(spec["col_letter"])
    return None


def _resolve_token(token: str, row: dict):
    """Resolve one {…} template token against a row: '#N' -> index, '@L' -> letter, else header name."""
    if token.startswith("#"):
        try:
            return row.get(int(token[1:]))
        except ValueError:
            return None
    if token.startswith("@"):
        try:
            return row.get(col_letter_to_index(token[1:]))
        except ValueError:
            return None
    return row.get(token)


def render_template(tpl: str, row: dict) -> str:
    """Replace each {Header}/{#N}/{@L} with its row value ('' if absent)."""
    return _TPL.sub(lambda m: str(v if (v := _resolve_token(m.group(1), row)) is not None else ""), tpl)


def _extract(rule: dict, value):
    """Optional post-processor: if the rule carries `extract` (a regex with one capture
    group), return group(1) of the first match against str(value); on no match, return the
    value unchanged. Used to normalize e.g. a sheet name 'IAP-M13' -> 'M13' for an id field."""
    pat = rule.get("extract")
    if not pat:
        return value
    m = re.search(pat, str(value))
    return m.group(1) if (m and m.groups()) else value


_RULE_KEYS = _COLUMN_KEYS + ("const", "template", "template_list")


def apply_mapping(row: dict, fields: dict) -> dict:
    """Build one entity dict from a row per the field rules. Empty column cells are
    omitted (so optional fields stay unset). Raises ValueError on an unknown or ambiguous rule.
    Any rule may add `extract` (a single-group regex) to post-process its resolved value."""
    out: dict = {}
    for fname, rule in fields.items():
        present = [k for k in _RULE_KEYS if k in rule]
        if len(present) > 1:
            # Silently picking one is how a mapping ends up quietly wrong -- e.g.
            # {"col_idx": 0, "template": "Pt{#0}"} would drop the template and emit a bare '0'.
            raise ValueError(f"field {fname!r} has {len(present)} rule keys {present}; use exactly "
                             "one of col/col_idx/col_letter/const/template/template_list "
                             "('extract' is a modifier and may accompany any of them)")
        if any(k in rule for k in _COLUMN_KEYS):
            v = row.get(_col_key(rule))
            if v is None or v == "":
                continue
            out[fname] = _extract(rule, v)
        elif "const" in rule:
            out[fname] = rule["const"]
        elif "template" in rule:
            out[fname] = _extract(rule, render_template(rule["template"], row))
        elif "template_list" in rule:
            out[fname] = [_extract(rule, render_template(rule["template_list"], row))]
        else:
            raise ValueError(f"field {fname!r} has an unknown rule {rule!r}; use one of "
                             "col/col_idx/col_letter/const/template/template_list")
    return out


def _sole_operator(spec: dict) -> tuple[str | None, str]:
    """(operator, error). Exactly one of equals / in / not_empty must be present."""
    ops = [o for o in _OPERATORS if o in spec]
    if len(ops) != 1:
        return None, ("must have exactly one operator: equals, in, or not_empty "
                      f"(got {sorted(ops) if ops else 'none'})")
    return ops[0], ""


def _matches(value, spec: dict, op: str) -> bool:
    if op == "equals":
        return value == spec["equals"]
    if op == "in":
        return value in spec["in"]
    return value not in (None, "")


def apply_filter(rows: list, flt: dict) -> tuple[bool, object]:
    """Return (True, kept_rows) or (False, error_msg). One column ref + one operator."""
    try:
        key = _col_key(flt)
    except ValueError as e:
        return False, str(e)
    op, err = _sole_operator(flt)
    if err:
        return False, "filter " + err
    return True, [r for r in rows if _matches(r.get(key), flt, op)]


def validate_per_cell(per_cell) -> str:
    """'' if `per_cell` is a usable list of cell rules, else an error message."""
    if not isinstance(per_cell, list) or not per_cell:
        return "per_cell must be a non-empty list of cell rules"
    for i, cell in enumerate(per_cell):
        if not isinstance(cell, dict):
            return f"per_cell[{i}] must be an object like {{'col_idx': 9, 'equals': 1, 'fields': {{...}}}}"
        if not any(k in cell for k in _COLUMN_KEYS):
            return f"per_cell[{i}] needs a column ref: one of col, col_idx, col_letter"
        try:
            _col_key(cell)
        except ValueError as e:
            return f"per_cell[{i}]: {e}"
        _, err = _sole_operator(cell)
        if err:
            return f"per_cell[{i}] {err}"
        fields = cell.get("fields", {})
        if not isinstance(fields, dict):
            return f"per_cell[{i}]['fields'] must be an object of field rules"
        bad = [k for k, v in fields.items() if not isinstance(v, dict)]
        if bad:
            return f"per_cell[{i}]['fields'] rule(s) {bad} must be objects like {{'const': 'ex_vivo'}}"
    return ""


def apply_mapping_cells(row: dict, fields: dict, per_cell: list) -> list[dict]:
    """Fan one row out to one entity per MATCHING cell (see the per_cell note in the module
    docstring). The shared row-level mapping is built once; each matching cell's own fields are
    layered over it, so a cell rule can override a shared default. Returns [] if no cell matches
    -- a row where every assay column reads 0 / 'n.d.' contributes nothing."""
    shared = apply_mapping(row, fields)
    out: list[dict] = []
    for cell in per_cell:
        op, err = _sole_operator(cell)
        if err:
            raise ValueError(f"per_cell rule {cell!r} {err}")
        if not _matches(row.get(_col_key(cell)), cell, op):
            continue
        item = dict(shared)
        item.update(apply_mapping(row, cell.get("fields", {})))
        out.append(item)
    return out


def _all_field_rules(fields: dict, per_cell: list | None):
    """Every field rule in play — the shared mapping plus each per_cell entry's own fields."""
    yield from fields.values()
    for cell in (per_cell or []):
        if isinstance(cell, dict) and isinstance(cell.get("fields"), dict):
            yield from cell["fields"].values()


def _index_refs(fields: dict, flt: dict | None, per_cell: list | None = None) -> set:
    """0-based column indices referenced positionally (col_idx/col_letter + {#N}/{@L} tokens)."""
    refs: set = set()

    def add(spec):
        if "col_idx" in spec or "col_letter" in spec:
            try:
                refs.add(_col_key(spec))
            except ValueError:
                pass

    for rule in _all_field_rules(fields, per_cell):
        add(rule)
        for key in ("template", "template_list"):
            if key in rule:
                for tok in _TPL.findall(rule[key]):
                    try:
                        if tok.startswith("#"):
                            refs.add(int(tok[1:]))
                        elif tok.startswith("@"):
                            refs.add(col_letter_to_index(tok[1:]))
                    except ValueError:
                        pass
    if flt:
        add(flt)
    for cell in (per_cell or []):        # the measured cell itself
        if isinstance(cell, dict):
            add(cell)
    return refs


def missing_columns(fields: dict, flt: dict | None, headers: list,
                    per_cell: list | None = None) -> list:
    """Sorted list of column refs (names or out-of-range positions) the mapping/filter
    references but the sheet lacks. Names are checked against `headers`; positional refs
    (col_idx/col_letter, {#N}/{@L}) are checked against the column count (len(headers))."""
    needed: set = set()
    for rule in _all_field_rules(fields, per_cell):
        if "col" in rule:
            needed.add(rule["col"])
        for key in ("template", "template_list"):
            if key in rule:
                needed.update(t for t in _TPL.findall(rule[key])
                              if not t.startswith(("#", "@")))
    if flt and "col" in flt:
        needed.add(flt["col"])
    for cell in (per_cell or []):
        if isinstance(cell, dict) and "col" in cell:
            needed.add(cell["col"])
    missing = [c for c in needed if c not in headers]
    width = len(headers)
    missing += [f"idx {i} (width {width})" for i in _index_refs(fields, flt, per_cell)
                if not (0 <= i < width)]
    return sorted(missing)
