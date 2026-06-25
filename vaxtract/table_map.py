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
"""
from __future__ import annotations

import re

_TPL = re.compile(r"\{([^}]+)\}")
_COLUMN_KEYS = ("col", "col_idx", "col_letter")


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


def apply_mapping(row: dict, fields: dict) -> dict:
    """Build one entity dict from a row per the field rules. Empty column cells are
    omitted (so optional fields stay unset). Raises ValueError on an unknown rule.
    Any rule may add `extract` (a single-group regex) to post-process its resolved value."""
    out: dict = {}
    for fname, rule in fields.items():
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


def apply_filter(rows: list, flt: dict) -> tuple[bool, object]:
    """Return (True, kept_rows) or (False, error_msg). One column ref + one operator."""
    try:
        key = _col_key(flt)
    except ValueError as e:
        return False, str(e)
    if "equals" in flt:
        return True, [r for r in rows if r.get(key) == flt["equals"]]
    if "in" in flt:
        vals = flt["in"]
        return True, [r for r in rows if r.get(key) in vals]
    if "not_empty" in flt:
        return True, [r for r in rows if r.get(key) not in (None, "")]
    return False, "filter must have exactly one operator: equals, in, or not_empty"


def _index_refs(fields: dict, flt: dict | None) -> set:
    """0-based column indices referenced positionally (col_idx/col_letter + {#N}/{@L} tokens)."""
    refs: set = set()

    def add(spec):
        if "col_idx" in spec or "col_letter" in spec:
            try:
                refs.add(_col_key(spec))
            except ValueError:
                pass

    for rule in fields.values():
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
    return refs


def missing_columns(fields: dict, flt: dict | None, headers: list) -> list:
    """Sorted list of column refs (names or out-of-range positions) the mapping/filter
    references but the sheet lacks. Names are checked against `headers`; positional refs
    (col_idx/col_letter, {#N}/{@L}) are checked against the column count (len(headers))."""
    needed: set = set()
    for rule in fields.values():
        if "col" in rule:
            needed.add(rule["col"])
        for key in ("template", "template_list"):
            if key in rule:
                needed.update(t for t in _TPL.findall(rule[key])
                              if not t.startswith(("#", "@")))
    if flt and "col" in flt:
        needed.add(flt["col"])
    missing = [c for c in needed if c not in headers]
    width = len(headers)
    missing += [f"idx {i} (width {width})" for i in _index_refs(fields, flt)
                if not (0 <= i < width)]
    return sorted(missing)
