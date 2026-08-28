"""Drop TCR clonotypes whose provenance is a not-antigen-sorted sheet.

POLICY: docs/POLICY_tcr_sheet_scope_2026-08-24.md (signed 2026-08-25).
A tumour-associated / TIL / unsorted-bulk clonotype cannot answer
per-(neoantigen × patient) specificity. The coverage gate already does not
REQUIRE those sheets; this module is the record-side counterpart — a
provenance-keyed filter into a NEW path, not a JSON hand-edit.

Primary Keskin marker is section_ref 'Supplementary Table 11a'. The policy
regex on section_ref / provenance locators is the corpus-wide backup.
quoted_text is NOT searched: a reactive clonotype's quote can mention TIL
without being a TIL-sheet row.
"""
from __future__ import annotations

import copy
import re

# Same patterns as tcr_source_census.OUT_OF_SCOPE_RULES, applied to locators only.
_TABLE_11A = re.compile(r"supplementary\s+table\s+11a\b", re.I)
_TUMOUR = re.compile(r"tumou?r[\s_\-]*(associated|infiltrat\w*)|\bTILs?\b", re.I)


def _locators(cl: dict) -> str:
    parts = [str(cl.get("section_ref") or "")]
    for p in cl.get("provenance") or []:
        if isinstance(p, dict):
            parts.append(str(p.get("section_ref") or ""))
            parts.append(str(p.get("locator") or ""))
    return " ".join(parts)


def is_out_of_scope_clonotype(cl: dict) -> bool:
    hay = _locators(cl)
    if _TABLE_11A.search(hay):
        return True
    if _TUMOUR.search(hay):
        return True
    return False


def strip_out_of_scope_clonotypes(rec: dict) -> tuple[dict, list]:
    """Return (shallow-copied record with TCR lane filtered, dropped clonotypes).

    Nested chains/observations ride the parent clonotype. neoantigen_tcr_flags
    are left untouched — they are a patient×antigen axis, not a clonotype list.
    """
    out = copy.copy(rec)
    dropped, kept = [], []
    for cl in rec.get("tcr_clonotypes") or []:
        (dropped if is_out_of_scope_clonotype(cl) else kept).append(cl)
    out["tcr_clonotypes"] = kept
    return out, dropped
