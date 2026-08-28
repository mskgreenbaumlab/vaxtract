"""H3 row gate for TCR clonotypes — schema-forced needs_review only.

DESIGN H3 / D9 (docs/website/tcr-on-site/DESIGN_tcr_on_atlas.md): public /tcrs
hides COALESCE(needs_review,false)=true so inferred_association cannot be
laundered into a hard TCR→antigen edge.

That is not a dump for "the adapter was conservative." Antigen-sorted table
rows (Keskin 8a/8b/11b/11c, Hu Dataset 8 tracks, Hu Dataset 9 cloned TCRs)
load public unless the schema itself requires the flag:

- specificity_evidence == inferred_association
- >1 observation AND observations_identity == not_asserted

neoantigen_tcr_flags have no schema-forced needs_review. A derived or
agent-stamped flag from an antigen-sorted sheet is public. Evidence and
epitopes are untouched (different gates).
"""
from __future__ import annotations

import copy


def clonotype_h3_needs_review(cl: dict) -> bool:
    if (cl.get("specificity_evidence") or "") == "inferred_association":
        return True
    obs = cl.get("observations") or []
    identity = cl.get("observations_identity") or "not_asserted"
    if len(obs) > 1 and identity == "not_asserted":
        return True
    return False


def flag_h3_needs_review(_flag: dict) -> bool:
    """Flags have no inferred_association field. Adapter/agent stamps are not H3."""
    return False


def _clear_lane(rows: list, must_fn, extra_keys: tuple[str, ...]) -> tuple[list, list[dict]]:
    cleared, new = [], []
    for row in rows:
        must = must_fn(row)
        was = bool(row.get("needs_review"))
        out = copy.deepcopy(row)
        out["needs_review"] = must
        if was and not must:
            info = {k: row.get(k) for k in extra_keys}
            info["section_ref"] = row.get("section_ref")
            cleared.append(info)
        new.append(out)
    return new, cleared


def apply_h3_to_clonotypes(rec: dict) -> tuple[dict, list[dict]]:
    """Return (deep-copied record, list of clonotypes whose flag was cleared)."""
    out = copy.deepcopy(rec)
    new, cleared = _clear_lane(
        out.get("tcr_clonotypes") or [],
        clonotype_h3_needs_review,
        ("paper_local_id", "specificity_evidence"),
    )
    out["tcr_clonotypes"] = new
    return out, cleared


def apply_h3_to_flags(rec: dict) -> tuple[dict, list[dict]]:
    """Return (deep-copied record, list of flags whose needs_review was cleared)."""
    out = copy.deepcopy(rec)
    new, cleared = _clear_lane(
        out.get("neoantigen_tcr_flags") or [],
        flag_h3_needs_review,
        ("patient_paper_id", "tcr_identified", "method_label"),
    )
    out["neoantigen_tcr_flags"] = new
    return out, cleared


def apply_h3(rec: dict) -> tuple[dict, list[dict], list[dict]]:
    """Clonotypes then flags. Returns (record, cleared_clonotypes, cleared_flags)."""
    out, cl_cleared = apply_h3_to_clonotypes(rec)
    out, fl_cleared = apply_h3_to_flags(out)
    return out, cl_cleared, fl_cleared
