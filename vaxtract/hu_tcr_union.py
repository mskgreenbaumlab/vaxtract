"""Union extra Hu TCR clonotypes onto a loader-built evidence record.

Base = rxload (17/279/11 evidence, 269 tet-spec tracks).
Donors = promptfix Dataset 9 clones (by chain identity not already in base)
         + old keep-file Dataset 10 clones (same).
Flags = promptfix neoantigen_tcr_flags whose IMP id exists on the base.
Dataset 10 epitopes are copied when the sequence is missing, parent IMP must
already exist. Provenance text is never rewritten.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path


def chain_key(cl: dict) -> tuple:
    chains = tuple(sorted(
        ((ch.get("junction_aa") or ch.get("cdr3_aa") or ""), ch.get("locus") or "")
        for ch in (cl.get("chains") or [])
    ))
    return (cl.get("patient_paper_id"), chains or (cl.get("paper_local_id"),))


def union_clonotypes(base: dict, *donors: dict) -> tuple[dict, dict]:
    """Return (new record, report). Does not finalize."""
    out = copy.deepcopy(base)
    seen = {chain_key(c) for c in out.get("tcr_clonotypes") or []}
    taken_ids = {c.get("paper_local_id") for c in out.get("tcr_clonotypes") or []}
    epi_seq = {(e.get("sequence") or "").upper(): e["paper_local_id"]
               for e in out.get("epitopes") or [] if e.get("sequence") and e.get("paper_local_id")}
    epi_ids = {e.get("paper_local_id") for e in out.get("epitopes") or []}
    imp_ids = {p.get("paper_local_id") for p in out.get("immunizing_peptides") or []}
    added, skipped, epitopes_added = [], [], []

    def adopt_epitope(cl: dict, donor: dict) -> str | None:
        eid = cl.get("epitope_paper_id")
        if not eid:
            return None
        if eid in epi_ids:
            return eid
        src = next((e for e in (donor.get("epitopes") or []) if e.get("paper_local_id") == eid), None)
        if not src:
            return None
        seq = (src.get("sequence") or "").upper()
        if seq and seq in epi_seq:
            cl["epitope_paper_id"] = epi_seq[seq]
            return epi_seq[seq]
        parents = [p for p in (src.get("parent_peptide_ids") or []) if p in imp_ids]
        if src.get("parent_peptide_ids") and not parents:
            return None
        new_e = copy.deepcopy(src)
        new_e["parent_peptide_ids"] = parents
        out.setdefault("epitopes", []).append(new_e)
        epi_ids.add(eid)
        if seq:
            epi_seq[seq] = eid
        epitopes_added.append(eid)
        return eid

    for donor in donors:
        for cl in donor.get("tcr_clonotypes") or []:
            k = chain_key(cl)
            if k in seen:
                skipped.append(cl.get("paper_local_id"))
                continue
            row = copy.deepcopy(cl)
            if row.get("specific_for_kind") == "epitope":
                if adopt_epitope(row, donor) is None:
                    skipped.append(row.get("paper_local_id"))
                    continue
            elif row.get("specific_for_kind") == "immunizing_peptide":
                if row.get("immunizing_peptide_paper_id") not in imp_ids:
                    skipped.append(row.get("paper_local_id"))
                    continue
            cid = row.get("paper_local_id")
            if cid in taken_ids:
                n = 2
                while f"{cid}_{n}" in taken_ids:
                    n += 1
                row["paper_local_id"] = f"{cid}_{n}"
            taken_ids.add(row["paper_local_id"])
            seen.add(k)
            out.setdefault("tcr_clonotypes", []).append(row)
            added.append(row["paper_local_id"])

    base_flags = {(f.get("patient_paper_id"), f.get("target_kind"),
                   f.get("immunizing_peptide_paper_id") or f.get("epitope_paper_id")
                   or f.get("candidate_paper_id"))
                  for f in out.get("neoantigen_tcr_flags") or []}
    flags_added = 0
    for donor in donors:
        for f in donor.get("neoantigen_tcr_flags") or []:
            tid = f.get("immunizing_peptide_paper_id") or f.get("epitope_paper_id") or f.get("candidate_paper_id")
            if f.get("target_kind") == "immunizing_peptide" and tid not in imp_ids:
                continue
            if f.get("target_kind") == "epitope" and tid not in epi_ids:
                continue
            key = (f.get("patient_paper_id"), f.get("target_kind"), tid)
            if key in base_flags:
                continue
            out.setdefault("neoantigen_tcr_flags", []).append(copy.deepcopy(f))
            base_flags.add(key)
            flags_added += 1

    if any("Dataset 10" in (c.get("section_ref") or "") or "Dataset 9" in (c.get("section_ref") or "")
           for c in out.get("tcr_clonotypes") or []):
        # tracks are bulk; cloned/multimer clones are single-cell/functional
        out["tcr_seq_status"] = "both"

    report = {"clonotypes_added": len(added), "clonotypes_skipped": len(skipped),
              "epitopes_added": epitopes_added, "flags_added": flags_added,
              "n_clonotypes": len(out.get("tcr_clonotypes") or [])}
    return out, report


def write_partial(rec: dict, dest: str, *, sidecar_from: str | None = None) -> None:
    from . import agent_core as ac
    import shutil
    p = Path(dest)
    p.parent.mkdir(parents=True, exist_ok=True)
    rec = copy.deepcopy(rec)
    rec.pop("finalize_overrides_used", None)
    Path(ac._partial_path(dest)).write_text(json.dumps(rec))
    Path(ac._lock_path(dest)).unlink(missing_ok=True)
    Path(ac._quarantine_path(dest)).unlink(missing_ok=True)
    if sidecar_from:
        for fn in (ac._ledger_path, ac._source_path):
            src, dst = Path(fn(sidecar_from)), Path(fn(dest))
            if src.is_file():
                shutil.copy2(src, dst)
            else:
                dst.unlink(missing_ok=True)
    else:
        Path(ac._ledger_path(dest)).unlink(missing_ok=True)
        Path(ac._source_path(dest)).unlink(missing_ok=True)
