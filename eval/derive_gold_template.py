#!/usr/bin/env python3
"""
derive_gold_template.py — emit a STARTER gold key from a validated ExtractedPaper
JSON, so a human corrects/confirms rather than authoring from scratch.

CRITICAL: the output has verified=False on purpose. It is machine-derived from a
PRIOR extraction, so it is NOT ground truth — using it to grade an agent without
human verification measures agreement with the prior model, not correctness.
A human MUST review every field against the source and set verified=true.

Usage:  python3 derive_gold_template.py EXTRACTED.json [GOLD_template.json]
"""
from __future__ import annotations
import sys, json, pathlib

src = json.loads(pathlib.Path(sys.argv[1]).read_text())
out = sys.argv[2] if len(sys.argv) > 2 else sys.argv[1].replace(".json", "_gold_TEMPLATE.json")

imp = {i["paper_local_id"]: i for i in src.get("immunizing_peptides", [])}
epi = {e["paper_local_id"]: e for e in src.get("epitopes", [])}
pool = {p["paper_local_id"]: p for p in src.get("pools", [])}

# which IMPs are INDIVIDUALLY immunogenic (per evidence).
# IMPORTANT: a pool-level positive response does NOT make each pool member
# individually immunogenic — the pool must be deconvoluted to the specific
# peptide(s). So only peptide-target and epitope-target evidence count here;
# pool-target evidence is pool-level (captured as qualitative magnitude), not
# per-peptide. (Counting pool members here over-marks immunogenicity.)
immuno = set(); mag_by_imp = {}
for e in src.get("evidence", []):
    if e.get("outcome") in ("immunogenic", "positive"):
        ids = []
        if e["target_kind"] == "immunizing_peptide" and e.get("immunizing_peptide_paper_id") in imp:
            ids = [e["immunizing_peptide_paper_id"]]
        elif e["target_kind"] == "epitope":
            ids = (epi.get(e.get("epitope_paper_id")) or {}).get("parent_peptide_ids", [])
        # NOTE: intentionally NOT expanding pool members into `immuno`.
        for pid in ids:
            immuno.add(pid)
            if e.get("magnitude"): mag_by_imp[pid] = e["magnitude"]

neoantigens = [{"patient": i.get("patient_paper_id"), "gene": i.get("gene_symbol"),
                "mutation": i.get("mutation"), "immunogenic": i["paper_local_id"] in immuno}
               for i in imp.values()]
epitopes = [{"sequence": e["sequence"], "mhc_class": e["mhc_class"], "allele": e.get("hla_allele")}
            for e in epi.values()]
magnitudes = []
for pid, m in mag_by_imp.items():
    i = imp.get(pid, {})
    magnitudes.append({"patient": i.get("patient_paper_id"), "gene": i.get("gene_symbol"),
                       "mutation": i.get("mutation"), "value": m.get("value"),
                       "unit": m.get("unit") if m.get("value") is not None else None,
                       "grade": m.get("grade")})
survival = [{"endpoint": s["endpoint"], "arm": s.get("arm_label", ""), "median": s.get("median_value"),
             "not_reached": s.get("not_reached", False), "hazard_ratio": s.get("hazard_ratio")}
            for s in src.get("survival_outcomes", [])]
ppi = {p["paper_local_id"]: p.get("n_peptides_immunogenic") for p in src.get("patients", [])}

key = {
    "pmid": src.get("pmid", ""), "paper": src.get("title", "")[:60],
    "source_refs": "", "verified": False, "verified_by": None,
    "cohort_size": src.get("cohort_size"), "n_enrolled": src.get("n_enrolled"),
    "per_patient_immunogenic": ppi,
    "neoantigens": neoantigens, "epitopes": epitopes, "magnitudes": magnitudes, "survival": survival,
    "notes": ("AUTO-DERIVED FROM A PRIOR EXTRACTION — NOT GROUND TRUTH. A human must verify every "
              "field against the source paper/supplement and set verified=true before using this to "
              "grade any agent."),
}
pathlib.Path(out).write_text(json.dumps(key, indent=2, default=str))
print(f"wrote {out}  (verified=FALSE)")
print(f"  neoantigens={len(neoantigens)} (immunogenic={sum(n['immunogenic'] for n in neoantigens)}) "
      f"epitopes={len(epitopes)} magnitudes={len(magnitudes)} survival={len(survival)}")
print("  >>> NEXT: a human verifies against the source, fixes errors, sets verified=true. <<<")
