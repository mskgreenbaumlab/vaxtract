#!/usr/bin/env python3
"""
score_extraction.py — grade an agent's ExtractedPaper JSON against a GOLD key.

Computes precision / recall / F1 per category (neoantigens, immunogenic
neoantigens, epitopes, magnitudes, survival) plus scalar checks (cohort_size,
per-patient immunogenic counts), and LISTS EVERY false positive and false
negative explicitly — nothing is summarised away.

Honesty rules baked in:
  - Refuses to grade against an unverified gold key (verified=False) unless
    --allow-unverified is passed (then it screams the caveat in the output).
  - A "false positive" against a possibly-incomplete key is labelled
    "agent-only (agent error OR gold gap)" — not silently scored as wrong.
  - All matching is set-based on normalised identities, so it's order-independent.

Neoantigen-identity adapter (2026-06-26): a neoantigen's identity is
(patient, gene, mutation). The agent's ImmunizingPeptide carries gene+mutation but
NOT the patient; the NeoantigenCandidate that selected the peptide does
(candidate.selected_peptide_id -> candidate.patient_paper_id). `build_peptide_patient_map`
recovers the patient via that link so agent neoantigens align with the gold grain.
Without it, the patient slot is empty and every match fails (tp=0) for a representation
reason, not a real miss.

Usage:
    python3 score_extraction.py GOLD_key.json AGENT_extracted.json [--report out.json] [--allow-unverified]
Exit 0 iff the key is usable AND every category meets its threshold.
"""
from __future__ import annotations
import sys, json, re, pathlib

# ----- tunable tolerances / thresholds (explicit, not hidden) -----------------
MAG_REL_TOL = 0.05      # numeric magnitude: within 5%
SURV_MEDIAN_ABS_TOL = 0.5   # months
HR_ABS_TOL = 0.02
THRESHOLDS = {          # min F1 (or accuracy) to pass each category
    "neoantigens_all": 0.90, "neoantigens_immunogenic": 0.90,
    "epitopes": 0.90, "magnitude_presence": 0.85, "survival": 0.90,
}

# ----- normalisation ----------------------------------------------------------
def nseq(s): return re.sub(r"\s+", "", (s or "").upper())
def nmut(s):
    if not s: return ""
    return re.sub(r"^P\.", "", (s or "").strip().upper())
def narm(s): return re.sub(r"\s+", " ", (s or "").strip().lower())
def nallele(s):
    if not s: return ""
    return (s or "").strip().upper().replace(" ", "")
# Longest alternative FIRST: regex alternation is first-match-wins, so "P" must come
# last or it shadows "PT"/"PATIENT" (e.g. "Pt1" -> "t1", "Patient 3" -> "atient 3").
def npt(s): return re.sub(r"^(patient|subject|case|pt|p)[\s._-]*", "", (s or "").strip(), flags=re.I)

# ----- neoantigen-identity adapter (importable + unit-tested) -----------------
def neo_identity(patient, gene, mutation):
    """Canonical (patient, gene, mutation) key, normalised for set matching."""
    return (npt(patient or ""), (gene or "").upper(), nmut(mutation))

def build_peptide_patient_map(agent):
    """peptide paper_local_id -> patient_paper_id.

    ImmunizingPeptide carries no patient field; recover it from the candidate that
    selected the peptide (candidate.selected_peptide_id -> candidate.patient_paper_id)."""
    m = {}
    for c in agent.get("candidates", []) or []:
        pid = c.get("selected_peptide_id"); pat = c.get("patient_paper_id")
        if pid and pat and pid not in m:
            m[pid] = pat
    return m

def imp_identity(i, pep_patient):
    """Neoantigen identity of an immunizing peptide, with patient recovered from
    the candidate link when the peptide itself does not carry one."""
    patient = i.get("patient_paper_id") or pep_patient.get(i.get("paper_local_id"), "")
    return neo_identity(patient, i.get("gene_symbol"), i.get("mutation"))

def agent_neoantigen_set(agent, pep_patient=None):
    """All (patient, gene, mutation) neoantigen identities the agent represents,
    unioned across candidates and immunizing peptides (patient recovered for peptides).
    Degenerate identities with no gene are dropped."""
    pep_patient = pep_patient if pep_patient is not None else build_peptide_patient_map(agent)
    out = set()
    for c in agent.get("candidates", []) or []:
        out.add(neo_identity(c.get("patient_paper_id"), c.get("gene_symbol"), c.get("mutation")))
    for i in agent.get("immunizing_peptides", []) or []:
        out.add(imp_identity(i, pep_patient))
    return {t for t in out if t[1]}   # require a gene symbol


def main(argv):
    # ----- args ---------------------------------------------------------------
    GOLD_PATH, AGENT_PATH = argv[1], argv[2]
    report_path = None; allow_unverified = False
    for i, a in enumerate(argv[3:]):
        if a == "--report": report_path = argv[4 + i]
        if a == "--allow-unverified": allow_unverified = True
    strict = "--strict" in argv

    gold = json.loads(pathlib.Path(GOLD_PATH).read_text())
    agent = json.loads(pathlib.Path(AGENT_PATH).read_text())

    print("=" * 78)
    print(f"GOLD  : {gold.get('paper','?')}  (pmid {gold.get('pmid','?')})  verified={gold.get('verified')}")
    print(f"AGENT : {AGENT_PATH}")
    print("=" * 78)
    if not gold.get("verified"):
        msg = ("GOLD KEY IS UNVERIFIED — it has NOT been human-checked against the source. "
               "Scores below are NOT a valid measure of agent quality.")
        if not allow_unverified:
            print("REFUSING TO GRADE:", msg); print("(pass --allow-unverified to score anyway, for harness testing only)")
            return 2
        print("!!! " + msg + " (scoring anyway: --allow-unverified) !!!\n")

    # ----- build agent-side index ---------------------------------------------
    imp = {i["paper_local_id"]: i for i in agent.get("immunizing_peptides", [])}
    epi = {e["paper_local_id"]: e for e in agent.get("epitopes", [])}
    pool = {p["paper_local_id"]: p for p in agent.get("pools", [])}
    pep_patient = build_peptide_patient_map(agent)

    agent_all_neo = agent_neoantigen_set(agent, pep_patient)

    # immunogenic neoantigens per agent: resolve immunogenic/positive evidence -> IMP identities
    def evidence_to_imps(e):
        out = []
        if e["target_kind"] == "immunizing_peptide":
            i = imp.get(e.get("immunizing_peptide_paper_id"));  out += [i] if i else []
        elif e["target_kind"] == "epitope":
            ep = epi.get(e.get("epitope_paper_id")) or {}
            for pid in ep.get("parent_peptide_ids", []):
                if pid in imp: out.append(imp[pid])
        elif e["target_kind"] == "pool":
            p = pool.get(e.get("pool_paper_id")) or {}
            for pid in p.get("member_peptide_ids", []):
                if pid in imp: out.append(imp[pid])
        return out

    agent_immuno_neo = set()
    agent_mag = {}   # (patient, gene, mut) -> magnitude dict
    for e in agent.get("evidence", []):
        if e.get("outcome") in ("immunogenic", "positive"):
            for i in evidence_to_imps(e):
                ident = imp_identity(i, pep_patient); agent_immuno_neo.add(ident)
                if e.get("magnitude"): agent_mag[ident] = e["magnitude"]

    agent_epi = {(nseq(e["sequence"]), e["mhc_class"]) for e in epi.values()}
    agent_surv = {(s["endpoint"], narm(s.get("arm_label", ""))): s for s in agent.get("survival_outcomes", [])}

    # ----- gold-side sets -----------------------------------------------------
    gold_all_neo = {neo_identity(n["patient"], n["gene"], n.get("mutation")) for n in gold.get("neoantigens", [])}
    gold_immuno_neo = {neo_identity(n["patient"], n["gene"], n.get("mutation")) for n in gold.get("neoantigens", []) if n.get("immunogenic")}
    gold_epi = {(nseq(e["sequence"]), e["mhc_class"]) for e in gold.get("epitopes", [])}
    gold_mag = {neo_identity(m["patient"], m.get("gene"), m.get("mutation")): m for m in gold.get("magnitudes", [])}
    gold_surv = {(s["endpoint"], narm(s.get("arm", ""))): s for s in gold.get("survival", [])}

    # ----- set scoring helper -------------------------------------------------
    results = {}
    def score_set(name, gold_set, agent_set, render=lambda x: x):
        tp = gold_set & agent_set; fn = gold_set - agent_set; fp = agent_set - gold_set
        P = len(tp) / (len(tp) + len(fp)) if (tp or fp) else 1.0
        R = len(tp) / (len(tp) + len(fn)) if (tp or fn) else 1.0
        F = 2 * P * R / (P + R) if (P + R) else 0.0
        results[name] = {"tp": len(tp), "fp": len(fp), "fn": len(fn), "precision": round(P, 4), "recall": round(R, 4), "f1": round(F, 4)}
        print(f"\n## {name}:  P={P:.3f}  R={R:.3f}  F1={F:.3f}   (tp={len(tp)} fp={len(fp)} fn={len(fn)})")
        if fn:
            print(f"  MISSED by agent (false negatives) — {len(fn)}:")
            for x in sorted(fn): print(f"    - {render(x)}")
        if fp:
            print(f"  AGENT-ONLY (false positives = agent error OR gold gap) — {len(fp)}:")
            for x in sorted(fp): print(f"    + {render(x)}")
        return F

    neo_r = lambda x: f"{x[0]} {x[1]} {x[2]}"
    score_set("neoantigens_all", gold_all_neo, agent_all_neo, neo_r)
    score_set("neoantigens_immunogenic", gold_immuno_neo, agent_immuno_neo, neo_r)
    score_set("epitopes", gold_epi, agent_epi, lambda x: f"{x[0]} (MHC-{x[1]})")

    # ----- magnitude accuracy (only over targets the gold has a magnitude for) -
    print("\n## magnitudes (accuracy over gold-keyed targets):")
    mag_keys = set(gold_mag); present = wrong = correct = 0; mag_detail = []
    for k in sorted(mag_keys):
        g = gold_mag[k]; a = agent_mag.get(k)
        if a is None:
            mag_detail.append(f"    - MISSING magnitude for {neo_r(k)} (gold={g.get('value') or g.get('grade')})"); continue
        present += 1
        ok = True
        if g.get("grade") is not None:
            ok = (a.get("grade") == g["grade"])
        elif g.get("value") is not None and a.get("value") is not None:
            ok = abs(a["value"] - g["value"]) <= MAG_REL_TOL * abs(g["value"]) and (g.get("unit") in (None, a.get("unit")))
        elif g.get("value") is not None and a.get("value") is None:
            ok = False  # gold had a number; agent gave none
        if ok: correct += 1
        else:
            wrong += 1
            mag_detail.append(f"    ~ WRONG for {neo_r(k)}: gold={g.get('value') or g.get('grade')} {g.get('unit') or ''} | agent={a.get('value') or a.get('grade')} {a.get('unit') or ''}")
    mag_presence = present / len(mag_keys) if mag_keys else 1.0
    mag_value_acc = correct / present if present else 1.0
    results["magnitude_presence"] = {"recall": round(mag_presence, 4), "value_accuracy": round(mag_value_acc, 4),
                                     "gold_keyed": len(mag_keys), "present": present, "correct": correct, "wrong": wrong}
    print(f"  presence recall = {mag_presence:.3f} ({present}/{len(mag_keys)} gold-keyed targets had an agent magnitude)")
    print(f"  value accuracy  = {mag_value_acc:.3f} ({correct}/{present} correct within tol; grades exact)")
    for d in mag_detail: print(d)

    # ----- survival -----------------------------------------------------------
    print("\n## survival outcomes:")
    gs, as_ = set(gold_surv), set(agent_surv)
    tp = gs & as_; fn = gs - as_; fp = as_ - gs
    sP = len(tp) / (len(tp) + len(fp)) if (tp or fp) else 1.0
    sR = len(tp) / (len(tp) + len(fn)) if (tp or fn) else 1.0
    sF = 2 * sP * sR / (sP + sR) if (sP + sR) else 0.0
    surv_value_issues = []
    for k in sorted(tp):
        g, a = gold_surv[k], agent_surv[k]
        if bool(g.get("not_reached")) != bool(a.get("not_reached")):
            surv_value_issues.append(f"    ~ {k}: not_reached gold={g.get('not_reached')} agent={a.get('not_reached')}")
        if g.get("median") is not None and a.get("median_value") is not None and abs(g["median"] - a["median_value"]) > SURV_MEDIAN_ABS_TOL:
            surv_value_issues.append(f"    ~ {k}: median gold={g['median']} agent={a['median_value']}")
        if g.get("hazard_ratio") is not None and a.get("hazard_ratio") is not None and abs(g["hazard_ratio"] - a["hazard_ratio"]) > HR_ABS_TOL:
            surv_value_issues.append(f"    ~ {k}: HR gold={g['hazard_ratio']} agent={a['hazard_ratio']}")
    results["survival"] = {"precision": round(sP, 4), "recall": round(sR, 4), "f1": round(sF, 4),
                           "tp": len(tp), "fp": len(fp), "fn": len(fn), "value_issues": len(surv_value_issues)}
    print(f"  P={sP:.3f} R={sR:.3f} F1={sF:.3f}  (tp={len(tp)} fp={len(fp)} fn={len(fn)})")
    for k in sorted(fn): print(f"    - MISSED arm: {k}")
    for k in sorted(fp): print(f"    + AGENT-ONLY arm: {k}")
    for d in surv_value_issues: print(d)

    # ----- scalars ------------------------------------------------------------
    print("\n## scalar fields:")
    def scalar(name, g, a):
        flag = "" if g == a else "   <-- MISMATCH"
        print(f"  {name}: gold={g} agent={a}{flag}"); return g == a
    scalar("cohort_size", gold.get("cohort_size"), agent.get("cohort_size"))
    scalar("n_enrolled", gold.get("n_enrolled"), agent.get("n_enrolled"))
    ppi_gold = {npt(k): v for k, v in (gold.get("per_patient_immunogenic") or {}).items()}
    ppi_agent = {npt(p["paper_local_id"]): p.get("n_peptides_immunogenic") for p in agent.get("patients", [])}
    ppi_mis = {k: (ppi_gold[k], ppi_agent.get(k)) for k in ppi_gold if ppi_gold[k] != ppi_agent.get(k)}
    print(f"  per_patient_immunogenic mismatches: {len(ppi_mis)}")
    for k, (g, a) in sorted(ppi_mis.items()): print(f"    ~ patient {k}: gold={g} agent={a}")

    # ----- STRICT hard rules (always printed; gate-affecting only with --strict)
    # Threshold F1s can mask a few SERIOUS errors among hundreds of items (e.g. one
    # fabricated neoantigen keeps precision ~0.996). These rules catch what averages hide.
    strict_violations = []
    if results["neoantigens_immunogenic"]["fp"] > 0:
        strict_violations.append(f"{results['neoantigens_immunogenic']['fp']} fabricated immunogenic neoantigen(s) (agent-only)")
    if results["survival"]["fn"] > 0:
        strict_violations.append(f"{results['survival']['fn']} missing survival arm(s)")
    if results["magnitude_presence"]["wrong"] > 0:
        strict_violations.append(f"{results['magnitude_presence']['wrong']} wrong magnitude value(s)")
    _miss = results["magnitude_presence"]["gold_keyed"] - results["magnitude_presence"]["present"]
    if _miss > 0:
        strict_violations.append(f"{_miss} missing magnitude(s) where gold had one")
    if results["survival"]["value_issues"] > 0:
        strict_violations.append(f"{results['survival']['value_issues']} survival value mismatch(es)")
    print("\n## STRICT hard-rule violations:", "; ".join(strict_violations) if strict_violations else "none")
    print(f"   (strict mode {'ON — these affect the verdict' if strict else 'OFF — reported only; pass --strict to enforce'})")

    # ----- verdict ------------------------------------------------------------
    print("\n" + "=" * 78)
    fails = []
    for cat, thr in THRESHOLDS.items():
        r = results.get(cat, {})
        score = r.get("f1", r.get("recall", 0.0))
        mark = "PASS" if score >= thr else "FAIL"
        if score < thr: fails.append(cat)
        print(f"  {mark}  {cat:28} score={score:.3f}  threshold={thr}")
    verdict = (not fails) and gold.get("verified") and (not (strict and strict_violations))
    print("=" * 78)
    if verdict:
        print("  VERDICT: PASS — agent meets gold on every category" + (" (strict)" if strict else "") + ".")
    elif not gold.get("verified"):
        print("  VERDICT: INDETERMINATE — gold key not verified; fix the key, then re-run.")
    elif strict and strict_violations:
        print(f"  VERDICT: FAIL (strict) — hard-rule violations: {'; '.join(strict_violations)}"
              + (f"; below threshold on: {', '.join(fails)}" if fails else ""))
    else:
        print(f"  VERDICT: FAIL — below threshold on: {', '.join(fails)}")
    print("=" * 78)

    if report_path:
        pathlib.Path(report_path).write_text(json.dumps(
            {"gold": gold.get("paper"), "agent": AGENT_PATH, "verified": gold.get("verified"),
             "results": results, "per_patient_mismatches": ppi_mis, "failed_categories": fails}, indent=2, default=str))
        print(f"  metrics written to {report_path}")

    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
