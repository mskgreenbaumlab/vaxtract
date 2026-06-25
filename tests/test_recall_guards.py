"""Prototype recall guards (2026-06-08): two block-once, overridable finalize nudges that catch the
silver-tier under-extractions diagnosed on the scale lane.

1. _peptide_recall_gap — anchors peptide count to the paper-stated n_selected_reported, closing the
   GAMEABLE strict reconciliation (39972124 lowered n_peptides_synthesized 108->3 to match 3 scraped
   peptides and passed).
2. _evidence_breadth_gap — flags when responses cover only a few of the vaccinated patients
   (33064988: evidence for 6 of 48 patients because the per-patient sheets weren't all iterated).
"""
import json
import pathlib
from types import SimpleNamespace as NS

import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]
META = {"pmid": "12345678", "journal": "Test J", "year": 2024, "title": "Recall guard test",
        "cohort_size": 4, "indication_summary": "melanoma"}


# ---- gap detectors (duck-typed) ----

def test_peptide_recall_no_anchor_no_nudge():
    rec = NS(n_selected_reported=None, immunizing_peptides=[1, 2, 3])
    assert agent_core._peptide_recall_gap(rec) is None

def test_peptide_recall_gap_flagged_when_under_anchor():
    rec = NS(n_selected_reported=108, immunizing_peptides=[1, 2, 3])  # the 39972124 shape
    msg = agent_core._peptide_recall_gap(rec)
    assert msg and "108" in msg and "3" in msg and "n_peptides_synthesized" in msg

def test_peptide_recall_clean_when_met_or_exceeded():
    assert agent_core._peptide_recall_gap(NS(n_selected_reported=2, immunizing_peptides=[1, 2])) is None
    assert agent_core._peptide_recall_gap(NS(n_selected_reported=2, immunizing_peptides=[1, 2, 3])) is None


# ---- companion-paper exemption (v2.11.4): a secondary-analysis paper that defers its manifest ----

def test_peptide_recall_exempt_when_companion_paper_set():
    # the 39972124 shape: 3 named neoepitopes vs a cited cohort count of 108, manifest in Rojas.
    rec = NS(n_selected_reported=108, immunizing_peptides=[1, 2, 3],
             companion_paper_ref="Rojas et al. Nature 2023; PMID 37165196")
    assert agent_core._peptide_recall_gap(rec) is None  # deferral declared -> not a miss

def test_peptide_recall_still_fires_without_companion_ref():
    # same shortfall but NO companion ref -> still a real miss (no escape hatch)
    rec = NS(n_selected_reported=108, immunizing_peptides=[1, 2, 3], companion_paper_ref=None)
    assert agent_core._peptide_recall_gap(rec) is not None


def _epos(pid, tgt):
    return NS(patient_paper_id=pid, outcome="immunogenic", immunizing_peptide_paper_id=tgt,
              epitope_paper_id=None, pool_paper_id=None, candidate_paper_id=None)

def _eneg(pid, tgt):
    return NS(patient_paper_id=pid, outcome="not_immunogenic", immunizing_peptide_paper_id=tgt,
              epitope_paper_id=None, pool_paper_id=None, candidate_paper_id=None)

def test_evidence_anchor_immunogenic_side_exempt_when_companion_set():
    # recorded 2 immunogenic vs cited 25 -> would normally nudge; companion ref relaxes it.
    rec = NS(evidence=[_epos("P1", "i1"), _epos("P2", "i2")],
             n_immunogenic_reported=25, n_tested_negative_reported=None,
             companion_paper_ref="PMID 37165196")
    assert agent_core._evidence_anchor_gap(rec) == []

def test_evidence_anchor_immunogenic_side_fires_without_companion():
    rec = NS(evidence=[_epos("P1", "i1"), _epos("P2", "i2")],
             n_immunogenic_reported=25, n_tested_negative_reported=None, companion_paper_ref=None)
    gaps = agent_core._evidence_anchor_gap(rec)
    assert gaps and "MISSED" in gaps[0]

def test_evidence_anchor_negative_exempt_under_companion():
    # v2.11.5: a secondary paper defers its NEGATIVE characterization to the companion too (39972124:
    # n_tested_negative_reported=83 vs 0 recorded). So even an over-the-tolerance negative count is
    # exempt under companion_paper_ref (was the residual evidence_count_mismatch on 39972124).
    rec = NS(evidence=[_eneg("P1", "i1"), _eneg("P2", "i2"), _eneg("P3", "i3")],
             n_immunogenic_reported=25, n_tested_negative_reported=1,
             companion_paper_ref="PMID 37165196")
    assert agent_core._evidence_anchor_gap(rec) == []  # both sides exempt under companion

def test_evidence_anchor_negative_over_fires_without_companion():
    # the same shape WITHOUT companion: 3 not_immunogenic vs 1 reported (3 > 1.5x) -> over-enumeration fires.
    rec = NS(evidence=[_eneg("P1", "i1"), _eneg("P2", "i2"), _eneg("P3", "i3")],
             n_immunogenic_reported=None, n_tested_negative_reported=1, companion_paper_ref=None)
    gaps = agent_core._evidence_anchor_gap(rec)
    assert gaps and any("OVER-enumerated" in g for g in gaps)


def _pat(pid, nsyn):
    return NS(paper_local_id=pid, n_peptides_synthesized=nsyn)

def _ev(pid):
    return NS(patient_paper_id=pid)

def test_breadth_small_cohort_never_fires():
    pats = [_pat(f"P{i}", 1) for i in range(5)]  # <6 vaccinated -> skip
    rec = NS(patients=pats, evidence=[_ev("P0")])
    assert agent_core._evidence_breadth_gap(rec) is None

def test_breadth_gap_flagged_when_coverage_thin():
    # the 33064988 shape (scaled): 12 vaccinated, evidence for 2 (<1/3) -> fires
    pats = [_pat(f"P{i}", 1) for i in range(12)]
    rec = NS(patients=pats, evidence=[_ev("P0"), _ev("P1")])
    msg = agent_core._evidence_breadth_gap(rec)
    assert msg and "2 of 12" in msg

def test_breadth_clean_at_plausible_response_rate():
    # Rojas/Keskin-like 50-62% coverage is real biology, NOT under-extraction -> no nudge
    pats = [_pat(f"P{i}", 1) for i in range(12)]
    rec = NS(patients=pats, evidence=[_ev(f"P{i}") for i in range(6)])  # 6/12 = 50%
    assert agent_core._evidence_breadth_gap(rec) is None

def test_breadth_only_counts_vaccinated_patients():
    # patients with 0 synthesized peptides are not "vaccinated" -> don't dilute the denominator
    pats = [_pat(f"P{i}", 2) for i in range(6)] + [_pat(f"S{i}", 0) for i in range(20)]
    rec = NS(patients=pats, evidence=[_ev(f"P{i}") for i in range(4)])  # 4/6 vaccinated >= 1/3
    assert agent_core._evidence_breadth_gap(rec) is None


# ---- finalize blocks once, then override proceeds ----

def _patient(pid, nsyn):
    return {"quoted_text": f"pt {pid}", "section_ref": "M", "paper_local_id": pid,
            "indication": "melanoma", "n_peptides_synthesized": nsyn, "n_peptides_immunogenic": 0}

def _imp(pid):
    return {"paper_local_id": pid, "sequence": "SLLQHLIGL", "is_neoantigen": True,
            "quoted_text": "q", "section_ref": "s"}

def test_finalize_peptide_recall_blocks_then_overrides(tmp_path):
    out = tmp_path / "r.json"
    meta = {**META, "cohort_size": 1, "n_selected_reported": 5}  # paper says 5 selected
    agent_core.init_partial(str(out), json.dumps(meta))
    agent_core.append_section(str(out), "patients", json.dumps([_patient("P1", 2)]))
    agent_core.append_section(str(out), "immunizing_peptides", json.dumps([_imp("i1"), _imp("i2")]))
    ok, msg = agent_core.finalize_partial(str(out), allow_missing_magnitudes=True)
    assert not ok and "peptide recall" in msg and "allow_peptide_count_mismatch=true" in msg
    ok2, msg2 = agent_core.finalize_partial(str(out), allow_missing_magnitudes=True,
                                            allow_peptide_count_mismatch=True)
    assert ok2 and "OVERRIDES USED" in msg2
    rec = json.loads(out.read_text())   # v2.11.3: the override is persisted for needs_review routing
    assert rec["finalize_overrides_used"] == ["allow_peptide_count_mismatch"]

def test_finalize_breadth_blocks_then_overrides(tmp_path):
    out = tmp_path / "r.json"
    agent_core.init_partial(str(out), json.dumps({**META, "cohort_size": 6}))
    agent_core.append_section(str(out), "patients",
                              json.dumps([_patient(f"P{i}", 1) for i in range(6)]))
    agent_core.append_section(str(out), "immunizing_peptides",
                              json.dumps([_imp(f"i{i}") for i in range(6)]))
    # one immunogenic response, covering only P0 of the 4 vaccinated
    agent_core.append_section(str(out), "evidence", json.dumps([{
        "patient_paper_id": "P0", "target_kind": "immunizing_peptide",
        "immunizing_peptide_paper_id": "i0", "assay": "elispot", "outcome": "immunogenic",
        "quoted_text": "q", "section_ref": "s"}]))
    ok, msg = agent_core.finalize_partial(str(out), allow_missing_magnitudes=True)
    assert not ok and "vaccinated patients" in msg and "allow_sparse_evidence=true" in msg
    ok2, msg2 = agent_core.finalize_partial(str(out), allow_missing_magnitudes=True,
                                            allow_sparse_evidence=True)
    assert ok2 and "OVERRIDES USED" in msg2
    rec = json.loads(out.read_text())
    # the immunogenic elispot row also has no magnitude -> both genuinely-suppressed guards recorded
    assert "allow_sparse_evidence" in rec["finalize_overrides_used"]


# ---- audit trail: a CLEAN finalize records no overrides (v2.11.3) ----

def test_clean_finalize_records_no_overrides(tmp_path):
    out = tmp_path / "r.json"
    agent_core.init_partial(str(out), json.dumps({**META, "cohort_size": 1}))
    agent_core.append_section(str(out), "patients", json.dumps([_patient("P1", 2)]))
    agent_core.append_section(str(out), "immunizing_peptides", json.dumps([_imp("i1"), _imp("i2")]))
    ok, msg = agent_core.finalize_partial(str(out), allow_missing_magnitudes=True)
    assert ok and "OVERRIDES USED" not in msg
    assert json.loads(out.read_text()).get("finalize_overrides_used", []) == []


# ---- companion paper finalizes CLEAN without needing the peptide-recall override (v2.11.4) ----

def test_companion_paper_finalizes_clean_no_override(tmp_path):
    # 39972124-shape: paper cites 5 selected but lists only the 2 neoepitopes it functionally re-tested;
    # companion_paper_ref declares the manifest lives in the prior paper -> recall anchor relaxed, so
    # the record finalizes CLEAN (extracted lane) with NO override, unlike the same record without the ref.
    out = tmp_path / "r.json"
    meta = {**META, "cohort_size": 1, "n_selected_reported": 5,
            "companion_paper_ref": "Rojas et al. Nature 2023; PMID 37165196"}
    agent_core.init_partial(str(out), json.dumps(meta))
    agent_core.append_section(str(out), "patients", json.dumps([_patient("P1", 2)]))
    agent_core.append_section(str(out), "immunizing_peptides", json.dumps([_imp("i1"), _imp("i2")]))
    ok, msg = agent_core.finalize_partial(str(out), allow_missing_magnitudes=True)
    assert ok and "OVERRIDES USED" not in msg
    rec = json.loads(out.read_text())
    assert rec["finalize_overrides_used"] == []
    assert rec["companion_paper_ref"].endswith("37165196")


# ---- the MCP tool forwards both new flags (the P19 tool-path lesson) ----

def test_tool_forwards_both_new_flags():
    import inspect
    params = [p for p in inspect.signature(agent_core.finalize_partial).parameters]
    assert "allow_peptide_count_mismatch" in params and "allow_sparse_evidence" in params
    src = (PKT / "vaxtract" / "extraction_agent.py").read_text()
    body = src[src.index('@tool("finalize"'):src.index("async def finalize")] + \
        src[src.index("async def finalize"):src.index("return", src.index("async def finalize"))]
    assert "allow_peptide_count_mismatch" in body and "allow_sparse_evidence" in body


# ---- override TIERING (v2.11.5): soft-only stays clean; any hard/unknown -> needs_review ----

def test_overrides_soft_only_true_for_known_soft():
    assert agent_core.overrides_are_soft_only(["allow_unknown_funnel_size", "allow_regimen_divergence"])

def test_overrides_soft_only_false_when_any_hard():
    # the 34903219 / 33064988 shapes: a real recall gap mixed with a soft one -> still needs_review
    assert not agent_core.overrides_are_soft_only(["allow_regimen_divergence", "allow_sparse_evidence"])
    assert not agent_core.overrides_are_soft_only(["allow_peptide_count_mismatch"])

def test_overrides_soft_only_false_for_empty():
    assert not agent_core.overrides_are_soft_only([])
    assert not agent_core.overrides_are_soft_only(None)

def test_overrides_soft_only_false_for_unknown():
    # conservative: an unrecognized override defaults to HARD -> needs_review
    assert not agent_core.overrides_are_soft_only(["allow_some_future_guard"])

def test_every_finalize_override_flag_is_tiered_exactly_once():
    # invariant: a future allow_* flag MUST be classified HARD or SOFT (this gate fails otherwise).
    import inspect
    flags = {p for p in inspect.signature(agent_core.finalize_partial).parameters if p.startswith("allow_")}
    assert agent_core.HARD_OVERRIDES.isdisjoint(agent_core.SOFT_OVERRIDES)
    assert flags == (agent_core.HARD_OVERRIDES | agent_core.SOFT_OVERRIDES)
