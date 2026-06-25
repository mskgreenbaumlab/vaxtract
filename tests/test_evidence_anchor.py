"""v2.11.1 P20.1 — evidence-completeness anchors. The fix for the run-to-run evidence-count
variance (diagnosed as enumeration-granularity, not duplicates): paper-stated n_immunogenic_reported
/ n_tested_negative_reported + a soft finalize nudge. Covers the gap detector (duck-typed), the
finalize block-once + override, and that the override flag is forwarded by the MCP tool."""
import json
import pathlib
from types import SimpleNamespace as NS

import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]
META = {"pmid": "12345678", "journal": "Test J", "year": 2024, "title": "Anchor test",
        "cohort_size": 1, "indication_summary": "melanoma"}


def _ev(pid, outcome, tgt="i1", tk="immunizing_peptide"):
    d = dict(patient_paper_id=pid, target_kind=tk, outcome=outcome,
             immunizing_peptide_paper_id=None, epitope_paper_id=None,
             pool_paper_id=None, candidate_paper_id=None)
    d[{"immunizing_peptide": "immunizing_peptide_paper_id", "epitope": "epitope_paper_id",
       "pool": "pool_paper_id", "candidate": "candidate_paper_id"}[tk]] = tgt
    return NS(**d)


# ---- gap detector (duck-typed) ----

def test_no_anchor_no_nudge():
    rec = NS(evidence=[_ev("P1", "immunogenic")], n_immunogenic_reported=None,
             n_tested_negative_reported=None)
    assert agent_core._evidence_anchor_gap(rec) == []

def test_immunogenic_recall_gap_flagged():
    rec = NS(evidence=[_ev("P1", "immunogenic", "i1")], n_immunogenic_reported=5,
             n_tested_negative_reported=None)
    msgs = agent_core._evidence_anchor_gap(rec)
    assert msgs and "MISSED" in msgs[0] and "5" in msgs[0]

def test_immunogenic_match_is_clean():
    rec = NS(evidence=[_ev("P1", "immunogenic", "i1"), _ev("P2", "immunogenic", "i2")],
             n_immunogenic_reported=2, n_tested_negative_reported=None)
    assert agent_core._evidence_anchor_gap(rec) == []

def test_over_recorded_immunogenic_not_flagged():  # pooled/extra is fine on the immunogenic side
    rec = NS(evidence=[_ev("P1", "immunogenic", "i1"), _ev("P2", "immunogenic", "i2")],
             n_immunogenic_reported=1, n_tested_negative_reported=None)
    assert agent_core._evidence_anchor_gap(rec) == []

def test_negative_over_enumeration_flagged():  # the grain blowup (a negative per blank cell)
    rec = NS(evidence=[_ev("P1", "not_immunogenic", "i1"), _ev("P1", "not_immunogenic", "i2"),
                       _ev("P1", "not_immunogenic", "i3")],
             n_immunogenic_reported=None, n_tested_negative_reported=1)  # 3 > 1.5x
    msgs = agent_core._evidence_anchor_gap(rec)
    assert msgs and "OVER-enumerated" in msgs[0]

def test_negative_under_recorded_is_clean():
    # v2.11.5: UNDER is the canonical grain rule (record only NAMED negatives, fewer than a cohort
    # total) -> must NOT fire. Was the dominant FALSE nudge on the iris batch (Rojas 0 vs 8 etc.).
    rec = NS(evidence=[_ev("P1", "not_immunogenic", "i1")],
             n_immunogenic_reported=None, n_tested_negative_reported=4)
    assert agent_core._evidence_anchor_gap(rec) == []

def test_negative_small_over_within_tolerance_clean():
    # +1 over (4 vs 3, 33%) is noise, below the 1.5x threshold -> clean (the live Keskin 4-vs-3 case).
    rec = NS(evidence=[_ev("P1", "not_immunogenic", f"i{i}") for i in range(4)],
             n_immunogenic_reported=None, n_tested_negative_reported=3)
    assert agent_core._evidence_anchor_gap(rec) == []

def test_negative_exact_match_clean():
    rec = NS(evidence=[_ev("P1", "not_immunogenic", "i1"), _ev("P2", "not_immunogenic", "i2")],
             n_immunogenic_reported=None, n_tested_negative_reported=2)
    assert agent_core._evidence_anchor_gap(rec) == []

def test_immunogenic_pool_expansion_counts_members():
    # v2.11.5: ONE immunogenic pool row stands in for its member peptides; the paper counts members
    # toward n_immunogenic_reported, so expand by member count. Distinct-target would be 1 (falsely
    # fire vs 7); member-expanded is 7 -> clean.
    pool = NS(paper_local_id="pl1", member_peptide_ids=[f"i{i}" for i in range(7)])
    rec = NS(evidence=[_ev("P1", "immunogenic", tgt="pl1", tk="pool")], pools=[pool],
             n_immunogenic_reported=7, n_tested_negative_reported=None)
    assert agent_core._evidence_anchor_gap(rec) == []

def test_companion_exempts_both_anchor_sides():
    # a companion paper defers BOTH immunogenic and negative characterization -> neither side fires,
    # even though without companion both would (imm 1<25 AND neg 9 > 1.5x*2).
    rec = NS(evidence=[_ev("P1", "immunogenic", "i1")]
                      + [_ev("P1", "not_immunogenic", f"n{i}") for i in range(9)],
             pools=[], n_immunogenic_reported=25, n_tested_negative_reported=2,
             companion_paper_ref="Rojas et al. 2023; PMID 37165196")
    assert agent_core._evidence_anchor_gap(rec) == []

def test_negative_outcome_does_not_count_against_tested_negative():
    # 'negative' (e.g. tumour non-recognition) is NOT a tested-non-immunogenic row -> excluded.
    # 4 not_immunogenic + 2 negative, paper reports 4 -> clean (the live-Keskin false-override fix).
    rec = NS(evidence=[_ev("P1", "not_immunogenic", f"i{i}") for i in range(4)]
                      + [_ev("P1", "negative", "t1"), _ev("P2", "negative", "t2")],
             n_immunogenic_reported=None, n_tested_negative_reported=4)
    assert agent_core._evidence_anchor_gap(rec) == []


# ---- finalize blocks once, then override proceeds ----

def _patient(pid="P1"):
    return {"quoted_text": f"pt {pid}", "section_ref": "M", "paper_local_id": pid,
            "indication": "PDAC", "n_peptides_synthesized": 2, "n_peptides_immunogenic": 1}

def _imp(pid="i1"):
    return {"paper_local_id": pid, "sequence": "SLLQHLIGL", "is_neoantigen": True,
            "patient_paper_id": "P1", "quoted_text": "q", "section_ref": "s"}

def _evd(outcome="not_immunogenic", tgt="i1"):
    return {"patient_paper_id": "P1", "target_kind": "immunizing_peptide",
            "immunizing_peptide_paper_id": tgt, "assay": "elispot", "outcome": outcome,
            "quoted_text": "q", "section_ref": "s"}

def test_finalize_blocks_then_overrides_on_anchor_gap(tmp_path):
    out = tmp_path / "r.json"
    meta = {**META, "n_tested_negative_reported": 1}   # paper says 1 negative
    agent_core.init_partial(str(out), json.dumps(meta))
    agent_core.append_section(str(out), "patients", json.dumps([_patient()]))
    agent_core.append_section(str(out), "immunizing_peptides", json.dumps([_imp("i1"), _imp("i2")]))
    # record TWO negatives -> disagrees with the stated 1
    agent_core.append_section(str(out), "evidence",
                              json.dumps([_evd(tgt="i1"), _evd(tgt="i2")]))
    ok, msg = agent_core.finalize_partial(str(out), allow_missing_magnitudes=True)
    assert not ok and "stated counts" in msg and "allow_evidence_count_mismatch=true" in msg
    ok2, _ = agent_core.finalize_partial(str(out), allow_missing_magnitudes=True,
                                         allow_evidence_count_mismatch=True)
    assert ok2


def test_tool_forwards_the_new_flag():
    import inspect
    params = [p for p in inspect.signature(agent_core.finalize_partial).parameters if p.startswith("allow_")]
    src = (PKT / "vaxtract" / "extraction_agent.py").read_text()
    body = src[src.index('@tool("finalize"'):src.index("@tool(", src.index('@tool("finalize"') + 10)]
    assert "allow_evidence_count_mismatch" in params and "allow_evidence_count_mismatch" in body
