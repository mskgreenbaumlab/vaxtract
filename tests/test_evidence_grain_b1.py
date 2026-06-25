"""v2.13 B1: a prediction/target MANIFEST's per-row 'No response' is the screening COUNT, not a
per-row evidence negative. Guidance lives in both prompts; the _evidence_anchor_gap over-enumeration
backstop catches a future run that emits one not_immunogenic row per manifest target (the Rojas ~225 shape)."""
import pathlib
from types import SimpleNamespace as NS

import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]


def _eneg(pid, tgt):
    return NS(patient_paper_id=pid, outcome="not_immunogenic", immunizing_peptide_paper_id=tgt,
              epitope_paper_id=None, pool_paper_id=None, candidate_paper_id=None)


def test_per_target_negative_blowup_flagged():
    # Rojas shape: ~200 manifest "No response" rows recorded as not_immunogenic vs a small stated
    # named-negative count -> over-enumeration backstop fires (routes the record to needs_review).
    rec = NS(evidence=[_eneg(f"P{i}", f"i{i}") for i in range(200)],
             n_immunogenic_reported=None, n_tested_negative_reported=2, companion_paper_ref=None)
    gaps = agent_core._evidence_anchor_gap(rec)
    assert gaps and any("OVER-enumerated" in g for g in gaps)


def test_b1_guidance_present_in_both_prompts():
    rules = (PKT / "lite_extract" / "RULES.md").read_text()
    pr = (PKT / "vaxtract" / "prompt_render.py").read_text()
    assert "MANIFEST SCREENING READOUT" in rules
    assert "MANIFEST SCREENING READOUT" in pr
