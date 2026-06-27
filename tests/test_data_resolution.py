"""v2.16 data_resolution: admit genuinely coarse-grained papers faithfully.

Covers derive_data_resolution (achieved grain), _is_faithfully_coarse (the derived + anchor-checked
gate), and that the per-sequence completeness anchors (_peptide_recall_gap, _evidence_breadth_gap)
are exempted for a faithfully-coarse paper but NOT for under-extraction or legacy records.
"""
from types import SimpleNamespace as NS

import agent_core


def _peptide(seq="MELANOMA"): return NS(sequence=seq)


def _rec(**kw):
    base = dict(
        immunizing_peptides=[], epitopes=[], neoantigen_mutations=[], screening_readouts=[],
        evidence=[], survival_outcomes=[], safety_summary=None, patients=[],
        n_immunogenic_reported=None, n_selected_reported=None,
        companion_paper_ref=None, data_resolution=None, peptide_manifest_present=None,
    )
    base.update(kw)
    return NS(**base)


# ----- derive_data_resolution (achieved grain) --------------------------------
def test_derive_per_sequence():
    assert agent_core.derive_data_resolution(_rec(immunizing_peptides=[_peptide()])) == "per_sequence"

def test_derive_per_mutation():
    assert agent_core.derive_data_resolution(_rec(neoantigen_mutations=[NS()])) == "per_mutation"

def test_derive_cohort_summary():
    assert agent_core.derive_data_resolution(_rec(n_immunogenic_reported=12)) == "cohort_summary"

def test_derive_clinical_only():
    assert agent_core.derive_data_resolution(_rec(survival_outcomes=[NS()])) == "clinical_only"

def test_derive_empty_is_none():
    assert agent_core.derive_data_resolution(_rec()) is None

def test_finest_layer_wins():
    # a per-mutation paper that also has survival is still per_mutation (finer wins)
    assert agent_core.derive_data_resolution(
        _rec(neoantigen_mutations=[NS()], survival_outcomes=[NS()])) == "per_mutation"


# ----- _is_faithfully_coarse (derived + anchor) -------------------------------
def test_coarse_without_manifest_is_faithfully_coarse():
    assert agent_core._is_faithfully_coarse(_rec(neoantigen_mutations=[NS()])) is True

def test_coarse_with_manifest_is_NOT_faithful():
    # a manifest WAS available -> a coarse result is under-extraction, not faithful
    assert agent_core._is_faithfully_coarse(
        _rec(neoantigen_mutations=[NS()], peptide_manifest_present=True)) is False

def test_per_sequence_is_not_coarse():
    assert agent_core._is_faithfully_coarse(_rec(immunizing_peptides=[_peptide()])) is False

def test_empty_record_is_not_coarse():
    assert agent_core._is_faithfully_coarse(_rec()) is False


# ----- guards exempt a faithfully-coarse paper --------------------------------
def test_recall_gap_exempts_coarse_paper():
    # 39762422 shape: 0 peptides, 40 mutations, n_selected_reported a cited count, no manifest
    rec = _rec(neoantigen_mutations=[NS() for _ in range(40)], n_selected_reported=1404)
    assert agent_core._peptide_recall_gap(rec) is None

def test_recall_gap_still_fires_when_manifest_present():
    rec = _rec(neoantigen_mutations=[NS() for _ in range(40)], n_selected_reported=1404,
               peptide_manifest_present=True)
    assert agent_core._peptide_recall_gap(rec) is not None

def test_recall_gap_still_fires_for_legacy_thin_per_sequence():
    # legacy per-sequence record (5 peptides) far below the stated 100 -> genuine miss, still flags
    rec = _rec(immunizing_peptides=[_peptide(f"PEP{i}") for i in range(5)], n_selected_reported=100)
    assert agent_core._peptide_recall_gap(rec) is not None

def test_breadth_gap_exempts_coarse_paper():
    pts = [NS(paper_local_id=f"P{i}", n_peptides_synthesized=10, cohort_kind="patient")
           for i in range(7)]
    rec = _rec(neoantigen_mutations=[NS() for _ in range(40)], patients=pts, evidence=[])
    assert agent_core._evidence_breadth_gap(rec) is None


# ----- prompt wiring: the agent is told to declare the grain (both prompts, in lockstep) ------
def test_data_resolution_guidance_in_both_prompts():
    import pathlib
    pkt = pathlib.Path(__file__).resolve().parents[1]
    pr = (pkt / "vaxtract" / "prompt_render.py").read_text()
    rules = (pkt / "lite_extract" / "RULES.md").read_text()
    for blob in (pr, rules):
        assert "data_resolution" in blob
        assert "peptide_manifest_present" in blob
