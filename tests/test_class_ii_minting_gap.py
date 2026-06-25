"""v2.12 P23: _class_ii_minting_gap anchor — Keskin regression (CD4 cue, 0 class-II) + class-I false-pos guard."""
from types import SimpleNamespace as NS
import pytest
import agent_core


def _epi(mhc_class, quoted="x", hla=None):
    return NS(mhc_class=mhc_class, quoted_text=quoted, hla_allele=hla)

def _ev(quoted="x", mhc_class=None, subset=None):
    return NS(quoted_text=quoted, mhc_class=mhc_class, t_cell_subset=subset)


def test_keskin_regression_cd4_cue_zero_class_ii_flagged():
    # the 30568305 shape: a CD4 / class-II response in prose but NO class-II record minted
    rec = NS(epitopes=[_epi("I", "HLA-A response")],
             evidence=[_ev("CD4+ response to GPC1 long peptide")])
    gaps = agent_core._class_ii_minting_gap(rec)
    assert gaps and any("CD4" in g for g in gaps)

def test_named_dr_allele_zero_class_ii_flagged():
    rec = NS(epitopes=[_epi("I", "x")],
             evidence=[_ev("HLA-DRB1*04:01-restricted reactivity")])
    gaps = agent_core._class_ii_minting_gap(rec)
    assert gaps and any("DR/DP/DQ" in g for g in gaps)

def test_clean_when_class_ii_epitope_present():
    # Rojas shape: class-II epitopes ARE minted -> no nudge even with CD4 cues
    rec = NS(epitopes=[_epi("II", "HLA-DRB1*04:01", hla="HLA-DRB1*04:01")],
             evidence=[_ev("CD4+ response")])
    assert agent_core._class_ii_minting_gap(rec) == []

def test_clean_when_evidence_typed_class_ii():
    rec = NS(epitopes=[_epi("I", "x")],
             evidence=[_ev("CD4 response", mhc_class="class_ii")])
    assert agent_core._class_ii_minting_gap(rec) == []

def test_inferred_class_ii_epitope_does_not_silence_nudge():
    # v2.13.1 (Fable): a GUESSED (mhc_class_inferred) class-II epitope must NOT satisfy the anti-guess
    # guard — only a reported (non-inferred) class-II record counts.
    inferred = NS(mhc_class="II", mhc_class_inferred=True, quoted_text="length-inferred", hla_allele=None)
    rec = NS(epitopes=[inferred], evidence=[_ev("CD4+ response to GPC1 long peptide")])
    gaps = agent_core._class_ii_minting_gap(rec)
    assert gaps and any("CD4" in g for g in gaps)

def test_class_i_only_paper_not_flagged():
    # FALSE-POSITIVE GUARD: a genuinely class-I-only paper (no CD4/class-II cue anywhere) -> no nudge
    rec = NS(epitopes=[_epi("I", "HLA-A*02:01 predicted 45 nM", hla="HLA-A*02:01")],
             evidence=[_ev("CD8+ cytotoxic response, HLA-A*02:01")])
    assert agent_core._class_ii_minting_gap(rec) == []

def test_override_is_hard():
    assert "allow_missing_class_ii" in agent_core.HARD_OVERRIDES
    assert not agent_core.overrides_are_soft_only(["allow_missing_class_ii"])

@pytest.mark.skip(reason="LIVE-ONLY: re-extract Keskin 30568305 -> must mint >=1 class_ii MinimalEpitope "
                         "(GPC1/SHANK2/SVEP1 CD4 hits); Rojas 37165196 stays clean. Run under RUN_LIVE.")
def test_live_keskin_mints_class_ii():
    pass
