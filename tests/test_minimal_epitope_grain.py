"""v2.12 P23 follow-up: _minimal_epitope_grain_gap — class-I long-peptide mislabel + dropped minimal layer.

Calibration is pinned: ZERO-fire on the 3 gold refs (whose class-I epitopes are 8-12mers with minimal
epitopes present) and FIRE on the lite-lane Li artifact (16 epitopes, all 19-29mers, 0 minimal).
"""
import json
import pathlib
from types import SimpleNamespace as NS

import pytest
from pydantic import ValidationError
import agent_core
from cancervac_packet.schema import ExtractedPaper, MinimalEpitope, Measurement

PKT = pathlib.Path(__file__).resolve().parents[1]


def _epi(seq, mhc="I", aff=False, quoted="x", section="T"):
    return NS(sequence=seq, mhc_class=mhc,
              predicted_affinity=({"unit": "unknown"} if aff else None),
              quoted_text=quoted, section_ref=section)


# ---- side (a) is now a HARD schema reject (v2.13.1, Fable #1): class-I >14 aa never constructs ----

def test_oversized_class_i_hard_rejected():
    with pytest.raises(ValidationError, match="minimal binder"):
        MinimalEpitope(quoted_text="x", section_ref="T", paper_local_id="E1",
                       sequence="GILARNLVPMVATVQGQNLK", is_neoantigen=True, mhc_class="I",
                       predicted_affinity=Measurement(unit="unknown", raw="n/a"))

def test_class_i_at_or_below_14_ok():
    # a 12mer class-I (Rojas's longest) constructs fine; the soft grain nudge (side b) also stays clean
    e = MinimalEpitope(quoted_text="x", section_ref="T", paper_local_id="E1",
                       sequence="SLLQHLIGLAAA", is_neoantigen=True, mhc_class="I",
                       predicted_affinity=Measurement(unit="unknown", raw="n/a"))
    assert e.mhc_class == "I"
    rec = NS(epitopes=[_epi("SLLQHLIGLAAA")], immunizing_peptides=[1] * 6, companion_paper_ref=None)
    assert agent_core._minimal_epitope_grain_gap(rec) == []


# ---- side (b): prediction cue present but the minimal-epitope layer was dropped ----

def test_dropped_minimal_layer_flagged():
    # the lite-Li shape: affinity cues on long mislabeled epitopes, 0 class-I 8-11mer minimal, many IMP
    rec = NS(epitopes=[_epi("GILARNLVPMVATVQGQNLK", aff=True) for _ in range(12)],
             immunizing_peptides=[1] * 50, companion_paper_ref=None)
    g = agent_core._minimal_epitope_grain_gap(rec)
    assert g and any("minimal-epitope layer" in m for m in g)

def test_dropped_layer_via_text_cue_flagged():
    # cue lives in section_ref (NetMHC), epitope itself class-II so not counted as class-I minimal
    rec = NS(epitopes=[_epi("AAAAAAAAAAAAAAA", mhc="II", section="NetMHCpan prediction")],
             immunizing_peptides=[1] * 6, companion_paper_ref=None)
    assert agent_core._minimal_epitope_grain_gap(rec)

def test_companion_exempts_dropped_layer():
    # same dropped-layer shape but the manifest is deferred to a companion paper -> side (b) exempt
    rec = NS(epitopes=[_epi("AAAAAAAAAAAAAAA", mhc="II", section="NetMHCpan prediction")],
             immunizing_peptides=[1] * 6, companion_paper_ref="PMID 37165196")
    assert agent_core._minimal_epitope_grain_gap(rec) == []

def test_no_cue_no_imp_clean():
    # no prediction cue and few IMP -> not the dropped-layer failure mode
    rec = NS(epitopes=[], immunizing_peptides=[1] * 3, companion_paper_ref=None)
    assert agent_core._minimal_epitope_grain_gap(rec) == []

def test_clean_when_minimal_present():
    # a real 8mer minimal class-I epitope -> no nudge even with a prediction cue
    rec = NS(epitopes=[_epi("SIINFEKL", aff=True)], immunizing_peptides=[1] * 6, companion_paper_ref=None)
    assert agent_core._minimal_epitope_grain_gap(rec) == []


# ---- real-data calibration: 0-fire on gold, fire on the lite-Li artifact ----

@pytest.mark.parametrize("name", ["keskin", "rojas", "li"])
def test_gold_refs_do_not_fire(name):
    rec = ExtractedPaper(**json.loads((PKT / "reference_records" / f"{name}_extracted.json").read_text()))
    assert agent_core._minimal_epitope_grain_gap(rec) == []

def test_lite_li_artifact_now_hard_rejected():
    # v2.13.1: the lite-Li record (12 class-I "epitopes" of 19-29 aa) is now INVALID outright — the
    # mislabeled long peptides are rejected at construction, a stronger guarantee than the old soft nudge.
    p = PKT / "outputs" / "classii_compare_2026-06-11" / "lite" / "33879241.json"
    if not p.exists():
        pytest.skip("lite-Li comparison artifact not present")
    with pytest.raises(ValidationError):
        ExtractedPaper(**json.loads(p.read_text()))


# ---- override tiering ----

def test_override_is_hard():
    assert "allow_missing_minimal_epitopes" in agent_core.HARD_OVERRIDES
    assert not agent_core.overrides_are_soft_only(["allow_missing_minimal_epitopes"])
