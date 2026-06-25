"""v2.13 A1: ExtractedPatient.cohort_kind — enum, optionality, gold-Li tagging, and the
'vaccinated patient' recall denominator excluding methodological arms."""
import json
import pathlib
from types import SimpleNamespace as NS
from typing import get_args

import pytest
from pydantic import ValidationError

import agent_core
from cancervac_packet.schema import ExtractedPatient, CohortKind, SCHEMA_VERSION

PKT = pathlib.Path(__file__).resolve().parents[1]


def _patient(pid="P1", **kw):
    base = {"quoted_text": f"pt {pid}", "section_ref": "Methods", "paper_local_id": pid,
            "indication": "melanoma", "n_peptides_synthesized": 0, "n_peptides_immunogenic": 0}
    base.update(kw)
    return base


def test_schema_version():
    assert SCHEMA_VERSION == "2.15.0"

def test_cohort_kind_enum_members():
    assert set(get_args(CohortKind)) == {
        "patient", "tumor_model", "model_antigen_validation", "healthy_donor", "other"}

def test_cohort_kind_defaults_none():
    assert ExtractedPatient(**_patient()).cohort_kind is None

def test_cohort_kind_non_methodological_accepts_freely():
    # patient / tumor_model / other assert nothing extra -> no cue required
    for k in ("patient", "tumor_model", "other"):
        assert ExtractedPatient(**_patient(cohort_kind=k)).cohort_kind == k

def test_cohort_kind_off_vocab_rejected():
    with pytest.raises(ValidationError):
        ExtractedPatient(**_patient(cohort_kind="cell_line"))


# ---- v2.13.1 (Fable #6): a methodological tag must be quotable (it removes the cohort from recall) ----

def test_methodological_tag_needs_cue():
    with pytest.raises(ValidationError, match="cohort_kind"):
        ExtractedPatient(**_patient(cohort_kind="model_antigen_validation"))  # 'melanoma' has no cue

def test_methodological_tag_ok_with_cue():
    p = ExtractedPatient(**_patient(cohort_kind="model_antigen_validation",
                                    indication="HLA-A2 model-antigen optimization (gp100)"))
    assert p.cohort_kind == "model_antigen_validation"

def test_healthy_donor_tag_needs_cue():
    with pytest.raises(ValidationError):
        ExtractedPatient(**_patient(cohort_kind="healthy_donor"))
    ok = ExtractedPatient(**_patient(cohort_kind="healthy_donor",
                                     quoted_text="PBMC from healthy donors"))
    assert ok.cohort_kind == "healthy_donor"


# ---- gold-Li refresh: the methodological arm is tagged, validates ----

def test_li_gold_cohort_kinds_tagged():
    from cancervac_packet.schema import ExtractedPaper
    rec = ExtractedPaper(**json.loads((PKT / "reference_records" / "li_extracted.json").read_text()))
    by = {p.paper_local_id: p.cohort_kind for p in rec.patients}
    assert by["P1"] == "model_antigen_validation"          # HLA-A2 model-antigen validation arm
    assert by["P2"] == "tumor_model" and by["P3"] == "tumor_model"   # E0771 / 4T1.2
    assert by["P4"] == "patient"                            # human PNET


# ---- the 'vaccinated patient' denominator excludes methodological arms ----

def _pat(pid, nsyn, kind=None):
    return NS(paper_local_id=pid, n_peptides_synthesized=nsyn, cohort_kind=kind)

def test_breadth_excludes_methodological_arms():
    # 5 real vaccinated cohorts (<6 -> no nudge) + 3 methodological. If the methodological arms were
    # counted, the denominator would be 8 (>=6) and a 0-coverage record would FIRE. Excluding them -> None.
    rec = NS(patients=[_pat(f"P{i}", 1, "patient") for i in range(5)]
                      + [_pat(f"M{i}", 1, "model_antigen_validation") for i in range(3)],
             evidence=[])
    assert agent_core._evidence_breadth_gap(rec) is None

def test_breadth_counts_tumor_model_and_legacy_none():
    # tumor_model + legacy-None DO count as vaccinated cohorts -> 8 vaccinated, 1 covered (<1/3) -> fires
    rec = NS(patients=[_pat(f"T{i}", 1, "tumor_model") for i in range(4)]
                      + [_pat(f"L{i}", 1, None) for i in range(4)],
             evidence=[NS(patient_paper_id="T0")])
    assert agent_core._evidence_breadth_gap(rec) is not None
