"""v2.10 P22 — response→outcome linkage. Covers (Move 1) per-patient `trial_setting` and
(Move 2) `ClinicalBenefitSignal`: per-object guards, per-patient + paper-level placement, the
paper-level SECTION_MODEL wiring + appendability, and back-compat with the signed-off fixtures."""
import json
import pathlib

import pytest
from pydantic import ValidationError

import agent_core
import schema
from schema import ExtractedPaper, ExtractedPatient, ClinicalBenefitSignal

PKT = pathlib.Path(__file__).resolve().parents[1]

META = {
    "pmid": "12345678", "journal": "Test J", "year": 2024, "title": "Outcome-linkage test",
    "cohort_size": 1, "indication_summary": "melanoma",
}


def _patient(pid="P1", **kw):
    base = {
        "quoted_text": f"patient {pid}", "section_ref": "Methods", "paper_local_id": pid,
        "indication": "PDAC", "n_peptides_synthesized": 0, "n_peptides_immunogenic": 0,
    }
    base.update(kw)
    return base


# ---- Move 1: trial_setting (per-patient) ----

def test_trial_setting_valid_values():
    for s in ("adjuvant", "neoadjuvant", "perioperative", "metastatic",
              "locally_advanced", "recurrent", "other"):
        p = ExtractedPatient(**_patient(trial_setting=s))
        assert p.trial_setting == s

def test_trial_setting_defaults_none():
    assert ExtractedPatient(**_patient()).trial_setting is None

def test_trial_setting_off_vocab_rejected():
    with pytest.raises(ValidationError):
        ExtractedPatient(**_patient(trial_setting="phase_2"))

def test_trial_setting_is_per_patient():
    # a combined-cohort paper carries different settings per subject (not collapsed)
    rec = ExtractedPaper(**{**META, "cohort_size": 2}, patients=[
        ExtractedPatient(**_patient("P1", indication="melanoma", trial_setting="metastatic")),
        ExtractedPatient(**_patient("P2", indication="PDAC", trial_setting="adjuvant")),
    ])
    assert [p.trial_setting for p in rec.patients] == ["metastatic", "adjuvant"]


# ---- Move 2: ClinicalBenefitSignal (per-object guards) ----

def test_benefit_signal_roundtrips():
    b = ClinicalBenefitSignal(readout="tumor_infiltration", direction="detected",
                              timepoint_phase="post_vaccine", associated_with_response=True)
    assert b.readout == "tumor_infiltration" and b.associated_with_response is True

def test_benefit_signal_off_vocab_readout_rejected():
    with pytest.raises(ValidationError):
        ClinicalBenefitSignal(readout="vibes", direction="increased")

def test_benefit_cleared_incoherent_with_antigen_loss():
    with pytest.raises(ValidationError):
        ClinicalBenefitSignal(readout="antigen_loss", direction="cleared")
    with pytest.raises(ValidationError):
        ClinicalBenefitSignal(readout="epitope_spreading", direction="cleared")

def test_benefit_cleared_ok_on_ctdna():
    assert ClinicalBenefitSignal(readout="ctdna_dynamics", direction="cleared").direction == "cleared"

def test_antigen_loss_uses_lost_direction():
    b = ClinicalBenefitSignal(readout="antigen_loss", direction="lost", timepoint_label="at relapse")
    assert b.direction == "lost"


# ---- placement: per-patient AND paper-level ----

def test_benefit_signal_per_patient_and_paper_level():
    rec = ExtractedPaper(**META, patients=[
        ExtractedPatient(**_patient("P1"), clinical_benefit_signals=[
            ClinicalBenefitSignal(readout="ctdna_dynamics", direction="cleared")]),
    ], clinical_benefit_signals=[
        ClinicalBenefitSignal(readout="tumor_infiltration", direction="detected")])
    assert rec.patients[0].clinical_benefit_signals[0].readout == "ctdna_dynamics"
    assert rec.clinical_benefit_signals[0].readout == "tumor_infiltration"


# ---- paper-level SECTION_MODEL wiring + appendability ----

def test_benefit_signal_in_section_model():
    assert agent_core.SECTION_MODEL.get("clinical_benefit_signals") == "ClinicalBenefitSignal"

def test_benefit_signal_appendable_paper_level(tmp_path):
    out = tmp_path / "r.json"
    agent_core.init_partial(str(out), json.dumps(META))
    ok, msg = agent_core.append_section(str(out), "clinical_benefit_signals", json.dumps(
        [{"readout": "epitope_spreading", "direction": "increased", "timepoint_label": "post-boost"}]))
    assert ok, msg
    part = json.loads((tmp_path / "r.json.partial.json").read_text())
    assert part["clinical_benefit_signals"][0]["readout"] == "epitope_spreading"


# ---- back-compat ----

def test_li_reference_back_compat_empty_linkage_fields():
    # li was NOT refreshed -> proves a record with the P22 defaults (None / []) still validates.
    rec = ExtractedPaper(**json.loads((PKT / "reference_records" / "li_extracted.json").read_text()))
    assert all(p.trial_setting is None for p in rec.patients)
    assert all(p.clinical_benefit_signals == [] for p in rec.patients)
    assert rec.clinical_benefit_signals == []
    assert schema.SCHEMA_VERSION == "2.16.0"

def test_refreshed_references_carry_p22_fields():
    # rojas + keskin refreshed 2026-06-07: trial_setting on every patient + >=1 benefit signal.
    for n, n_sig in (("rojas", 1), ("keskin", 4)):
        rec = ExtractedPaper(**json.loads((PKT / "reference_records" / f"{n}_extracted.json").read_text()))
        assert all(p.trial_setting == "adjuvant" for p in rec.patients), n
        got = sum(len(p.clinical_benefit_signals) for p in rec.patients)
        assert got == n_sig, f"{n}: expected {n_sig} benefit signals, got {got}"
        assert all(b.readout == "tumor_infiltration" for p in rec.patients for b in p.clinical_benefit_signals)
