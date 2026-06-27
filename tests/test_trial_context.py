"""Trial-context axes (schema v2.15): concomitant_therapy, safety_summary, tmb/msi.

Mirrors the repo convention (conftest puts PKT + PKT/cancervac_packet on sys.path, so
`import schema` / `import agent_core` resolve). ExtractedPatient inherits _Extracted ->
quoted_text + section_ref are REQUIRED, so the helpers supply them.
"""
import json
import pathlib

import pytest
from pydantic import ValidationError

import schema
from schema import (
    SCHEMA_VERSION, ExtractedPatient, ExtractedPaper,
    ConcomitantTherapy, SafetySummary,
)
import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]


def _patient(**kw):
    base = dict(
        quoted_text="patient P1", section_ref="Methods", paper_local_id="P1",
        indication="PDAC", vaccine_platform="rna",
        n_peptides_synthesized=10, n_peptides_immunogenic=2,
    )
    base.update(kw)
    return ExtractedPatient(**base)


def _paper(patients=None, **kw):
    base = dict(
        pmid="37165196", journal="Nature", year=2023, title="t",
        cohort_size=1, indication_summary="PDAC",
        patients=patients if patients is not None else [_patient()],
    )
    base.update(kw)
    return ExtractedPaper(**base)


def test_version_bumped():
    assert SCHEMA_VERSION == "2.16.0"


# --- vocab lockstep (the import-time guard already enforces this; assert membership too) ---
def test_concomitant_drug_class_vocab():
    import vocab
    assert set(vocab.CONCOMITANT_DRUG_CLASSES) == {
        "checkpoint_inhibitor", "chemotherapy", "radiotherapy", "targeted", "other"}
    assert set(vocab.CONCOMITANT_THERAPY_TIMING) == {"concurrent", "sequential", "unknown"}
    assert set(vocab.MSI_STATUSES) == {"mss", "msi_high", "unknown"}


# --- axis #1: concomitant_therapy ---
def test_concomitant_therapy_accepts_valid():
    p = _patient(concomitant_therapy=[
        ConcomitantTherapy(drug_class="checkpoint_inhibitor", agent="atezolizumab",
                           timing="concurrent"),
        ConcomitantTherapy(drug_class="chemotherapy", agent="mFOLFIRINOX", timing="sequential"),
    ])
    assert len(p.concomitant_therapy) == 2
    assert p.concomitant_therapy[0].drug_class == "checkpoint_inhibitor"


def test_concomitant_therapy_default_empty():
    assert _patient().concomitant_therapy == []


def test_concomitant_drug_class_rejects_unknown():
    with pytest.raises(ValidationError):
        ConcomitantTherapy(drug_class="immunotherapy")  # not in the Literal


# --- axis #2: safety_summary ---
def test_safety_summary_accepts_valid():
    s = SafetySummary(max_related_grade=3, any_grade3plus_related=True,
                      n_patients_with_related_ae=5, irae_present=True, raw="grade 3 in 5 pts")
    assert _paper(safety_summary=s).safety_summary.max_related_grade == 3


def test_safety_summary_default_none():
    assert _paper().safety_summary is None


def test_safety_grade_range_enforced():
    with pytest.raises(ValidationError):
        SafetySummary(max_related_grade=6)
    with pytest.raises(ValidationError):
        SafetySummary(max_related_grade=0)


def test_safety_grade_vs_grade3plus_contradiction():
    with pytest.raises(ValidationError):
        SafetySummary(max_related_grade=4, any_grade3plus_related=False)


# --- axis #3: tmb / msi ---
def test_tmb_msi_accepts_valid():
    p = _patient(tmb_value=12.0, tmb_raw="12 mut/Mb", msi_status="mss")
    assert p.tmb_value == 12.0 and p.msi_status == "mss"


def test_tmb_msi_defaults_none():
    p = _patient()
    assert p.tmb_value is None and p.tmb_raw is None and p.msi_status is None


def test_tmb_negative_rejected():
    with pytest.raises(ValidationError):
        _patient(tmb_value=-1)


def test_msi_status_rejects_unknown_value():
    with pytest.raises(ValidationError):
        _patient(msi_status="microsatellite_stable")


# --- the agent write-path for the paper-level safety_summary scalar (root-cause fix) ---
def _init_partial(tmp_path):
    out = str(tmp_path / "rec.json")
    meta = dict(pmid="99999999", journal="J", year=2024, title="t",
                cohort_size=1, indication_summary="x")
    ok, msg = agent_core.init_partial(out, json.dumps(meta))
    assert ok, msg
    return out


def _partial(out):
    import pathlib
    return json.loads(pathlib.Path(agent_core._partial_path(out)).read_text())


def test_safety_summary_not_settable_via_add_entities(tmp_path):
    # it is a scalar, not a SECTION_MODEL list -> add_entities must reject it
    out = _init_partial(tmp_path)
    ok, msg = agent_core.append_section(out, "safety_summary", "[{}]")
    assert not ok and "unknown section" in msg


def test_set_safety_summary_writes_and_validates(tmp_path):
    out = _init_partial(tmp_path)
    ok, msg = agent_core.set_safety_summary(
        out, json.dumps({"max_related_grade": 3, "any_grade3plus_related": True,
                         "n_patients_with_related_ae": 5, "raw": "grade 3 in 5 pts"}))
    assert ok, msg
    assert _partial(out)["safety_summary"]["max_related_grade"] == 3


def test_set_safety_summary_enforces_subschema(tmp_path):
    out = _init_partial(tmp_path)
    assert not agent_core.set_safety_summary(out, json.dumps({"max_related_grade": 9}))[0]      # >5
    # the grade<->grade3plus contradiction guard fires through the setter too
    assert not agent_core.set_safety_summary(
        out, json.dumps({"max_related_grade": 4, "any_grade3plus_related": False}))[0]


def test_set_safety_summary_clear(tmp_path):
    out = _init_partial(tmp_path)
    agent_core.set_safety_summary(out, json.dumps({"max_related_grade": 2}))
    ok, _ = agent_core.set_safety_summary(out, "null")
    assert ok and _partial(out)["safety_summary"] is None


# --- back-compat: the 4 signed-off reference records still validate ---
@pytest.mark.parametrize("name", [
    "keskin_extracted", "rojas_extracted", "li_extracted", "rojas_extracted.refresh"])
def test_reference_records_still_pass(name):
    path = PKT / "reference_records" / f"{name}.json"
    ok, msg = agent_core.validate_record(path.read_text())
    assert ok, f"{name} no longer validates under v2.15: {msg}"
