"""v2.14 #4: ScreeningReadout — the screening bucket. Entity validation + exactly-one-target +
finalize precedence (named evidence supersedes a screening row) + section wiring + back-compat."""
import json
import pathlib
from typing import get_args

import pytest
from pydantic import ValidationError

import agent_core
from cancervac_packet.schema import (
    ScreeningReadout, ManifestOutcome, ExtractedPaper, SCHEMA_VERSION,
)

PKT = pathlib.Path(__file__).resolve().parents[1]
META = {"pmid": "12345678", "journal": "Test J", "year": 2024, "title": "screening test",
        "cohort_size": 1, "indication_summary": "melanoma"}


def _scr(**kw):
    base = dict(quoted_text="No response", section_ref="MOESM4", assay="elispot",
                manifest_outcome="no_response", target_kind="candidate", candidate_paper_id="C1")
    base.update(kw)
    return ScreeningReadout(**base)


def test_schema_version():
    assert SCHEMA_VERSION == "2.16.0"

def test_manifest_outcome_enum():
    assert set(get_args(ManifestOutcome)) == {"response", "no_response", "not_evaluable"}

def test_screening_roundtrips():
    s = _scr(patient_paper_id="P25")
    assert s.manifest_outcome == "no_response" and s.candidate_paper_id == "C1"

def test_patient_optional():
    assert _scr().patient_paper_id is None   # a cohort-level manifest may have no per-patient grain

def test_exactly_one_target_enforced():
    with pytest.raises(ValidationError):
        _scr(epitope_paper_id="E1")          # candidate target_kind but TWO ids set
    with pytest.raises(ValidationError):
        ScreeningReadout(quoted_text="x", section_ref="T", assay="elispot",
                         manifest_outcome="response", target_kind="candidate")  # no id set

def test_section_model_wired():
    assert agent_core.SECTION_MODEL.get("screening_readouts") == "ScreeningReadout"

def test_paper_screening_default_empty_and_backcompat():
    rec = ExtractedPaper(**json.loads((PKT / "reference_records" / "li_extracted.json").read_text()))
    assert rec.screening_readouts == []   # legacy record validates with the new optional list


# ---- finalize PRECEDENCE: a named evidence row supersedes a screening row at the same key ----

def test_dedup_drops_screening_covered_by_evidence():
    rec = {"evidence": [{"patient_paper_id": "P25", "candidate_paper_id": "C1", "assay": "elispot",
                         "outcome": "immunogenic"}],
           "screening_readouts": [
               {"patient_paper_id": "P25", "candidate_paper_id": "C1", "assay": "elispot",
                "manifest_outcome": "no_response"},                              # collides -> dropped
               {"patient_paper_id": "P25", "candidate_paper_id": "C2", "assay": "elispot",
                "manifest_outcome": "no_response"}]}                              # distinct -> kept
    dropped = agent_core._drop_screening_covered_by_evidence(rec)
    assert dropped == 1
    assert [s["candidate_paper_id"] for s in rec["screening_readouts"]] == ["C2"]

def test_dedup_noop_when_no_overlap():
    rec = {"evidence": [{"patient_paper_id": "P1", "candidate_paper_id": "C9", "assay": "elispot",
                         "outcome": "immunogenic"}],
           "screening_readouts": [{"patient_paper_id": "P25", "candidate_paper_id": "C1",
                                   "assay": "elispot", "manifest_outcome": "no_response"}]}
    assert agent_core._drop_screening_covered_by_evidence(rec) == 0
    assert len(rec["screening_readouts"]) == 1
