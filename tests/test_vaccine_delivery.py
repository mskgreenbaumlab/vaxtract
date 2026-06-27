"""v2.9 P21 VaccineDelivery — per-patient delivery covariates + the cross-patient
regimen-consistency nudge. Confirms the per-object guards, the soft nudge logic, that the
nudge is reachable + overridable through finalize_partial, and back-compat."""
import json
import pathlib
from types import SimpleNamespace as NS

import pytest
from pydantic import ValidationError

import agent_core
import schema
from schema import ExtractedPaper, VaccineDelivery

PKT = pathlib.Path(__file__).resolve().parents[1]

META = {
    "pmid": "12345678", "journal": "Test J", "year": 2024, "title": "Delivery test",
    "cohort_size": 2, "indication_summary": "melanoma",
}


def _patient(pid, platform="synthetic_long_peptide", delivery=None):
    p = {
        "quoted_text": f"patient {pid}", "section_ref": "Methods", "paper_local_id": pid,
        "indication": "melanoma", "vaccine_platform": platform,
        "n_peptides_synthesized": 0, "n_peptides_immunogenic": 0,  # 0 -> no peptide-count guard in this delivery test
    }
    if delivery is not None:
        p["vaccine_delivery"] = delivery
    return p


POLY = {"adjuvant": "poly_iclc", "adjuvant_detail": "poly-ICLC 1 mg", "n_priming_doses": 5,
        "n_boost_doses": 2, "dose_per_peptide_ug": 300, "dose_basis": "per_peptide"}


# ---- per-object guards (schema) ----

def test_valid_delivery_roundtrips():
    vd = VaccineDelivery(**POLY, dose_amount_raw="300 µg/peptide", n_doses_received=7,
                         weeks_surgery_to_first_dose=18.5)
    assert vd.adjuvant == "poly_iclc" and vd.dose_per_peptide_ug == 300.0

def test_received_exceeds_planned_raises():
    with pytest.raises(ValidationError):
        VaccineDelivery(n_priming_doses=3, n_boost_doses=2, n_doses_received=8)

def test_adjuvant_none_with_detail_raises():
    with pytest.raises(ValidationError):
        VaccineDelivery(adjuvant="none", adjuvant_detail="poly-ICLC")

def test_rna_takes_adjuvant_none_and_lipoplex_in_formulation():
    vd = VaccineDelivery(adjuvant="none", formulation_detail="RNA-lipoplex, IV",
                         weeks_surgery_to_first_dose=6.0)
    assert vd.adjuvant == "none" and "lipoplex" in vd.formulation_detail

def test_bad_adjuvant_value_rejected():
    with pytest.raises(ValidationError):
        VaccineDelivery(adjuvant="lipid_nanoparticle")  # removed by Decision #1


# ---- cross-patient regimen-divergence nudge (unit, duck-typed) ----

def _vd(**kw):  # duck object with the regimen fields the helper reads
    base = dict(adjuvant=None, adjuvant_detail=None, formulation_detail=None, dose_amount_raw=None,
                dose_per_peptide_ug=None, dose_basis="unspecified", n_priming_doses=None,
                n_boost_doses=None, schedule_detail=None,
                weeks_surgery_to_first_dose=None, n_doses_received=None)
    base.update(kw)
    return NS(**base)

def _pt(platform, vd):
    return NS(vaccine_platform=platform, vaccine_delivery=vd)

def test_divergence_flagged_within_an_arm():
    rec = NS(patients=[_pt("synthetic_long_peptide", _vd(adjuvant="poly_iclc")),
                       _pt("synthetic_long_peptide", _vd(adjuvant="montanide"))])
    assert agent_core._regimen_divergence(rec) == ["synthetic_long_peptide"]

def test_identical_regimen_is_clean_even_if_perpatient_fields_differ():
    rec = NS(patients=[
        _pt("synthetic_long_peptide", _vd(adjuvant="poly_iclc", weeks_surgery_to_first_dose=10, n_doses_received=7)),
        _pt("synthetic_long_peptide", _vd(adjuvant="poly_iclc", weeks_surgery_to_first_dose=14, n_doses_received=3)),
    ])
    assert agent_core._regimen_divergence(rec) == []   # per-patient fields exempt

def test_different_platforms_are_not_compared():
    rec = NS(patients=[_pt("rna", _vd(adjuvant="none")),
                       _pt("dna", _vd(adjuvant="poly_iclc"))])
    assert agent_core._regimen_divergence(rec) == []

def test_single_patient_arm_cannot_diverge():
    rec = NS(patients=[_pt("rna", _vd(adjuvant="none")),
                       _pt("rna", _vd(adjuvant=None))])  # 2 same-arm but adjuvant None vs none differ
    assert agent_core._regimen_divergence(rec) == ["rna"]  # None != 'none' is a real divergence


# ---- the nudge is reachable + overridable through finalize_partial ----

def _finalize_two_patients(tmp_path, d1, d2, **flags):
    out = tmp_path / "r.json"
    agent_core.init_partial(str(out), json.dumps(META))
    agent_core.append_section(str(out), "patients",
                              json.dumps([_patient("P1", delivery=d1), _patient("P2", delivery=d2)]))
    return agent_core.finalize_partial(str(out), **flags)

def test_finalize_blocks_then_overrides_on_divergent_regimen(tmp_path):
    ok, msg = _finalize_two_patients(tmp_path, POLY, {**POLY, "adjuvant": "montanide"})
    assert not ok and "DIVERGENT" in msg and "allow_regimen_divergence=true" in msg
    ok2, _ = _finalize_two_patients(tmp_path, POLY, {**POLY, "adjuvant": "montanide"},
                                    allow_regimen_divergence=True)
    assert ok2

def test_finalize_clean_on_consistent_regimen(tmp_path):
    ok, msg = _finalize_two_patients(tmp_path, {**POLY, "weeks_surgery_to_first_dose": 10},
                                     {**POLY, "weeks_surgery_to_first_dose": 20})
    assert ok, msg   # same regimen, different latency -> no nudge


# ---- CohortLatency (v2.9.1): the paper-level home for unlabelled per-patient latency ----

def test_cohort_latency_roundtrips_with_median_and_benchmark():
    from schema import CohortLatency
    cl = CohortLatency(metric="surgery_to_first_vaccine", n_patients=16, median_value=9.4,
                       time_unit="weeks", benchmark_value=9.0, raw="median 9.4 wk, range 7.4-11")
    assert cl.metric == "surgery_to_first_vaccine" and cl.benchmark_value == 9.0

def test_cohort_latency_requires_median_or_raw():
    from schema import CohortLatency
    with pytest.raises(ValidationError):
        CohortLatency(metric="surgery_to_first_vaccine")        # neither median nor raw

def test_cohort_latency_raw_only_is_allowed():
    from schema import CohortLatency
    assert CohortLatency(metric="other", raw="time-to-vaccine list, no median given").raw

def test_cohort_latency_in_section_model_and_appendable(tmp_path):
    assert agent_core.SECTION_MODEL.get("cohort_latencies") == "CohortLatency"
    out = tmp_path / "r.json"
    agent_core.init_partial(str(out), json.dumps(META))
    ok, msg = agent_core.append_section(str(out), "cohort_latencies", json.dumps(
        [{"metric": "surgery_to_first_vaccine", "n_patients": 16, "median_value": 9.4,
          "benchmark_value": 9.0, "raw": "median 9.4 wk"}]))
    assert ok, msg
    part = json.loads((tmp_path / "r.json.partial.json").read_text())
    assert part["cohort_latencies"][0]["median_value"] == 9.4


# ---- back-compat ----

def test_li_reference_back_compat_no_delivery():
    # li was NOT refreshed -> proves a record with vaccine_delivery=None + cohort_latencies=[] validates.
    rec = ExtractedPaper(**json.loads((PKT / "reference_records" / "li_extracted.json").read_text()))
    assert all(p.vaccine_delivery is None for p in rec.patients)
    assert rec.cohort_latencies == []
    assert schema.SCHEMA_VERSION == "2.16.0"

def test_refreshed_references_carry_delivery_and_latency():
    # rojas + keskin refreshed 2026-06-07: per-patient vaccine_delivery + paper-level cohort_latencies,
    # with a NON-divergent (uniform) trial-constant regimen.
    import agent_core
    for n, adj, n_lat in (("rojas", "none", 2), ("keskin", "poly_iclc", 1)):
        rec = ExtractedPaper(**json.loads((PKT / "reference_records" / f"{n}_extracted.json").read_text()))
        assert all(p.vaccine_delivery is not None for p in rec.patients), n
        assert all(p.vaccine_delivery.adjuvant == adj for p in rec.patients), n
        assert len(rec.cohort_latencies) == n_lat, n
        assert agent_core._regimen_divergence(rec) == [], f"{n} regimen should be uniform"
