"""v2.17 TCR extraction MVP — value objects, entities, exactly-one-target discriminators,
functional-avidity placement (Option B), paper-level cross-references, tool-path wiring, and the
finalize QC gateway anchor (§5.7). Facts drawn from the two dry-run papers: 39972124 (PDAC) and
33479501 (Hu 2021 melanoma). Design: docs/TCRSEQ_EXTRACTION_DESIGN_2026-07-02_fable_rev.md."""
import json
import pathlib

import pytest
from pydantic import ValidationError

import agent_core
import schema
from schema import (
    TcrChain, TcrFunctionalAvidity, TcrSeqMethod, DataDeposition,
    NeoantigenTcrFlag, TcrClonotype, ExtractedEvidence, ExtractedPaper,
)

META = {"pmid": "39972124", "journal": "Nature", "year": 2025, "title": "TCR MVP test",
        "cohort_size": 1, "indication_summary": "PDAC"}
PROV = {"quoted_text": "q", "section_ref": "Methods"}


# ---- TcrChain (AIRR field names; lossless) ----

def test_chain_paired_cdr3_from_hu2021():
    # Hu 2021 lists full paired-ab CDR3 + V/J, e.g. TRBV9,TRBJ2-7,YASSVGAGTYEQYF
    c = TcrChain(locus="TRB", cdr3_aa="YASSVGAGTYEQYF", v_call="TRBV9", j_call="TRBJ2-7")
    assert c.cdr3_aa == "YASSVGAGTYEQYF" and c.v_call == "TRBV9"

def test_chain_junction_vs_cdr3_are_separate_fields():
    # the comparability crux: junction includes the conserved C (CVSS...), cdr3 excludes (YASS...)
    c = TcrChain(locus="TRB", junction_aa="CVSSPGTLAAGEQYF", cdr3_scheme="imgt_junction_aa")
    assert c.junction_aa.startswith("C") and c.cdr3_aa is None

def test_chain_lowercase_cdr3_normalized():
    assert TcrChain(cdr3_aa="yassvgagtyeqyf").cdr3_aa == "YASSVGAGTYEQYF"

def test_chain_raw_only_allowed():
    # 39972124: clones are opaque IDs, sequences only in the GEO deposit
    assert TcrChain(raw="TCR_1-2").raw == "TCR_1-2"

def test_chain_empty_rejected():
    with pytest.raises(ValidationError):
        TcrChain(locus="TRB")

def test_chain_bad_cdr3_alphabet_rejected():
    with pytest.raises(ValidationError):
        TcrChain(cdr3_aa="ZZZ1")


# ---- TcrFunctionalAvidity (concentration units; Option B) ----

def test_avidity_mut_wt_from_39972124():
    # P14 TCR2: NA 0.02 uM, WT 508.2 uM
    a = TcrFunctionalAvidity(ec50_mut_value=0.02, ec50_wt_value=508.2, unit="uM",
                             fold_difference=2.6e-5)
    assert a.ec50_mut_value == 0.02 and a.unit == "uM"

def test_avidity_non_reactive_sentinel_flagged():
    # 39972124 assigns EC50=1000 uM to non-reactive TCRs; non_reactive marks the censored value
    a = TcrFunctionalAvidity(ec50_wt_value=1000.0, ec50_mut_value=1.1, unit="uM", non_reactive=True)
    assert a.non_reactive is True

def test_avidity_parsed_value_needs_known_unit():
    with pytest.raises(ValidationError):
        TcrFunctionalAvidity(ec50_mut_value=0.7, unit="unknown")

def test_avidity_empty_rejected():
    with pytest.raises(ValidationError):
        TcrFunctionalAvidity(unit="uM")

def test_avidity_raw_only_allowed():
    assert TcrFunctionalAvidity(raw="EC50 sub-nanomolar").raw


# ---- TcrSeqMethod (one per distinct recipe) ----

def test_method_bulk_immunoseq_nt_v_j():
    # 39972124 bulk: immunoSEQ, TRB nt-CDR3 + V + J
    m = TcrSeqMethod(method_local_id="M_bulk", modality="bulk_tcr_seq",
                     platform="adaptive_immunoseq", software=["immunoseq_analyzer"],
                     clonotype_definition="trb_v_cdr3nt_j", pipeline_name="CloneTrack", **PROV)
    assert m.clonotype_definition == "trb_v_cdr3nt_j" and m.platform == "adaptive_immunoseq"

def test_method_single_cell_plate_based_with_patient_scope():
    # Hu 2021: plate-based single-cell for Pts 3-5 (a recipe covering specific patients)
    m = TcrSeqMethod(method_local_id="M_sc2", modality="single_cell_tcr_seq",
                     platform="plate_based_targeted", clonotype_definition="paired_ab_cdr3",
                     patient_scope=["P3", "P4", "P5"], **PROV)
    assert m.patient_scope == ["P3", "P4", "P5"]


# ---- DataDeposition (per-repo accession patterns) ----

@pytest.mark.parametrize("repo,acc", [
    ("geo", "GSE222011"),           # 39972124
    ("dbgap", "phs001451.v2.p1"),   # Hu 2021
    ("dbgap", "phs001451"),
    ("sra", "PRJNA123456"),
    ("immuneaccess", "some/token-url"),  # no strict pattern -> any non-empty
])
def test_deposition_valid_accessions(repo, acc):
    d = DataDeposition(repository=repo, accession=acc, access="controlled", **PROV)
    assert d.accession == acc

@pytest.mark.parametrize("repo,acc", [
    ("geo", "222011"),        # missing GSE prefix
    ("dbgap", "phs12"),       # too few digits
    ("geo", "GSE"),           # no number
])
def test_deposition_bad_accession_rejected(repo, acc):
    with pytest.raises(ValidationError):
        DataDeposition(repository=repo, accession=acc, **PROV)


# ---- NeoantigenTcrFlag (the flagship; three-way discriminator, no pool) ----

def test_flag_cd8_epitope():
    f = NeoantigenTcrFlag(target_kind="epitope", epitope_paper_id="E1",
                          patient_paper_id="P1", tcr_identified="cd8", **PROV)
    assert f.tcr_identified == "cd8"

def test_flag_pooled_patient_none_allowed():
    f = NeoantigenTcrFlag(target_kind="candidate", candidate_paper_id="C1",
                          tcr_identified="both", **PROV)
    assert f.patient_paper_id is None

def test_flag_two_targets_rejected():
    with pytest.raises(ValidationError):
        NeoantigenTcrFlag(target_kind="epitope", epitope_paper_id="E1",
                          immunizing_peptide_paper_id="I1", tcr_identified="cd4", **PROV)

def test_flag_missing_own_target_rejected():
    with pytest.raises(ValidationError):
        NeoantigenTcrFlag(target_kind="epitope", tcr_identified="cd4", **PROV)

def test_flag_has_no_pool_target():
    with pytest.raises(ValidationError):
        NeoantigenTcrFlag(target_kind="pool", tcr_identified="cd4", **PROV)  # 'pool' not in TcrTarget


# ---- TcrClonotype (specificity discriminator + strength discipline) ----

def test_clonotype_minimal_shell():
    c = TcrClonotype(paper_local_id="TCR_1-2", patient_paper_id="P1",
                     receptor_modality="single_cell_tcr_seq", **PROV)
    assert c.chains == [] and c.specificity_evidence == "not_established"

def test_clonotype_functional_clone_with_specificity():
    c = TcrClonotype(paper_local_id="C1", specific_for_kind="epitope", epitope_paper_id="E1",
                     specificity_evidence="functional_tcr_clone",
                     validation_evidence_ref="EV_hero1",
                     restricting_hla="hla-a*02:01", **PROV)
    assert c.epitope_paper_id == "E1" and c.restricting_hla == "HLA-A*02:01"

def test_clonotype_inferred_association_forces_review():
    with pytest.raises(ValidationError):
        TcrClonotype(paper_local_id="C1", specific_for_kind="epitope", epitope_paper_id="E1",
                     specificity_evidence="inferred_association", needs_review=False, **PROV)

def test_clonotype_inferred_association_ok_when_flagged():
    c = TcrClonotype(paper_local_id="C1", specific_for_kind="epitope", epitope_paper_id="E1",
                     specificity_evidence="inferred_association", needs_review=True, **PROV)
    assert c.needs_review is True

def test_clonotype_strength_without_target_rejected():
    with pytest.raises(ValidationError):
        TcrClonotype(paper_local_id="C1", specificity_evidence="functional_tcr_clone", **PROV)

def test_clonotype_target_id_without_kind_rejected():
    with pytest.raises(ValidationError):
        TcrClonotype(paper_local_id="C1", epitope_paper_id="E1", **PROV)

def test_clonotype_two_specificity_targets_rejected():
    with pytest.raises(ValidationError):
        TcrClonotype(paper_local_id="C1", specific_for_kind="epitope", epitope_paper_id="E1",
                     candidate_paper_id="X1", specificity_evidence="multimer_sort", **PROV)


# ---- TcrObservation (v2.18): longitudinal / per-tissue trajectory ----

from schema import TcrObservation  # noqa: E402

def test_observation_frequency_as_reported():
    # Hu tet-spec track: fraction 0.0003 at Week 12 in blood
    o = TcrObservation(tissue="blood", sorted_compartment="whole_pbmc", timepoint_phase="on_treatment",
                       timepoint_label="Week 12", frequency_value=0.0003, frequency_basis="fraction")
    assert o.frequency_value == 0.0003 and o.frequency_basis == "fraction"

def test_observation_percent_not_coerced():
    # a percent stays a percent (basis-tagged), never silently converted to a fraction
    o = TcrObservation(tissue="tumor_relapsed", frequency_value=8.5, frequency_basis="percent")
    assert o.frequency_value == 8.5 and o.frequency_basis == "percent"

def test_observation_reuses_evidence_vocab():
    o = TcrObservation(assay_stimulation="in_vitro_expanded", timepoint_phase="memory", raw="Long-term")
    assert o.assay_stimulation == "in_vitro_expanded" and o.timepoint_phase == "memory"

def test_observation_empty_rejected():
    with pytest.raises(ValidationError):
        TcrObservation(tissue="blood", timepoint_label="Week 4")  # no frequency/count/raw

def test_clonotype_carries_observations():
    c = TcrClonotype(paper_local_id="TCR_track1", chain_pairing="trb_only",
                     chains=[{"locus": "TRB", "junction_aa": "CASSRGLAGYNEQFF"}],
                     observations_identity="paper_explicit",
                     observations=[{"tissue": "blood", "timepoint_label": "Pre-Vax", "frequency_value": 0.0,
                                    "frequency_basis": "fraction", "raw": "pre-vax"},
                                   {"tissue": "blood", "timepoint_label": "Week 16", "frequency_value": 0.085,
                                    "frequency_basis": "fraction"}], **PROV)
    assert len(c.observations) == 2 and c.observations_identity == "paper_explicit"

def test_single_observation_not_asserted_ok():
    # one observation asserts no grouping -> not_asserted is fine, no needs_review needed
    c = TcrClonotype(paper_local_id="C1",
                     observations=[{"tissue": "blood", "frequency_value": 0.1, "frequency_basis": "fraction"}],
                     **PROV)
    assert c.observations_identity == "not_asserted" and c.needs_review is False

def test_multi_observation_not_asserted_forces_review():
    # decision 5.1: grouping >1 observation the paper didn't tie together must set needs_review
    with pytest.raises(ValidationError):
        TcrClonotype(paper_local_id="C1", observations_identity="not_asserted",
                     observations=[{"tissue": "blood", "frequency_value": 0.1, "frequency_basis": "fraction"},
                                   {"tissue": "tumor", "frequency_value": 0.2, "frequency_basis": "fraction"}],
                     needs_review=False, **PROV)

def test_multi_observation_not_asserted_ok_when_flagged():
    c = TcrClonotype(paper_local_id="C1", observations_identity="not_asserted", needs_review=True,
                     observations=[{"tissue": "blood", "frequency_value": 0.1, "frequency_basis": "fraction"},
                                   {"tissue": "tumor", "frequency_value": 0.2, "frequency_basis": "fraction"}],
                     **PROV)
    assert c.needs_review is True

def test_multi_observation_paper_explicit_no_review_needed():
    c = TcrClonotype(paper_local_id="C1", observations_identity="paper_explicit",
                     observations=[{"tissue": "blood", "frequency_value": 0.1, "frequency_basis": "fraction"},
                                   {"tissue": "tumor", "frequency_value": 0.2, "frequency_basis": "fraction"}],
                     **PROV)
    assert c.needs_review is False and len(c.observations) == 2

def test_clonotype_without_observations_backcompat():
    c = TcrClonotype(paper_local_id="C1", **PROV)
    assert c.observations == [] and c.observations_identity == "not_asserted"


# ---- ExtractedEvidence: functional_avidity + evidence_local_id (Option B) ----

def test_evidence_carries_functional_avidity():
    ev = ExtractedEvidence(patient_paper_id="P1", target_kind="epitope", epitope_paper_id="E1",
                           assay="tcr_reporter", outcome="immunogenic", evidence_local_id="EV_hero1",
                           functional_avidity={"ec50_mut_value": 0.7, "unit": "uM"}, **PROV)
    assert ev.functional_avidity.ec50_mut_value == 0.7 and ev.evidence_local_id == "EV_hero1"


# ---- ExtractedPaper cross-references ----

def _paper(**kw):
    base = dict(
        patients=[{"paper_local_id": "P1", "indication": "PDAC", "n_peptides_synthesized": 1,
                   "n_peptides_immunogenic": 1, **PROV}],
        epitopes=[{"paper_local_id": "E1", "sequence": "SLLQHLIGL", "is_neoantigen": True,
                   "mhc_class": "I", "predicted_affinity": {"unit": "unknown", "raw": "n/a"}, **PROV}],
    )
    base.update(kw)
    return ExtractedPaper(**META, **base)

def test_paper_with_full_tcr_bundle_validates():
    rec = _paper(
        tcr_seq_status="both",
        data_depositions=[{"repository": "geo", "accession": "GSE222011", **PROV}],
        tcr_seq_methods=[{"method_local_id": "M1", "platform": "adaptive_immunoseq", **PROV}],
        neoantigen_tcr_flags=[{"target_kind": "epitope", "epitope_paper_id": "E1",
                               "patient_paper_id": "P1", "tcr_identified": "cd8", **PROV}],
        evidence=[{"patient_paper_id": "P1", "target_kind": "epitope", "epitope_paper_id": "E1",
                   "assay": "tcr_reporter", "outcome": "immunogenic", "evidence_local_id": "EV1",
                   "functional_avidity": {"ec50_mut_value": 0.7, "unit": "uM"}, **PROV}],
        tcr_clonotypes=[{"paper_local_id": "TCR1", "patient_paper_id": "P1",
                         "specific_for_kind": "epitope", "epitope_paper_id": "E1",
                         "specificity_evidence": "functional_tcr_clone", "method_ref": "M1",
                         "validation_evidence_ref": "EV1", **PROV}],
    )
    assert rec.tcr_seq_status == "both" and rec.tcr_clonotypes[0].validation_evidence_ref == "EV1"

def test_paper_pooled_clonotype_patient_none_ok():
    rec = _paper(tcr_clonotypes=[{"paper_local_id": "TCR1", **PROV}])  # patient None = cohort-level
    assert rec.tcr_clonotypes[0].patient_paper_id is None

def test_paper_flag_unknown_target_rejected():
    with pytest.raises(ValidationError):
        _paper(neoantigen_tcr_flags=[{"target_kind": "epitope", "epitope_paper_id": "NOPE",
                                      "tcr_identified": "cd8", **PROV}])

def test_paper_clonotype_unknown_patient_rejected():
    with pytest.raises(ValidationError):
        _paper(tcr_clonotypes=[{"paper_local_id": "TCR1", "patient_paper_id": "GHOST", **PROV}])

def test_paper_clonotype_dangling_method_ref_rejected():
    with pytest.raises(ValidationError):
        _paper(tcr_clonotypes=[{"paper_local_id": "TCR1", "method_ref": "M_missing", **PROV}])

def test_paper_validation_ref_must_be_tcr_reporter():
    # points at an evidence row that exists but is NOT a tcr_reporter row -> reject
    with pytest.raises(ValidationError):
        _paper(
            evidence=[{"patient_paper_id": "P1", "target_kind": "epitope", "epitope_paper_id": "E1",
                       "assay": "elispot", "outcome": "immunogenic", "evidence_local_id": "EV1", **PROV}],
            tcr_clonotypes=[{"paper_local_id": "TCR1", "validation_evidence_ref": "EV1", **PROV}],
        )

def test_paper_duplicate_method_id_rejected():
    with pytest.raises(ValidationError):
        _paper(tcr_seq_methods=[{"method_local_id": "M1", **PROV},
                                {"method_local_id": "M1", **PROV}])

def test_paper_duplicate_clonotype_id_rejected():
    with pytest.raises(ValidationError):
        _paper(tcr_clonotypes=[{"paper_local_id": "TCR1", **PROV},
                               {"paper_local_id": "TCR1", **PROV}])


# ---- tool-path wiring (P19-lesson regression: an entity unreachable via SECTION_MODEL is dead) ----

@pytest.mark.parametrize("section,model", [
    ("data_depositions", "DataDeposition"),
    ("tcr_seq_methods", "TcrSeqMethod"),
    ("neoantigen_tcr_flags", "NeoantigenTcrFlag"),
    ("tcr_clonotypes", "TcrClonotype"),
])
def test_tcr_sections_in_section_model(section, model):
    assert agent_core.SECTION_MODEL.get(section) == model

def test_deposition_appendable_through_tool_path(tmp_path):
    out = tmp_path / "r.json"
    agent_core.init_partial(str(out), json.dumps(META))
    ok, msg = agent_core.append_section(str(out), "data_depositions", json.dumps([
        {"repository": "dbgap", "accession": "phs001451.v2.p1", "data_type": "wes",
         "access": "controlled", "quoted_text": "dbGaP phs001451.v2.p1", "section_ref": "Data availability"}]))
    assert ok, msg
    part = json.loads((tmp_path / "r.json.partial.json").read_text())
    assert part["data_depositions"][0]["accession"] == "phs001451.v2.p1"


# ---- QC gateway anchor (§5.7): tcr_seq_status asserts TCR-seq but nothing extracted ----

def test_gateway_gap_silent_when_no_tcr_claim():
    assert agent_core._tcr_gateway_gap(_paper()) == ""                    # status None
    assert agent_core._tcr_gateway_gap(_paper(tcr_seq_status="none")) == ""

def test_gateway_gap_fires_when_claimed_but_empty():
    assert agent_core._tcr_gateway_gap(_paper(tcr_seq_status="both"))
    assert agent_core._tcr_gateway_gap(_paper(tcr_seq_status="mentioned_only"))

def test_gateway_gap_satisfied_by_deposition():
    rec = _paper(tcr_seq_status="both",
                 data_depositions=[{"repository": "geo", "accession": "GSE222011", **PROV}])
    assert agent_core._tcr_gateway_gap(rec) == ""

def test_gateway_gap_satisfied_by_method():
    rec = _paper(tcr_seq_status="both",
                 tcr_seq_methods=[{"method_local_id": "M1", "platform": "adaptive_immunoseq", **PROV}])
    assert agent_core._tcr_gateway_gap(rec) == ""

def test_gateway_override_is_hard_routes_to_needs_review():
    # a declared-but-unfilled Tier-A gap is a real recall gap: HARD override -> needs_review.
    assert "allow_missing_tcr_gateway" in agent_core.HARD_OVERRIDES
    assert agent_core.overrides_are_soft_only(["allow_missing_tcr_gateway"]) is False


def _init_min_paper(tmp_path, **meta_extra):
    # 0 declared peptides so the peptide-count gate is satisfied with no immunizing_peptides —
    # isolates the TCR gateway anchor as the only thing that can fire.
    out = tmp_path / "r.json"
    agent_core.init_partial(str(out), json.dumps({**META, **meta_extra}))
    ok, msg = agent_core.append_section(str(out), "patients", json.dumps([
        {"paper_local_id": "P1", "indication": "PDAC", "n_peptides_synthesized": 0,
         "n_peptides_immunogenic": 0, **PROV}]))
    assert ok, msg
    return out

def test_finalize_blocks_on_tcr_gateway(tmp_path):
    out = _init_min_paper(tmp_path, tcr_seq_status="both")
    ok, msg = agent_core.finalize_partial(str(out))
    assert not ok and "TCR gateway" in msg

def test_finalize_override_records_audit_and_routes_review(tmp_path):
    out = _init_min_paper(tmp_path, tcr_seq_status="mentioned_only")
    ok, msg = agent_core.finalize_partial(str(out), allow_missing_tcr_gateway=True)
    assert ok, msg
    rec = json.loads(pathlib.Path(str(out)).read_text())
    assert "allow_missing_tcr_gateway" in rec.get("finalize_overrides_used", [])
    assert agent_core.overrides_are_soft_only(rec["finalize_overrides_used"]) is False

def test_finalize_clean_when_gateway_satisfied(tmp_path):
    out = _init_min_paper(tmp_path, tcr_seq_status="both")
    agent_core.append_section(str(out), "data_depositions", json.dumps([
        {"repository": "geo", "accession": "GSE222011", "data_type": "paired_scvdj",
         "access": "open", **PROV}]))
    ok, msg = agent_core.finalize_partial(str(out))
    assert ok, msg
    rec = json.loads(pathlib.Path(str(out)).read_text())
    assert rec.get("finalize_overrides_used", []) == []


# ---- status-consistency nudge (§5.3): TCR content extracted but tcr_seq_status wrong/missing ----

def test_status_consistency_silent_without_tcr_content():
    assert agent_core._tcr_status_consistency_gap(_paper(tcr_seq_status=None)) == ""

def test_status_consistency_fires_when_content_but_status_none():
    # the Phase-4 Rojas miss: methods present, status left None
    rec = _paper(tcr_seq_methods=[{"method_local_id": "M1", "modality": "single_cell_tcr_seq", **PROV}])
    assert agent_core._tcr_status_consistency_gap(rec)

def test_status_consistency_fires_on_understated_modality():
    # the Phase-4 Hu miss: status=single_cell but a bulk method is also present -> should be 'both'
    rec = _paper(tcr_seq_status="single_cell", tcr_seq_methods=[
        {"method_local_id": "M1", "modality": "single_cell_tcr_seq", **PROV},
        {"method_local_id": "M2", "modality": "bulk_tcr_seq", **PROV}])
    assert "expected 'both'" in agent_core._tcr_status_consistency_gap(rec)

def test_status_consistency_silent_when_matching():
    rec = _paper(tcr_seq_status="both", tcr_seq_methods=[
        {"method_local_id": "M1", "modality": "single_cell_tcr_seq", **PROV},
        {"method_local_id": "M2", "modality": "bulk_tcr_seq", **PROV}])
    assert agent_core._tcr_status_consistency_gap(rec) == ""

def test_status_consistency_mentioned_only_contradicts_content():
    rec = _paper(tcr_seq_status="mentioned_only",
                 tcr_seq_methods=[{"method_local_id": "M1", "modality": "bulk_tcr_seq", **PROV}])
    assert agent_core._tcr_status_consistency_gap(rec)

def test_status_mismatch_override_is_hard():
    assert "allow_tcr_status_mismatch" in agent_core.HARD_OVERRIDES
    assert agent_core.overrides_are_soft_only(["allow_tcr_status_mismatch"]) is False

def test_finalize_blocks_on_status_mismatch(tmp_path):
    out = _init_min_paper(tmp_path)  # tcr_seq_status unset
    agent_core.append_section(str(out), "tcr_seq_methods", json.dumps([
        {"method_local_id": "M1", "modality": "bulk_tcr_seq", "platform": "adaptive_immunoseq", **PROV}]))
    ok, msg = agent_core.finalize_partial(str(out))
    assert not ok and "TCR status" in msg


# ---- large-xlsx streaming (Task 2 fix): big workbooks read read-only so CDR3 sheets are reachable ----

def test_read_table_streams_large_files(tmp_path, monkeypatch):
    openpyxl = pytest.importorskip("openpyxl")
    p = tmp_path / "clonotypes.xlsx"
    wb = openpyxl.Workbook(); ws = wb.active; ws.title = "Suppl 8 CDR3b"
    ws.append(["CDR3b", "V", "J"]); ws.append(["CASSLGNQPQHF", "TRBV9", "TRBJ2-7"]); wb.save(str(p))
    # force the >threshold streaming branch on this tiny file
    monkeypatch.setattr(agent_core, "_BIG_XLSX_BYTES", 10)
    out = agent_core.read_table_rows(str(p))
    assert "CASSLGNQPQHF" in out and "TRBV9" in out          # rows read via read_only stream
    assert "Suppl 8 CDR3b" in out                            # sheet name surfaced

def test_read_table_normal_path_unchanged(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    p = tmp_path / "small.xlsx"
    wb = openpyxl.Workbook(); ws = wb.active; ws.append(["a", "b"]); ws.append([1, 2]); wb.save(str(p))
    out = agent_core.read_table_rows(str(p))   # small file -> normal (non-read-only) path
    assert "1" in out and "2" in out
