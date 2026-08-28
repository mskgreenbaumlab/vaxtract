"""CD4 evidence must target the long immunizing peptide, not a nested class-I 8-11mer.

Keskin 30568305 (archived: outputs/archive/2026-08-22_classii_keskin/classii_keskin_2026-08-22) typed t_cell_subset=cd4 honestly then hung those
rows on invented EPI_RESP_* aliases of Table-S5 9-mers. Gold already points CD4 at IMPs.
"""
import json
import pathlib
from types import SimpleNamespace as NS

import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]
META = {"pmid": "12345678", "journal": "Test J", "year": 2024, "title": "CD4 grain test",
        "cohort_size": 1, "indication_summary": "glioblastoma"}

IMP_SEQ = "SEMEENLANRSHAKLETALRDSSRVL"  # Keskin GPC1 26-mer
MER9 = "HAKLETALR"
AFF = {"unit": "unknown", "raw": "n/a"}


def _patient():
    return {"quoted_text": "pt P1", "section_ref": "M", "paper_local_id": "P1",
            "indication": "glioblastoma", "n_peptides_synthesized": 1, "n_peptides_immunogenic": 1}


def _imp():
    return {"paper_local_id": "IMP_GPC1", "sequence": IMP_SEQ, "is_neoantigen": True,
            "gene_symbol": "GPC1", "quoted_text": "immunizing peptide GPC1", "section_ref": "Table S5"}


def _class_i_epitope():
    return {"paper_local_id": "EPI_HAKLETALR", "sequence": MER9, "is_neoantigen": True,
            "mhc_class": "I", "hla_allele": "HLA-A*68:01", "gene_symbol": "GPC1",
            "parent_peptide_ids": ["IMP_GPC1"],
            "predicted_affinity": AFF,
            "quoted_text": "predicted MHC class I neoepitope HAKLETALR", "section_ref": "Table S5"}


def _class_ii_epitope():
    return {"paper_local_id": "EPI_DR", "sequence": "SEMEENLANRSHAKL", "is_neoantigen": True,
            "mhc_class": "II", "hla_allele": "HLA-DRB1*04:01", "gene_symbol": "GPC1",
            "quoted_text": "HLA-DRB1*04:01-restricted class II epitope", "section_ref": "Results"}


def _ev_cd4_imp():
    return {"patient_paper_id": "P1", "target_kind": "immunizing_peptide",
            "immunizing_peptide_paper_id": "IMP_GPC1", "assay": "elispot", "outcome": "immunogenic",
            "quoted_text": "CD4+ T cell responses against mutated GPC1 neoepitope",
            "section_ref": "Figure 2b", "t_cell_subset": "cd4",
            "magnitude": {"raw": "positive"}}


def _ev_cd4_9mer():
    return {"patient_paper_id": "P1", "target_kind": "epitope",
            "epitope_paper_id": "EPI_HAKLETALR", "assay": "elispot", "outcome": "immunogenic",
            "quoted_text": "CD4+ T cell responses against mutated GPC1 neoepitope",
            "section_ref": "Figure 2b", "t_cell_subset": "cd4",
            "magnitude": {"raw": "positive"}}


def _ev_cd8_9mer():
    return {"patient_paper_id": "P1", "target_kind": "epitope",
            "epitope_paper_id": "EPI_HAKLETALR", "assay": "elispot", "outcome": "immunogenic",
            "quoted_text": "CD8+ T cell responses to predicted class I epitopes",
            "section_ref": "Results", "t_cell_subset": "cd8", "mhc_class": "class_i",
            "magnitude": {"raw": "positive"}}


def _ev_cd4_class_ii_epi():
    return {"patient_paper_id": "P1", "target_kind": "epitope",
            "epitope_paper_id": "EPI_DR", "assay": "elispot", "outcome": "immunogenic",
            "quoted_text": "CD4+ HLA-DRB1*04:01-restricted response",
            "section_ref": "Results", "t_cell_subset": "cd4", "mhc_class": "class_ii",
            "magnitude": {"raw": "positive"}}


def _assemble(out, *, epitopes, evidence):
    agent_core.init_partial(str(out), json.dumps(META))
    agent_core.append_section(str(out), "patients", json.dumps([_patient()]))
    agent_core.append_section(str(out), "immunizing_peptides", json.dumps([_imp()]))
    if epitopes:
        agent_core.append_section(str(out), "epitopes", json.dumps(epitopes))
    agent_core.append_section(str(out), "evidence", json.dumps(evidence))


def test_gap_empty_for_cd4_on_imp():
    rec = NS(
        epitopes=[NS(paper_local_id="EPI_HAKLETALR", sequence=MER9, mhc_class="I",
                     parent_peptide_ids=["IMP_GPC1"])],
        evidence=[NS(t_cell_subset="cd4", target_kind="immunizing_peptide",
                     epitope_paper_id=None, immunizing_peptide_paper_id="IMP_GPC1")],
    )
    assert agent_core._cd4_class_i_minimal_target_gap(rec) == []


def test_gap_fires_for_cd4_on_class_i_9mer():
    rec = NS(
        epitopes=[NS(paper_local_id="EPI_HAKLETALR", sequence=MER9, mhc_class="I",
                     parent_peptide_ids=["IMP_GPC1"])],
        evidence=[NS(t_cell_subset="cd4", target_kind="epitope",
                     epitope_paper_id="EPI_HAKLETALR")],
    )
    msgs = agent_core._cd4_class_i_minimal_target_gap(rec)
    assert msgs and "EPI_HAKLETALR" in msgs[0] and "IMP_GPC1" in msgs[0]


def test_gap_empty_for_cd8_on_class_i_9mer():
    rec = NS(
        epitopes=[NS(paper_local_id="EPI_HAKLETALR", sequence=MER9, mhc_class="I",
                     parent_peptide_ids=["IMP_GPC1"])],
        evidence=[NS(t_cell_subset="cd8", target_kind="epitope",
                     epitope_paper_id="EPI_HAKLETALR")],
    )
    assert agent_core._cd4_class_i_minimal_target_gap(rec) == []


def test_gap_empty_for_cd4_on_class_ii_epitope():
    rec = NS(
        epitopes=[NS(paper_local_id="EPI_DR", sequence="SEMEENLANRSHAKL", mhc_class="II",
                     parent_peptide_ids=[])],
        evidence=[NS(t_cell_subset="cd4", target_kind="epitope", epitope_paper_id="EPI_DR")],
    )
    assert agent_core._cd4_class_i_minimal_target_gap(rec) == []


def test_finalize_cd4_on_imp_ok(tmp_path):
    out = tmp_path / "r.json"
    _assemble(out, epitopes=[_class_i_epitope()], evidence=[_ev_cd4_imp()])
    ok, msg = agent_core.finalize_partial(str(out))
    assert ok, msg
    rec = json.loads(out.read_text())
    assert rec["evidence"][0]["target_kind"] == "immunizing_peptide"
    assert "allow_cd4_on_class_i_minimal" not in (rec.get("finalize_overrides_used") or [])


def test_finalize_blocks_cd4_on_class_i_9mer(tmp_path):
    out = tmp_path / "r.json"
    _assemble(out, epitopes=[_class_i_epitope()], evidence=[_ev_cd4_9mer()])
    ok, msg = agent_core.finalize_partial(str(out))
    assert not ok and "CD4 evidence grain" in msg and "allow_cd4_on_class_i_minimal=true" in msg
    assert "IMP_GPC1" in msg


def test_finalize_override_records_cd4_grain(tmp_path):
    out = tmp_path / "r.json"
    _assemble(out, epitopes=[_class_i_epitope()], evidence=[_ev_cd4_9mer()])
    ok, msg = agent_core.finalize_partial(str(out), allow_cd4_on_class_i_minimal=True)
    assert ok, msg
    rec = json.loads(out.read_text())
    assert rec["finalize_overrides_used"] == ["allow_cd4_on_class_i_minimal"]


def test_finalize_cd8_on_class_i_9mer_ok(tmp_path):
    out = tmp_path / "r.json"
    _assemble(out, epitopes=[_class_i_epitope()], evidence=[_ev_cd8_9mer()])
    ok, msg = agent_core.finalize_partial(str(out))
    assert ok, msg


def test_finalize_cd4_on_class_ii_epitope_ok(tmp_path):
    out = tmp_path / "r.json"
    _assemble(out, epitopes=[_class_ii_epitope()], evidence=[_ev_cd4_class_ii_epi()])
    ok, msg = agent_core.finalize_partial(str(out))
    assert ok, msg


def test_tool_schema_exposes_cd4_grain_override():
    from vaxtract.tool_registry import finalize_override_names, get_tool
    spec = get_tool("finalize")
    assert "allow_cd4_on_class_i_minimal" in spec.input_schema["properties"]
    assert "allow_cd4_on_class_i_minimal" in finalize_override_names()


def test_override_is_hard():
    assert "allow_cd4_on_class_i_minimal" in agent_core.HARD_OVERRIDES
    assert not agent_core.overrides_are_soft_only(["allow_cd4_on_class_i_minimal"])
