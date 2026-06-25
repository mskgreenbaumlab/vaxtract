"""v2.11 P20 — clonality + gene-level antigen dynamics. Piece A (NeoantigenMutation entity, incl. the
tool-path wiring test that the P19 funnel taught us to add) and Piece B (clonality on _PeptideCore),
plus back-compat with the refreshed references."""
import json
import pathlib

import pytest
from pydantic import ValidationError

import agent_core
import schema
from schema import NeoantigenMutation, VafPoint, ImmunizingPeptide, ExtractedPaper

PKT = pathlib.Path(__file__).resolve().parents[1]
META = {"pmid": "12345678", "journal": "Test J", "year": 2024, "title": "P20 test",
        "cohort_size": 1, "indication_summary": "melanoma"}


# ---- VafPoint ----

def test_vafpoint_roundtrips():
    v = VafPoint(timepoint_label="recurrent", timepoint_phase="post_vaccine", value=0.166)
    assert v.value == 0.166 and v.timepoint_label == "recurrent"

def test_vafpoint_raw_only_allowed():
    assert VafPoint(raw="0.17 (range 0.1-0.2)").raw

def test_vafpoint_empty_rejected():
    with pytest.raises(ValidationError):
        VafPoint()

def test_vafpoint_out_of_range_rejected():
    with pytest.raises(ValidationError):
        VafPoint(value=1.5)


# ---- NeoantigenMutation (Piece A) — Table S6 shape ----

def _mut(**kw):
    base = dict(paper_local_id="MUT1", patient_paper_id="P1", gene_symbol="RLF",
                genomic_change="chr1:40705080 A>T", quoted_text="RLF row", section_ref="Table S6")
    base.update(kw)
    return NeoantigenMutation(**base)

def test_table_s6_emerged_mutation():
    m = _mut(status="emerged", clonality="subclonal", cancer_cell_fraction=0.5,
             hla_restrictions=["HLA-A*11:01", "HLA-DRB1*08:03"],
             vaf=[VafPoint(timepoint_label="primary", value=0.0),
                  VafPoint(timepoint_label="recurrent", value=0.166)])
    assert m.status == "emerged" and m.clonality == "subclonal"
    assert [v.value for v in m.vaf] == [0.0, 0.166]
    assert m.hla_restrictions == ["HLA-A*11:01", "HLA-DRB1*08:03"]

def test_lost_status_for_swanton_case():
    assert _mut(status="lost").status == "lost"

def test_mutation_needs_gene_or_genomic():
    with pytest.raises(ValidationError):
        NeoantigenMutation(paper_local_id="M", patient_paper_id="P1",
                           quoted_text="q", section_ref="s")  # neither gene nor genomic_change

def test_genomic_only_mutation_ok():  # gene unknown but coords given
    assert NeoantigenMutation(paper_local_id="M", patient_paper_id="P1",
                              genomic_change="chr2:29016794 T>G", quoted_text="q", section_ref="s")

def test_off_vocab_clonality_and_status_rejected():
    with pytest.raises(ValidationError):
        _mut(clonality="branchy")
    with pytest.raises(ValidationError):
        _mut(status="vanished")

def test_hla_restriction_pair_canonicalizes():
    m = _mut(hla_restrictions=["hla-dpa1*01:03/dpb1*02:01"])
    assert m.hla_restrictions[0] == "HLA-DPA1*01:03/DPB1*02:01"


# ---- Piece B: clonality on _PeptideCore ----

def test_peptide_carries_clonality():
    ip = ImmunizingPeptide(paper_local_id="i1", sequence="SLLQHLIGL", is_neoantigen=True,
                           clonality="clonal", cancer_cell_fraction=0.98, wgd_timing="pre_wgd",
                           quoted_text="q", section_ref="s")
    assert ip.clonality == "clonal" and ip.cancer_cell_fraction == 0.98 and ip.wgd_timing == "pre_wgd"

def test_peptide_clonality_defaults_none():
    ip = ImmunizingPeptide(paper_local_id="i1", sequence="SLLQHLIGL", is_neoantigen=True,
                           quoted_text="q", section_ref="s")
    assert ip.clonality is None and ip.cancer_cell_fraction is None and ip.wgd_timing is None

def test_peptide_ccf_out_of_range_rejected():
    with pytest.raises(ValidationError):
        ImmunizingPeptide(paper_local_id="i1", sequence="SLLQHLIGL", is_neoantigen=True,
                          cancer_cell_fraction=2.0, quoted_text="q", section_ref="s")


# ---- wiring: SECTION_MODEL + reachable through the tool path (P19-lesson regression) ----

def test_mutation_in_section_model():
    assert agent_core.SECTION_MODEL.get("neoantigen_mutations") == "NeoantigenMutation"

def test_mutation_appendable_through_tool_path(tmp_path):
    out = tmp_path / "r.json"
    agent_core.init_partial(str(out), json.dumps(META))
    ok, msg = agent_core.append_section(str(out), "neoantigen_mutations", json.dumps([
        {"paper_local_id": "MUT1", "patient_paper_id": "P1", "gene_symbol": "RLF",
         "genomic_change": "chr1:40705080 A>T", "status": "emerged", "clonality": "subclonal",
         "hla_restrictions": ["HLA-A*11:01"],
         "vaf": [{"timepoint_label": "primary", "value": 0.0},
                 {"timepoint_label": "recurrent", "value": 0.166}],
         "quoted_text": "RLF", "section_ref": "Table S6"}]))
    assert ok, msg
    part = json.loads((tmp_path / "r.json.partial.json").read_text())
    assert part["neoantigen_mutations"][0]["status"] == "emerged"


def test_peptide_grain_clonality_persists_through_tool_path(tmp_path):
    """Piece B regression: a peptide-grain CCF/Clonal sheet (39910301 In Vitro) must land
    clonality + CCF ON THE PEPTIDE through the tool path (the 2026-06-08 re-run dropped the
    Clonal 0/1 column). Mirrors row 16097-101-3: CCF=0.18, Clonal=0 -> subclonal."""
    out = tmp_path / "r.json"
    agent_core.init_partial(str(out), json.dumps(META))
    ok, msg = agent_core.append_section(str(out), "immunizing_peptides", json.dumps([
        {"paper_local_id": "16097-101-3", "sequence": "STRDPLSEITKQEKDFLWSHRHY",
         "gene_symbol": "PIK3CA", "is_neoantigen": True, "mutation": "PIK3CA:E545K",
         "clonality": "subclonal", "cancer_cell_fraction": 0.18,
         "quoted_text": "PIK3CA|p.E545K CCF 0.18 Clonal 0", "section_ref": "Supp. In Vitro"}]))
    assert ok, msg
    part = json.loads((tmp_path / "r.json.partial.json").read_text())
    pep = part["immunizing_peptides"][0]
    assert pep["clonality"] == "subclonal" and pep["cancer_cell_fraction"] == 0.18


# ---- v2.11.2 (P20.2): Clonal 0/1 coercion (the add_table DSL has no value transform) ----

@pytest.mark.parametrize("raw,expected", [
    (1, "clonal"), (0, "subclonal"), (1.0, "clonal"), (0.0, "subclonal"),
    ("1", "clonal"), ("0", "subclonal"), (" 1 ", "clonal"),
    (True, "clonal"), (False, "subclonal"),
    ("clonal", "clonal"), ("subclonal", "subclonal"), ("unknown", "unknown"), (None, None),
])
def test_clonality_coercion_on_peptide(raw, expected):
    ip = ImmunizingPeptide(paper_local_id="P1", sequence="SIINFEKL", is_neoantigen=True,
                           clonality=raw, quoted_text="q", section_ref="s")
    assert ip.clonality == expected

def test_clonality_coercion_on_mutation():
    m = NeoantigenMutation(paper_local_id="M1", patient_paper_id="P1", gene_symbol="RLF",
                           clonality=0, quoted_text="q", section_ref="s")
    assert m.clonality == "subclonal"

@pytest.mark.parametrize("bad", [2, 7, "2", "yes", "clonalish"])
def test_clonality_coercion_leaves_offvocab_to_fail(bad):
    # only an exact 0/1 flag is coerced; anything else still hits the Literal and errors
    with pytest.raises(ValidationError):
        ImmunizingPeptide(paper_local_id="P1", sequence="SIINFEKL", is_neoantigen=True,
                          clonality=bad, quoted_text="q", section_ref="s")

def test_clonal_column_mapped_directly_via_add_table_dsl():
    """End-to-end of the real fix: add_table maps {'col':'Clonal'} (raw 0/1) and the schema
    coerces -- this is what 39910301's In Vitro sheet needs and what 3 live runs failed to do."""
    import table_map
    row = {"Peptide_ID": "16097-101-10", "Vaccine_Peptide": "KNQEVTIKALKEKIREYE",
           "Hugo_Symbol": "CUX1", "CCF": 0.95, "Clonal": 1}
    mapped = table_map.apply_mapping(row, {
        "paper_local_id": {"col": "Peptide_ID"}, "sequence": {"col": "Vaccine_Peptide"},
        "gene_symbol": {"col": "Hugo_Symbol"}, "is_neoantigen": {"const": True},
        "cancer_cell_fraction": {"col": "CCF"}, "clonality": {"col": "Clonal"}})
    ip = ImmunizingPeptide(quoted_text="q", section_ref="s", **mapped)
    assert ip.clonality == "clonal" and ip.cancer_cell_fraction == 0.95


# ---- back-compat ----

def test_references_validate_with_empty_mutations():
    for n in ("rojas", "keskin", "li"):
        rec = ExtractedPaper(**json.loads((PKT / "reference_records" / f"{n}_extracted.json").read_text()))
        assert rec.neoantigen_mutations == []
        assert all(p_.clonality is None for p_ in rec.immunizing_peptides)  # Piece B defaults
    assert schema.SCHEMA_VERSION == "2.15.0"
