"""v2.7.1 — MhcAllele accepts a class-II HLA-DP/DQ heterodimer PAIR.

Root-cause fix for the Rojas P28/P29 loss: DP/DQ restriction is an alpha+beta
heterodimer, reported as a pair (HLA-DPA1*01:03/DPB1*02:01). The single-allele
pattern couldn't hold it, so 17 pairs were dropped and atomic add_table took 12
valid single-DRB neighbours with them. The paired branch stores the pair faithfully.
"""
import pathlib
import sys

import pytest

PKT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PKT / "cancervac_packet"))
import schema  # noqa: E402


def _epitope(allele, mhc_class="II"):
    kw = dict(paper_local_id="EPI_X", mhc_class=mhc_class, sequence="SLLQHLIGLAAAA",
              is_neoantigen=True, hla_allele=allele, quoted_text="q", section_ref="s")
    if mhc_class == "I":
        kw["predicted_affinity"] = {"unit": "unknown", "raw": "predicted, no IC50"}
    return schema.MinimalEpitope(**kw)


def test_version_bumped():
    # v2.12.0 (P23) — CD4/CD8 t_cell_subset + evidence mhc_class axes + class-II epitope anchor.
    assert schema.SCHEMA_VERSION == "2.15.0"


@pytest.mark.parametrize("allele", [
    "HLA-DPA1*01:03/DPB1*02:01",   # real Rojas P28
    "HLA-DQA1*01:01/DQB1*05:01",   # real Rojas P29
    "HLA-DPA1*01:03/DPB1*04:02",
    "HLA-DQA1*01:03/DQB1*06:01",
])
def test_real_rojas_pairs_accepted_and_preserved(allele):
    assert _epitope(allele).hla_allele == allele  # stored verbatim, no loss


def test_pair_canonicalized_spaces_and_case():
    # spaces around the slash and lowercase normalize to the canonical slash-joined form
    assert _epitope("HLA-DPA1*01:03 / DPB1*04:02").hla_allele == "HLA-DPA1*01:03/DPB1*04:02"
    assert _epitope("hla-dqa1*01:03/dqb1*06:01").hla_allele == "HLA-DQA1*01:03/DQB1*06:01"


@pytest.mark.parametrize("allele", [
    "HLA-DRB3*02:02", "HLA-DRB1*01:01", "HLA-DRB5*01:02",   # the 12 single-DRB neighbours
])
def test_single_drb_still_valid(allele):
    assert _epitope(allele).hla_allele == allele


def test_class1_single_and_murine_still_valid():
    assert _epitope("HLA-A*02:01", "I").hla_allele == "HLA-A*02:01"
    assert _epitope("H-2Db", "I").hla_allele == "H-2Db"


@pytest.mark.parametrize("bad", [
    "HLA-DPB1*02:01/junk",          # malformed second chain
    "HLA-DPA1*01:03/DPA1*02:01",    # two alpha chains (must be A then B)
    "HLA-DPA1*01:03/DPB1*02:01/X",  # triple
    "DPA1*01:03/DPB1*02:01",        # missing HLA- prefix
])
def test_malformed_pairs_rejected(bad):
    with pytest.raises(Exception):
        _epitope(bad)


def test_pair_valid_on_patient_hla_alleles_list_too():
    # the patient-level list field uses the same MhcAllele type
    p = schema.ExtractedPatient(
        paper_local_id="P1", indication="x", species="human", vaccine_platform="rna",
        n_peptides_synthesized=1, n_peptides_immunogenic=0,
        hla_alleles=["HLA-DPA1*01:03/DPB1*02:01", "HLA-A*02:01"],
        quoted_text="q", section_ref="s")
    assert p.hla_alleles == ["HLA-DPA1*01:03/DPB1*02:01", "HLA-A*02:01"]
