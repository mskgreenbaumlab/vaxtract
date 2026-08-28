"""Union of extra Hu TCR clonotypes onto a loader record."""
from vaxtract.hu_tcr_union import chain_key, union_clonotypes


def _cl(pid, patient, cdr3, *, section="Supplementary Dataset 9", kind="immunizing_peptide",
        imp="IMP1", epi=None):
    return {
        "paper_local_id": pid,
        "patient_paper_id": patient,
        "section_ref": section,
        "quoted_text": "q",
        "specific_for_kind": kind,
        "immunizing_peptide_paper_id": imp if kind == "immunizing_peptide" else None,
        "epitope_paper_id": epi,
        "specificity_evidence": "functional_tcr_clone",
        "chains": [{"locus": "TRB", "junction_aa": cdr3, "quoted_text": "q", "section_ref": section}],
    }


def test_union_adds_only_new_chain_identities():
    base = {
        "immunizing_peptides": [{"paper_local_id": "IMP1", "sequence": "AAAAAAAA"}],
        "epitopes": [],
        "tcr_clonotypes": [_cl("A", "Pt2", "CASSXXX", section="Dataset 8")],
        "neoantigen_tcr_flags": [],
    }
    donor = {
        "immunizing_peptides": [{"paper_local_id": "IMP1", "sequence": "AAAAAAAA"}],
        "epitopes": [],
        "tcr_clonotypes": [
            _cl("A", "Pt2", "CASSXXX", section="Dataset 8"),          # duplicate identity
            _cl("B", "Pt2", "CASSYYY", section="Supplementary Dataset 9"),
        ],
        "neoantigen_tcr_flags": [{
            "patient_paper_id": "Pt2", "target_kind": "immunizing_peptide",
            "immunizing_peptide_paper_id": "IMP1", "tcr_identified": "cd4",
            "quoted_text": "q", "section_ref": "Dataset 9",
        }],
    }
    out, rep = union_clonotypes(base, donor)
    assert rep["clonotypes_added"] == 1
    assert len(out["tcr_clonotypes"]) == 2
    assert out["tcr_seq_status"] == "both"
    assert rep["flags_added"] == 1


def test_dataset10_epitope_is_copied_when_parent_imp_exists():
    base = {
        "immunizing_peptides": [{"paper_local_id": "IMP_MLL", "sequence": "SRLQTRKNKKLALSSTPSNIAPSD"}],
        "epitopes": [],
        "tcr_clonotypes": [],
        "neoantigen_tcr_flags": [],
    }
    donor = {
        "epitopes": [{
            "paper_local_id": "EPII_Pt6_MLL_a", "sequence": "KNKKLALSSTPSNIA",
            "mhc_class": "II", "is_neoantigen": True,
            "parent_peptide_ids": ["IMP_MLL"],
            "quoted_text": "q", "section_ref": "Supplementary Dataset 10",
        }],
        "tcr_clonotypes": [_cl("TUM01", "Pt6", "CASSZZZ",
                               section="Supplementary Dataset 10; Fig. 4d",
                               kind="epitope", imp=None, epi="EPII_Pt6_MLL_a")],
        "neoantigen_tcr_flags": [],
        "immunizing_peptides": [],
    }
    out, rep = union_clonotypes(base, donor)
    assert rep["clonotypes_added"] == 1
    assert "EPII_Pt6_MLL_a" in {e["paper_local_id"] for e in out["epitopes"]}
    assert out["tcr_clonotypes"][0]["epitope_paper_id"] == "EPII_Pt6_MLL_a"


def test_designated_hu_keep_file_has_loader_evidence_and_unioned_tcr():
    import json, pathlib
    p = (pathlib.Path(__file__).resolve().parents[1]
         / "outputs/tcr_hu_keep_2026-08-25/33479501_extracted.json")
    if not p.is_file():
        return
    rec = json.loads(p.read_text())
    ev = rec["evidence"]
    n4a = sum(1 for e in ev if "4a" in (e.get("section_ref") or ""))
    n4b = sum(1 for e in ev if "4b" in (e.get("section_ref") or ""))
    n11 = sum(1 for e in ev if "11b" in (e.get("section_ref") or "").lower())
    assert n4a == 17 and n4b == 279 and n11 == 11
    assert all(e.get("target_kind") == "epitope" for e in ev)
    assert len(rec["tcr_clonotypes"]) == 333
    assert rec.get("tcr_seq_status") == "both"
    assert rec.get("finalize_overrides_used") == ["allow_missing_magnitudes"]


def test_chain_key_ignores_paper_local_id_when_cdr3_present():
    a = _cl("ONE", "Pt2", "CASSABC")
    b = _cl("TWO", "Pt2", "CASSABC")
    assert chain_key(a) == chain_key(b)
