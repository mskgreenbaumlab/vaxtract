"""H3 clonotype row gate: schema-forced needs_review only."""
import json
import pathlib

from vaxtract.tcr_h3_gate import (
    apply_h3,
    apply_h3_to_clonotypes,
    clonotype_h3_needs_review,
    flag_h3_needs_review,
)

_HU = (pathlib.Path(__file__).resolve().parents[1]
       / "outputs/tcr_hu_keep_2026-08-25/33479501_extracted.json")
_KESKIN = (pathlib.Path(__file__).resolve().parents[1]
           / "outputs/classii_keskin_tcrscope_2026-08-25/30568305_extracted.json")


def test_inferred_association_stays_flagged():
    assert clonotype_h3_needs_review({
        "specificity_evidence": "inferred_association",
        "observations": [],
        "observations_identity": "not_asserted",
    })


def test_functional_clone_with_no_obs_is_public():
    assert not clonotype_h3_needs_review({
        "specificity_evidence": "functional_tcr_clone",
        "observations": [],
        "observations_identity": "not_asserted",
        "needs_review": True,
    })


def test_not_established_adapter_row_is_public():
    assert not clonotype_h3_needs_review({
        "specificity_evidence": "not_established",
        "observations": [],
        "observations_identity": "not_asserted",
        "needs_review": True,
    })


def test_grouped_obs_without_paper_identity_stay_flagged():
    assert clonotype_h3_needs_review({
        "specificity_evidence": "not_established",
        "observations_identity": "not_asserted",
        "observations": [{"frequency_value": 0.1}, {"frequency_value": 0.2}],
    })


def test_paper_explicit_tracks_are_public():
    assert not clonotype_h3_needs_review({
        "specificity_evidence": "not_established",
        "observations_identity": "paper_explicit",
        "observations": [{"frequency_value": 0.1}, {"frequency_value": 0.2}],
    })


def test_apply_clears_only_unstamped_schema_rows_and_copies():
    rec = {
        "immunizing_peptides": [{"sequence": "AAAAAAAAB"}],
        "tcr_clonotypes": [
            {"paper_local_id": "a", "specificity_evidence": "inferred_association",
             "needs_review": True, "observations": []},
            {"paper_local_id": "b", "specificity_evidence": "functional_tcr_clone",
             "needs_review": True, "observations": [],
             "observations_identity": "not_asserted"},
        ],
    }
    out, cleared = apply_h3_to_clonotypes(rec)
    assert [c["paper_local_id"] for c in cleared] == ["b"]
    assert out["tcr_clonotypes"][0]["needs_review"] is True
    assert out["tcr_clonotypes"][1]["needs_review"] is False
    assert rec["tcr_clonotypes"][1]["needs_review"] is True
    assert out["immunizing_peptides"] == rec["immunizing_peptides"]


def test_hu_keep_clears_exactly_the_20_dataset_9_clones():
    if not _HU.is_file():
        return
    rec = json.loads(_HU.read_text())
    out, cleared = apply_h3_to_clonotypes(rec)
    assert len(cleared) == 20
    assert all("Dataset 9" in (c.get("section_ref") or "") for c in cleared)
    assert all(c.get("specificity_evidence") == "functional_tcr_clone" for c in cleared)
    public = [c for c in out["tcr_clonotypes"] if not c.get("needs_review")]
    assert len(public) == 333
    assert len(out["tcr_clonotypes"]) == 333
    assert len(out["evidence"]) == len(rec["evidence"]) == 307


def test_keskin_tcrscope_clears_all_436_antigen_sorted():
    if not _KESKIN.is_file():
        return
    rec = json.loads(_KESKIN.read_text())
    out, cleared = apply_h3_to_clonotypes(rec)
    assert len(cleared) == 436
    assert len(out["tcr_clonotypes"]) == 436
    assert all(not c.get("needs_review") for c in out["tcr_clonotypes"])
    assert all(c.get("specificity_evidence") == "not_established" for c in out["tcr_clonotypes"])
    assert len(out["evidence"]) == len(rec["evidence"])


def test_adapter_stamped_flag_is_public():
    assert not flag_h3_needs_review({
        "tcr_identified": "cd8",
        "method_label": "add_reactive_tcr (single-cell)",
        "needs_review": True,
    })


def test_keskin_tcrscope_clears_four_stamped_flags():
    if not _KESKIN.is_file():
        return
    rec = json.loads(_KESKIN.read_text())
    out, cl_cleared, fl_cleared = apply_h3(rec)
    assert len(cl_cleared) == 436
    assert len(fl_cleared) == 4
    assert len(out["neoantigen_tcr_flags"]) == 7
    assert all(not f.get("needs_review") for f in out["neoantigen_tcr_flags"])
    assert {f.get("tcr_identified") for f in fl_cleared} == {"cd8", "cd4"}
