"""TCR sheet-scope record filter (signed 2026-08-25)."""
import json
import pathlib

from vaxtract.tcr_scope_filter import is_out_of_scope_clonotype, strip_out_of_scope_clonotypes

_KESKIN = (pathlib.Path(__file__).resolve().parents[1]
           / "outputs/classii_keskin_tcrscope_2026-08-25/30568305_extracted.json")
_REFINAL = (pathlib.Path(__file__).resolve().parents[1]
            / "outputs/classii_keskin_opus5_1m_2026-08-22_refinal/30568305_extracted.json")


def test_table_11a_locator_is_out_of_scope():
    assert is_out_of_scope_clonotype({"section_ref": "Supplementary Table 11a"})
    assert is_out_of_scope_clonotype({"section_ref": "MOESM10 (a) Tumor-associated Tcells"})
    assert not is_out_of_scope_clonotype({"section_ref": "Supplementary Table 11b"})
    assert not is_out_of_scope_clonotype({"section_ref": "Supplementary Table 8a"})


def test_quoted_text_til_does_not_exclude_an_antigen_sorted_row():
    # locators only — a reactive clonotype's quote can mention TIL without being the TIL sheet
    assert not is_out_of_scope_clonotype({
        "section_ref": "Supplementary Table 11b",
        "quoted_text": "tumour-associated comparison; TIL repertoire as denominator",
    })


def test_strip_drops_only_the_11a_rows():
    rec = {
        "immunizing_peptides": [{"sequence": "AAAAAAAAB"}],
        "tcr_clonotypes": [
            {"section_ref": "Supplementary Table 11a", "paper_local_id": "a"},
            {"section_ref": "Supplementary Table 11b", "paper_local_id": "b"},
            {"section_ref": "Supplementary Table 8a", "paper_local_id": "c"},
        ],
    }
    out, dropped = strip_out_of_scope_clonotypes(rec)
    assert [c["paper_local_id"] for c in dropped] == ["a"]
    assert [c["paper_local_id"] for c in out["tcr_clonotypes"]] == ["b", "c"]
    assert rec["tcr_clonotypes"][0]["paper_local_id"] == "a"  # source not mutated in place of list? copy.copy
    assert out["immunizing_peptides"] == rec["immunizing_peptides"]


def test_keskin_tcrscope_keepfile_is_436_with_no_11a():
    if not _KESKIN.is_file():
        return
    rec = json.loads(_KESKIN.read_text())
    assert len(rec["tcr_clonotypes"]) == 436
    assert all("11a" not in (c.get("section_ref") or "").lower() for c in rec["tcr_clonotypes"])
    if _REFINAL.is_file():
        src = json.loads(_REFINAL.read_text())
        assert len(src["tcr_clonotypes"]) == 664
        assert len(src["immunizing_peptides"]) == len(rec["immunizing_peptides"])
        assert len(src["epitopes"]) == len(rec["epitopes"])
        assert len(src["evidence"]) == len(rec["evidence"])
