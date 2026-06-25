import pytest
from table_map import (apply_mapping, render_template, apply_filter, missing_columns,
                       col_letter_to_index)

ROW = {"Patient number": 1, "Neoantigen number": 7, "Gene": "KRAS",
       "Substitution": "G12V", "Mutant Neoantigen Sequence": "AAAVGVGKSAL",
       "MHC-II Mutant Epitope": "", "ELISpot Response": "De novo response"}

# A row keyed by BOTH header name and 0-based position, as agent_core builds for sheets
# with merged/blank/duplicated headers. Cols: 0=Patient 5=HLA 6=mut-seq 7=affinity 11=IMP-ID.
POS_ROW = {0: 1, 5: "B08:01", 6: "TGKWKLSDL", 7: 412, 11: "14362-002-IMP22", 12: "",
           "Patient ID": 1, "HLA allele": "B08:01"}


def test_render_template_interpolates_columns():
    assert render_template("P{Patient number}-N{Neoantigen number}", ROW) == "P1-N7"
    assert render_template("{Gene}:{Substitution}", ROW) == "KRAS:G12V"


def test_apply_mapping_col_const_template_and_list():
    fields = {
        "paper_local_id": {"template": "P{Patient number}-N{Neoantigen number}"},
        "sequence": {"col": "Mutant Neoantigen Sequence"},
        "gene_symbol": {"col": "Gene"},
        "is_neoantigen": {"const": True},
        "mutation": {"template": "{Gene}:{Substitution}"},
        "parent_peptide_ids": {"template_list": "P{Patient number}-N{Neoantigen number}"},
        "predicted_affinity": {"const": {"value": None, "unit": "unknown"}},
    }
    out = apply_mapping(ROW, fields)
    assert out["paper_local_id"] == "P1-N7"
    assert out["sequence"] == "AAAVGVGKSAL"
    assert out["is_neoantigen"] is True
    assert out["mutation"] == "KRAS:G12V"
    assert out["parent_peptide_ids"] == ["P1-N7"]
    assert out["predicted_affinity"] == {"value": None, "unit": "unknown"}


def test_apply_mapping_omits_empty_col_values():
    out = apply_mapping(ROW, {"x": {"col": "MHC-II Mutant Epitope"}})
    assert "x" not in out  # empty cell -> field omitted (so optionals stay unset)


def test_apply_mapping_raises_on_unknown_rule():
    with pytest.raises(ValueError):
        apply_mapping(ROW, {"x": {"bogus": 1}})


def test_apply_filter_equals_in_not_empty():
    rows = [{"r": "De novo response", "e": "x"}, {"r": "no", "e": ""}]
    ok, kept = apply_filter(rows, {"col": "r", "in": ["De novo response"]})
    assert ok and kept == [rows[0]]
    ok, kept = apply_filter(rows, {"col": "e", "not_empty": True})
    assert ok and kept == [rows[0]]
    ok, kept = apply_filter(rows, {"col": "r", "equals": "no"})
    assert ok and kept == [rows[1]]


def test_apply_filter_rejects_unknown_operator():
    ok, msg = apply_filter([{"r": 1}], {"col": "r", "weird": 1})
    assert not ok and "operator" in msg.lower()


def test_missing_columns_detects_unknown_headers():
    headers = ["Gene", "Patient number"]
    fields = {"a": {"col": "Gene"}, "b": {"template": "{Patient number}-{Nope}"}}
    assert missing_columns(fields, {"col": "AlsoNope", "in": [1]}, headers) == ["AlsoNope", "Nope"]
    assert missing_columns({"a": {"col": "Gene"}}, None, headers) == []


# ---- positional (by-index / by-letter) column addressing --------------------

def test_col_letter_to_index():
    assert col_letter_to_index("A") == 0
    assert col_letter_to_index("L") == 11
    assert col_letter_to_index("aa") == 26   # case-insensitive
    assert col_letter_to_index("AB") == 27
    with pytest.raises(ValueError):
        col_letter_to_index("12")


def test_apply_mapping_by_index_and_letter():
    # the IMP-ID and affinity columns are nameless -> reach them by position
    fields = {
        "parent_imp_id": {"col_idx": 11},
        "affinity_nM": {"col_letter": "H"},        # H -> index 7
        "sequence": {"col_idx": 6},
        "hla": {"col": "HLA allele"},              # name still works in the same mapping
    }
    out = apply_mapping(POS_ROW, fields)
    assert out == {"parent_imp_id": "14362-002-IMP22", "affinity_nM": 412,
                   "sequence": "TGKWKLSDL", "hla": "B08:01"}


def test_apply_mapping_omits_empty_positional_cell():
    out = apply_mapping(POS_ROW, {"x": {"col_idx": 12}})  # blank cell
    assert "x" not in out


def test_apply_mapping_rejects_bad_col_idx():
    with pytest.raises(ValueError):
        apply_mapping(POS_ROW, {"x": {"col_idx": "notint"}})


def test_template_index_and_letter_tokens():
    assert render_template("{#0}:{#11}", POS_ROW) == "1:14362-002-IMP22"
    assert render_template("{@G}", POS_ROW) == "TGKWKLSDL"   # G -> index 6
    assert render_template("{Patient ID}-{#5}", POS_ROW) == "1-B08:01"


def test_apply_filter_by_index_and_letter():
    rows = [{0: "A", 5: "x"}, {0: "B", 5: ""}]
    ok, kept = apply_filter(rows, {"col_idx": 0, "in": ["A"]})
    assert ok and kept == [rows[0]]
    ok, kept = apply_filter(rows, {"col_letter": "F", "not_empty": True})  # F -> 5
    assert ok and kept == [rows[0]]


def test_missing_columns_flags_out_of_range_index():
    headers = ["Gene", "Patient number"]           # width 2 -> valid indices 0,1
    fields = {"a": {"col_idx": 1}, "b": {"col_idx": 9}, "c": {"col": "Gene"}}
    miss = missing_columns(fields, None, headers)
    assert miss == ["idx 9 (width 2)"]             # in-range index + valid name not flagged
    # template positional token out of range is caught too
    assert missing_columns({"a": {"template": "{#5}"}}, None, headers) == ["idx 5 (width 2)"]
