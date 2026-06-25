"""Positional column addressing end-to-end: a Keskin-S5-shaped sheet whose two-row
merged header leaves the key columns nameless/duplicated, so add_table can only reach
them by position. Guards the cost fix (route bulk through add_table instead of the
output-expensive add_entities path)."""
import json
import pathlib
import openpyxl
import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]
REF = json.loads((PKT / "reference_records" / "rojas_extracted.json").read_text())
META = {k: v for k, v in REF.items() if not isinstance(v, list)}


def _write_merged_header_sheet(path):
    """Mimic Keskin Table S5: a leading TITLE row, a group-header row with blanks under
    the merged groups, a sub-label row whose names DUPLICATE ('Sequence' x2), then data.
    The mutant-peptide sequence (col 6) and IMP id (col 7) are only addressable by index."""
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(["Supplementary Table 5. Summary"])                                   # row 0: title
    ws.append(["Patient ID", "Pool", "Gene", "Change", "Len", "HLA allele",
               "Mutated peptide", "", "Immunizing peptide", ""])                    # row 1: groups (blanks)
    ws.append(["", "", "", "", "", "", "Sequence", "ID", "Sequence", "ID"])         # row 2: dup sub-labels
    ws.append([1, "A", "ATP10B", "p.R821Q", 9, "A01:01",
               "QTQKHLDLY", "EPI01", "VPDINMEKK", "IMP03"])                         # row 3: data
    ws.append([1, "B", "DOK7", "p.Q151K", 9, "B08:01",
               "TGKWKLSDL", "EPI02", "LARDIPPAV", "IMP22"])                         # row 4: data
    wb.save(path)


EPITOPE_FIELDS = {
    # header_row=2 is the duplicated sub-label row -> every data column is reached by position
    "sequence": {"col_idx": 6},                          # 'Sequence' is duplicated -> must use index
    "gene_symbol": {"col_letter": "C"},                  # C -> index 2 (clean col, by letter)
    "paper_local_id": {"col_letter": "H"},               # H -> index 7 (the EPI id sub-label is 'ID', dup)
    "parent_peptide_ids": {"template_list": "{#9}"},     # IMP id by index
    "is_neoantigen": {"const": True},
    "mhc_class": {"const": "I"},
    "predicted_affinity": {"const": {"value": None, "unit": "unknown", "raw": "not reported"}},
    "quoted_text": {"const": "epitope listed in Supp Table 5"},
    "section_ref": {"const": "Supplementary Table 5"},
}


def test_read_xlsx_dicts_exposes_positional_keys(tmp_path):
    x = tmp_path / "s5.xlsx"; _write_merged_header_sheet(x)
    headers, rows = agent_core._read_xlsx_dicts(str(x), header_row=2)
    # header row 2 has duplicated 'Sequence' -> name keys collide, but positions don't
    assert rows[0][6] == "QTQKHLDLY"
    assert rows[0][9] == "IMP03"
    assert rows[1][6] == "TGKWKLSDL"


def test_add_table_maps_nameless_columns_by_position(tmp_path):
    out = tmp_path / "r.json"; x = tmp_path / "s5.xlsx"; _write_merged_header_sheet(x)
    agent_core.init_partial(str(out), json.dumps(META))
    mapping = {"fields": EPITOPE_FIELDS}
    ok, msg = agent_core.table_to_entities(str(out), "epitopes", str(x), json.dumps(mapping), header_row=2)
    assert ok, msg
    part = json.loads((tmp_path / "r.json.partial.json").read_text())
    eps = part["epitopes"]
    assert len(eps) == 2
    assert eps[0]["sequence"] == "QTQKHLDLY"
    assert eps[0]["parent_peptide_ids"] == ["IMP03"]
    assert eps[1]["sequence"] == "TGKWKLSDL"
    assert eps[1]["parent_peptide_ids"] == ["IMP22"]


def test_add_table_out_of_range_index_reported(tmp_path):
    out = tmp_path / "r.json"; x = tmp_path / "s5.xlsx"; _write_merged_header_sheet(x)
    agent_core.init_partial(str(out), json.dumps(META))
    mapping = {"fields": {"sequence": {"col_idx": 99}}}
    ok, msg = agent_core.table_to_entities(str(out), "epitopes", str(x), json.dumps(mapping), header_row=2)
    assert not ok and "not found" in msg.lower() and "idx 99" in msg


def test_read_table_positional_filter_and_columns(tmp_path):
    x = tmp_path / "s5.xlsx"; _write_merged_header_sheet(x)
    # filter by the nameless pool column (index 1) and project positional columns
    txt = agent_core.read_table_rows(str(x), header_row=2,
                                     row_filter={"col_idx": 1, "equals": "B"},
                                     columns=[6, 9])
    assert "1/2 data rows" in txt
    assert "TGKWKLSDL" in txt and "IMP22" in txt
    assert "QTQKHLDLY" not in txt           # pool-A row filtered out


def test_read_table_plain_preview_has_no_positional_keys(tmp_path):
    """The plain (no columns/filter) preview must be byte-for-byte unchanged: raw matrix,
    no leaked integer keys."""
    x = tmp_path / "s5.xlsx"; _write_merged_header_sheet(x)
    txt = agent_core.read_table_rows(str(x))
    assert '"QTQKHLDLY"' in txt
    assert '"6":' not in txt                 # no positional keys in the raw preview
