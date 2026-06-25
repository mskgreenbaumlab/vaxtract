"""read_table_rows inspection args (row_filter / columns).

These let the agent slice a big table IN-TOOL instead of grepping the SDK's spilled
read_table result file (RUNS.md run 4: 31 host Greps). The host file tools are also
hard-denied in extraction_agent.py; this is the sanctioned replacement.
"""
import json
import pathlib

import openpyxl

import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]


def _write_sheet(path):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(["LID", "SEQ", "GENE", "RESP"])
    ws.append(["P1_N1", "SLLQHLIGL", "TP53", "De novo response"])
    ws.append(["P1_N2", "KIIGNRTLA", "KRAS", "No data"])
    ws.append(["P2_N1", "RTLAKIIGN", "ARHGAP35", "De novo response"])
    ws.append([None, None, None, None])  # fully-empty row -> skipped
    wb.save(path)


def test_no_filter_returns_raw_matrix_unchanged(tmp_path):
    x = tmp_path / "t.xlsx"; _write_sheet(x)
    out = agent_core.read_table_rows(str(x))
    # header row + 3 data rows present in the matrix preview (back-compat)
    body = json.loads(out.split(":\n", 1)[1])
    assert body[0] == ["LID", "SEQ", "GENE", "RESP"]
    assert len(body) == 4  # header + 3 (empty row dropped by openpyxl values pass)


def test_row_filter_in_keeps_matches_and_reports_total(tmp_path):
    x = tmp_path / "t.xlsx"; _write_sheet(x)
    out = agent_core.read_table_rows(
        str(x), row_filter={"col": "RESP", "in": ["De novo response"]})
    rows = json.loads(out.split(":\n", 1)[1])
    assert {r["LID"] for r in rows} == {"P1_N1", "P2_N1"}
    assert "2/3 data rows" in out  # matched/total surfaced so the agent knows the size


def test_row_filter_equals_and_columns_projection(tmp_path):
    x = tmp_path / "t.xlsx"; _write_sheet(x)
    out = agent_core.read_table_rows(
        str(x), row_filter={"col": "RESP", "equals": "No data"}, columns=["LID", "GENE"])
    rows = json.loads(out.split(":\n", 1)[1])
    assert rows == [{"LID": "P1_N2", "GENE": "KRAS"}]


def test_columns_only_projects_all_rows(tmp_path):
    x = tmp_path / "t.xlsx"; _write_sheet(x)
    out = agent_core.read_table_rows(str(x), columns=["GENE"])
    rows = json.loads(out.split(":\n", 1)[1])
    assert rows == [{"GENE": "TP53"}, {"GENE": "KRAS"}, {"GENE": "ARHGAP35"}]


def test_unknown_filter_column_reports_headers(tmp_path):
    x = tmp_path / "t.xlsx"; _write_sheet(x)
    out = agent_core.read_table_rows(str(x), row_filter={"col": "NOPE", "not_empty": True})
    assert "not found" in out and "LID" in out


def test_unknown_projection_column_reports_headers(tmp_path):
    x = tmp_path / "t.xlsx"; _write_sheet(x)
    out = agent_core.read_table_rows(str(x), columns=["GENE", "MISSING"])
    assert "['MISSING']" in out and "not found" in out


def test_bad_filter_operator_surfaces_dsl_error(tmp_path):
    x = tmp_path / "t.xlsx"; _write_sheet(x)
    out = agent_core.read_table_rows(str(x), row_filter={"col": "RESP"})
    assert "exactly one operator" in out
