"""read_docx — .docx supplement reader (tables + prose). Fixtures are built in-test with
python-docx so there's no dependency on gitignored data/raw. Covers the table summary, caption
association, specific-table preview (reusing read_table's header-skip/row_filter/columns/byte-cap
via _render_matrix), underline reveal, prose paging, and out-of-range handling."""
import json
import pathlib

import docx
import pytest

import agent_core


def _doc(path, captions_and_rows, prose=None):
    """captions_and_rows: list of (caption, rows[][]) -> paragraph caption then a table."""
    d = docx.Document()
    for cap, rows in captions_and_rows:
        if cap:
            d.add_paragraph(cap)
        t = d.add_table(rows=len(rows), cols=len(rows[0]))
        for r, row in enumerate(rows):
            for c, val in enumerate(row):
                t.rows[r].cells[c].text = str(val)
    for para in (prose or []):
        d.add_paragraph(para)
    d.save(str(path))
    return str(path)


# ---- summary + caption association ----

def test_summary_lists_tables_dims_and_captions(tmp_path):
    p = _doc(tmp_path / "s.docx", [
        ("Supplementary Table S6. New neoantigen mutations in recurrent tumor",
         [["Gene", "Mutation"], ["TP53", "R175H"], ["KRAS", "G12D"]]),
        ("Supplementary Table S1. Baseline characteristics", [["Age", "Sex"], ["61", "M"]]),
    ])
    out = agent_core.read_docx_from(p)  # no table_index -> summary
    assert "2 table(s)" in out
    assert "S6" in out and "S1" in out and "table_index" in out
    assert "3x2" in out and "2x2" in out  # dims rows×cols


def test_no_tables_doc_points_to_prose(tmp_path):
    p = _doc(tmp_path / "n.docx", [], prose=["Supplementary methods.", "Patients were…"])
    out = agent_core.read_docx_from(p)
    assert "0 tables" in out and "text_offset" in out


# ---- specific table ----

def test_read_specific_table_preview(tmp_path):
    p = _doc(tmp_path / "t.docx", [
        ("Table S6.", [["Gene", "Mutation"], ["TP53", "R175H"], ["KRAS", "G12D"]])])
    out = agent_core.read_docx_from(p, table_index=0)
    assert "showing" in out and "table[0]" in out
    body = json.loads(out[out.index("\n", out.index(":\n")) + 1:]) if False else json.loads(out.split(":\n", 1)[1])
    assert ["Gene", "Mutation"] in body and ["TP53", "R175H"] in body


def test_table_index_out_of_range(tmp_path):
    p = _doc(tmp_path / "t.docx", [("c", [["a"], ["b"]])])
    out = agent_core.read_docx_from(p, table_index=5)
    assert "out of range" in out


def test_columns_projection_with_title_row_skip(tmp_path):
    # first row is a single-cell title -> _header_index should skip it and use row 1 as header
    p = _doc(tmp_path / "t.docx", [
        ("cap", [["Table S6. neoantigens", ""], ["Gene", "Mutation"],
                 ["TP53", "R175H"], ["KRAS", "G12D"]])])
    out = agent_core.read_docx_from(p, table_index=0, columns=["Gene"])
    assert "data rows" in out
    body = json.loads(out.split(":\n", 1)[1])
    assert {"Gene": "TP53"} in body and {"Gene": "KRAS"} in body


def test_row_filter(tmp_path):
    p = _doc(tmp_path / "t.docx", [
        ("cap", [["Gene", "Mutation"], ["TP53", "R175H"], ["KRAS", "G12D"]])])
    out = agent_core.read_docx_from(p, table_index=0, row_filter={"col": "Gene", "equals": "KRAS"})
    body = json.loads(out.split(":\n", 1)[1])
    assert len(body) == 1 and body[0]["Mutation"] == "G12D"


# ---- underline reveal (minimal-epitope marking) ----

def test_underline_reveal(tmp_path):
    d = docx.Document()
    d.add_paragraph("Table S2.")
    t = d.add_table(rows=2, cols=1)
    t.rows[0].cells[0].text = "Peptide"
    cell = t.rows[1].cells[0]
    cell.paragraphs[0].text = ""
    r1 = cell.paragraphs[0].add_run("QLAS")
    r2 = cell.paragraphs[0].add_run("STYTAYIV"); r2.font.underline = True
    r3 = cell.paragraphs[0].add_run("GYV")
    p = str(tmp_path / "u.docx"); d.save(p)
    out = agent_core.read_docx_from(p, table_index=0, underline=True)
    assert "<u>STYTAYIV</u>" in out


# ---- prose paging ----

def test_prose_paging(tmp_path):
    long = "X" * 100
    p = _doc(tmp_path / "p.docx", [], prose=[long, long, long])  # ~302 chars joined
    win = agent_core.read_docx_from(p, text_offset=0, max_chars=120)
    assert "docx prose 0-120" in win and "text_offset=120" in win
    rest = agent_core.read_docx_from(p, text_offset=120, max_chars=1000)
    assert "120-" in rest and "text_offset=" not in rest.split("chars:\n", 1)[1]  # no further page


# ---- discovery: .docx is now an accepted source suffix ----

def test_docx_is_a_discovered_suffix():
    src = (pathlib.Path(__file__).resolve().parents[1] / "vaxtract" / "extraction_agent.py").read_text()
    assert '".pdf", ".xlsx", ".docx"' in src
