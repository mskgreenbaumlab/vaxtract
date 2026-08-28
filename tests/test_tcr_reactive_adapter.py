"""Deterministic paired single-cell reactive-TCR loader (stability lab, 2026-07-06) — the last
non-deterministic lane. Hermetic: synthesizes the Keskin MOESM7/10 sheet shape (two-row header, one row
per cell, unique CDR3a/CDR3b pair = one clonotype) and checks grouping, the min_cells expansion cutoff,
schema-validity, and the merge/idempotency of the agent_core tool path.
"""
import json
import pathlib
import tempfile

import openpyxl

import agent_core as ac
from vaxtract import tcr_reactive_adapter as tra
from cancervac_packet.schema import TcrClonotype

META = json.dumps({"pmid": "30568305", "title": "t", "journal": "j", "year": 2018,
                   "cohort_size": 8, "indication_summary": "melanoma"})

HEADER = ["Cell designation", "TRAV", "TRAJ", "Alpha CDR3 DNA seq", "Alpha CDR3 amino acid seq",
          "TRBV", "TRBJ", "TRBC", "Beta CDR3 DNA seq", "Beta CDR3 amino acid seq", "Clone", "Clonotype"]
# clonotype A appears in 3 cells (expanded); B in 1 cell (singleton)
CELLS = [
    ["c1", "TRAV1", "TRAJ1", "gca", "CAVNDYKLSF", "TRBV1", "TRBJ1", "TRBC1", "gcc", "CASSLGQPQHF", "1", "A"],
    ["c2", "TRAV1", "TRAJ1", "gca", "CAVNDYKLSF", "TRBV1", "TRBJ1", "TRBC1", "gcc", "CASSLGQPQHF", "1", "A"],
    ["c3", "TRAV1", "TRAJ1", "gca", "CAVNDYKLSF", "TRBV1", "TRBJ1", "TRBC1", "gcc", "CASSLGQPQHF", "1", "A"],
    ["c4", "TRAV2", "TRAJ2", "gcg", "CAGPRDSNYQLIW", "TRBV2", "TRBJ2", "TRBC2", "gct", "CASSQDLGEKLFF", "2", "B"],
]


def _sheet():
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "(a)Pt8 mutSHANK2-react Tcells"
    ws.append(["Pt8 mutantSHANK2 reactive T cells"])           # title row (row 0)
    ws.append(HEADER)                                          # column header (row 1)
    for c in CELLS:
        ws.append(c)
    p = pathlib.Path(tempfile.mkdtemp()) / "tcr.xlsx"
    wb.save(p)
    return str(p), ws.title


def test_groups_cells_into_clonotypes_by_cdr3_pair():
    path, sheet = _sheet()
    cl = tra.load_reactive_clonotypes(path, sheet, patient_paper_id="Pt8", antigen="mutant SHANK2",
                                      section_ref="S", min_cells=1)
    assert len(cl) == 2 and all(isinstance(c, TcrClonotype) for c in cl)
    a = next(c for c in cl if any(ch.junction_aa == "CAVNDYKLSF" for ch in c.chains))
    assert a.chain_pairing == "paired_ab"
    assert {ch.locus for ch in a.chains} == {"TRA", "TRB"}
    assert a.receptor_modality == "single_cell_tcr_seq" and a.needs_review is False


def test_min_cells_drops_singletons():
    path, sheet = _sheet()
    expanded = tra.load_reactive_clonotypes(path, sheet, patient_paper_id="Pt8", antigen="x",
                                            section_ref="S", min_cells=2)
    assert len(expanded) == 1                                  # only clonotype A (3 cells) survives


def test_deterministic():
    path, sheet = _sheet()
    a = tra.load_reactive_clonotypes(path, sheet, patient_paper_id="Pt8", antigen="x", section_ref="S")
    b = tra.load_reactive_clonotypes(path, sheet, patient_paper_id="Pt8", antigen="x", section_ref="S")
    assert [c.model_dump() for c in a] == [c.model_dump() for c in b]


def test_tool_path_merges_and_is_idempotent():
    path, sheet = _sheet()
    out = str(pathlib.Path(tempfile.mkdtemp()) / "o.json")
    ok, _ = ac.init_partial(out, META)
    assert ok
    ok, msg = ac.load_reactive_tcr_into(out, path, sheet, "Pt8", "mutant SHANK2", min_cells=1)
    assert ok and "loaded 2" in msg
    rec = json.loads(pathlib.Path(out + ".partial.json").read_text())
    assert len(rec["tcr_clonotypes"]) == 2 and rec["tcr_seq_status"] == "single_cell"
    ok2, msg2 = ac.load_reactive_tcr_into(out, path, sheet, "Pt8", "mutant SHANK2", min_cells=1)
    assert ok2 and "loaded 0" in msg2                          # chain-key dedupe -> nothing new
    rec2 = json.loads(pathlib.Path(out + ".partial.json").read_text())
    assert len(rec2["tcr_clonotypes"]) == 2
