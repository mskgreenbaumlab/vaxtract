"""Deterministic reactivity-matrix loader (`build_reactivity_evidence`)."""
import json
import os
import pathlib

import openpyxl
import pytest

import agent_core as ac
from vaxtract.reactivity_adapter import classify_assay_column, subset_from_sheet

META = json.dumps({"pmid": "33479501", "title": "t", "journal": "j", "year": 2021,
                   "cohort_size": 2, "indication_summary": "melanoma"})
PROV = {"quoted_text": "q", "section_ref": "Suppl Dataset 4b"}
IMP_SEQ = "PQVDGEIPLHRSDRVKVLSIGEGGF"
# Distinct 15-mers nested in IMP_SEQ — one per data row so (patient, epitope, assay)
# keys do not collide. The loader is one-row-per-positive-CELL, not per sheet-row.
ASPS = [IMP_SEQ[i:i + 15] for i in range(8)]
SHEET = "Suppl Dataset 4b. CD4 T cells"
# Four assay columns that vary. Positives: see _react_rows.
HEADER = ["Patient ID", "Gene", "Immunizing peptide", "Assay peptide",
          "W16 ex vivo peptide pulsed", "W16 after pre-stimulation",
          "Yr3 ex vivo", "Yr3 pre-stim"]


def _react_rows():
    # 8 rows x mixed 1/0 across 4 assay cols. 21 positives.
    pattern = [
        (1, 1, 1, 0),
        (1, 0, 1, 0),
        (1, 1, 0, 1),
        (0, 1, 1, 0),
        (1, 1, 0, 0),
        (1, 0, 0, 1),
        (0, 1, 1, 1),
        (1, 1, 1, 1),
    ]
    out = []
    for i, calls in enumerate(pattern):
        pt = str(2 if i < 4 else 3) if i in (0, 4) else ""   # continuation rows
        out.append([pt, "SHANK2", IMP_SEQ, ASPS[i], *calls])
    return out


def _book(tmp_path, header=HEADER, rows=None, sheet=SHEET):
    paper = tmp_path / "paper"
    paper.mkdir(exist_ok=True)
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    ws = wb.create_sheet(sheet[:31])
    ws.append(["Supplementary Dataset 4b. CD4+ T cell reactivity by IFN-g ELISPOT"])
    ws.append(header)
    for r in rows or _react_rows():
        ws.append(r)
    wb.save(paper / "supp4.xlsx")
    return paper


def _record(tmp_path, paper, patients=("Pt2", "Pt3")):
    out = str(tmp_path / "o.json")
    ac.set_source_dir(out, str(paper))
    assert ac.init_partial(out, META)[0]
    for i, pid in enumerate(patients):
        ac.append_section(out, "patients", json.dumps([{
            "paper_local_id": pid, "indication": "melanoma",
            "n_peptides_synthesized": 1, "n_peptides_immunogenic": 1, **PROV}]))
        ac.append_section(out, "immunizing_peptides", json.dumps([{
            "paper_local_id": f"IMP{i+1}", "sequence": IMP_SEQ, "gene_symbol": "SHANK2",
            "is_neoantigen": True, "patient_paper_id": pid, **PROV}]))
    return out


def test_classify_assay_column_reads_hu_style_labels():
    w16 = classify_assay_column("W16 ex vivo peptide pulsed")
    assert w16["timepoint_phase"] == "post_vaccine"
    assert w16["assay_stimulation"] == "ex_vivo"
    assert w16["assay_antigen_format"] == "peptide_pulsed"
    yr = classify_assay_column("Yr3 after pre-stimulation minigene")
    assert yr["timepoint_phase"] == "memory"
    assert yr["assay_stimulation"] == "in_vitro_expanded"
    assert yr["assay_antigen_format"] == "minigene"
    assert subset_from_sheet("Supp11b NetMHCpan2.4 CD4 T cell") == "cd4"
    assert subset_from_sheet("Suppl Dataset 4a. CD8+ T cells") == "cd8"


def test_loader_emits_one_row_per_positive_cell(tmp_path):
    paper = _book(tmp_path)
    from vaxtract import tcr_source_census as census
    n_expected = census.census_tcr_sheets(paper, with_structure=True).reactivity[0].n_expected
    assert n_expected == 21
    out = _record(tmp_path, paper)
    ok, msg = ac.build_reactivity_evidence(out, paper_dir=str(paper))
    assert ok, msg
    rec = json.loads(pathlib.Path(out).read_text() if pathlib.Path(out).exists()
                     else pathlib.Path(out + ".partial.json").read_text())
    ev = rec["evidence"]
    pos = [e for e in ev if e["outcome"] == "immunogenic"]
    neg = [e for e in ev if e["outcome"] == "not_immunogenic"]
    assert len(pos) == 21, msg
    assert len(neg) == 11, msg
    assert all(e["target_kind"] == "epitope" for e in ev)
    assert all(e["t_cell_subset"] == "cd4" for e in ev)
    assert {e["timepoint_phase"] for e in ev} == {"post_vaccine", "memory"}
    assert {e["assay_stimulation"] for e in ev} == {"ex_vivo", "in_vitro_expanded"}
    assert all(e.get("magnitude", {}).get("raw") == "1" for e in pos)
    assert all(e.get("magnitude", {}).get("value") is None for e in pos)
    assert all("grade" not in (e.get("magnitude") or {}) or e["magnitude"].get("grade") is None
               for e in pos)
    assert all(e.get("magnitude", {}).get("raw") == "0" for e in neg)
    # minted the 15-mer once
    assert {e.get("sequence") for e in rec["epitopes"]} >= set(ASPS)


def test_loader_is_idempotent(tmp_path):
    paper = _book(tmp_path)
    out = _record(tmp_path, paper)
    assert ac.build_reactivity_evidence(out, paper_dir=str(paper))[0]
    ok, msg = ac.build_reactivity_evidence(out, paper_dir=str(paper))
    rec = json.loads(pathlib.Path(out + ".partial.json").read_text())
    assert len(rec["evidence"]) == 32
    assert sum(1 for e in rec["evidence"] if e["outcome"] == "immunogenic") == 21
    assert "already present" in msg or rec["evidence"]


def test_unknown_patient_is_skipped_not_invented(tmp_path):
    paper = _book(tmp_path)
    out = _record(tmp_path, paper, patients=("Pt2",))  # Pt3 rows have no patient in record
    ok, msg = ac.build_reactivity_evidence(out, paper_dir=str(paper))
    rec = json.loads(pathlib.Path(out + ".partial.json").read_text())
    pts = {e["patient_paper_id"] for e in rec["evidence"]}
    assert pts <= {"Pt2"}
    assert "patient not in record" in msg


def test_ledger_satisfies_the_reactivity_gate(tmp_path):
    paper = _book(tmp_path)
    out = _record(tmp_path, paper)
    assert ac.build_reactivity_evidence(out, paper_dir=str(paper))[0]
    rec = ac._load_partial(out)[0]
    from vaxtract.schema import ExtractedPaper
    # Gate reads the ledger, not the pydantic object, for coverage.
    paper_obj = type("P", (), {"pmid": "33479501"})()
    gap = ac._reactivity_source_coverage_gap(paper_obj, out)
    assert gap == ""


def test_cd8_mutated_9mer_is_epitope_grain(tmp_path):
    """Hu 4a: CD8 ELISPOT against a 9-mer in 'Mutated peptide Sequence'."""
    paper = tmp_path / "paper"
    paper.mkdir(exist_ok=True)
    imp = IMP_SEQ
    header = ["Patient ID", "Immunizing peptide Sequence", "Mutated peptide Sequence",
              "W16 peptide pulsed"]
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    ws = wb.create_sheet("Suppl Dataset 4a. CD8+ T cells")
    ws.append(["CD8+ T cell reactivity by IFN-g ELISPOT against peptide-pulsed APC"])
    ws.append(header)
    for i in range(12):
        ws.append(["1", imp, imp[i:i + 9], 1 if i < 11 else 0])
    wb.save(paper / "supp4.xlsx")
    out = _record(tmp_path, paper, patients=("Pt1",))
    ok, msg = ac.build_reactivity_evidence(out, paper_dir=str(paper))
    assert ok, msg
    rec = json.loads(pathlib.Path(out + ".partial.json").read_text())
    ev = rec["evidence"]
    pos = [e for e in ev if e["outcome"] == "immunogenic"]
    assert len(pos) == 11, msg
    assert all(e["target_kind"] == "epitope" for e in ev)
    assert all(e["t_cell_subset"] == "cd8" for e in ev)
    assert all(e.get("sequence") and 8 <= len(e["sequence"]) <= 11
               for e in rec["epitopes"] if e["paper_local_id"] in {x["epitope_paper_id"] for x in ev})


def test_inventory_names_the_loader(tmp_path):
    paper = _book(tmp_path)
    txt = ac.survey_sources(str(paper))
    line = next(l for l in txt.splitlines() if "REACTIVITY MATRIX" in l)
    assert "build_reactivity_evidence" in line
    assert "per_cell" in line


ATLAS = pathlib.Path(os.environ["VAXTRACT_ATLAS_33479501"]) if os.environ.get("VAXTRACT_ATLAS_33479501") else pathlib.Path("/nonexistent")
KEEP = pathlib.Path(os.environ["VAXTRACT_KEEP_33479501"]) if os.environ.get("VAXTRACT_KEEP_33479501") else pathlib.Path("/nonexistent")


@pytest.mark.skipif(not ATLAS.is_dir() or not KEEP.is_file(),
                    reason="set VAXTRACT_ATLAS_33479501 and VAXTRACT_KEEP_33479501 to run source-sheet check")
def test_hu_supp11b_loads_eleven_calls_from_the_source_sheet(tmp_path):
    src = json.loads(KEEP.read_text())
    out = str(tmp_path / "hu.json")
    meta = {k: src[k] for k in ("pmid", "title", "journal", "year",
                                "cohort_size", "indication_summary") if k in src}
    ac.set_source_dir(out, str(ATLAS))
    assert ac.init_partial(out, json.dumps(meta))[0]
    ac.append_section(out, "patients", json.dumps([
        {k: p[k] for k in p if k in ("paper_local_id", "indication",
                                     "n_peptides_synthesized", "n_peptides_immunogenic",
                                     "quoted_text", "section_ref")}
        for p in src["patients"]
    ]))
    imps = []
    for p in src["immunizing_peptides"]:
        imps.append({k: p[k] for k in ("paper_local_id", "sequence", "gene_symbol",
                                       "is_neoantigen", "patient_paper_id",
                                       "quoted_text", "section_ref") if k in p})
    # batches of 50
    for i in range(0, len(imps), 40):
        ac.append_section(out, "immunizing_peptides", json.dumps(imps[i:i + 40]))
    ok, msg = ac.build_reactivity_evidence(
        out, paper_dir=str(ATLAS), sheet="Supp11b NetMHCpan2.4 CD4 T cell")
    assert ok, msg
    rec = json.loads(pathlib.Path(out + ".partial.json").read_text())
    n11 = sum(1 for e in rec["evidence"]
              if e.get("outcome") == "immunogenic"
              and ("11b" in (e.get("section_ref") or "").lower()
                   or "Supp11b" in (e.get("section_ref") or "")))
    assert n11 == 11, msg


def test_nd_blank_none_emit_no_row(tmp_path):
    rows = _react_rows()
    # Extra row of n.d. / blank / None must not add evidence.
    rows.append(["3", "SHANK2", IMP_SEQ, ASPS[0], "n.d.", "", None, "nd"])
    paper = _book(tmp_path, rows=rows)
    out = _record(tmp_path, paper)
    ok, msg = ac.build_reactivity_evidence(out, paper_dir=str(paper))
    assert ok, msg
    rec = json.loads(pathlib.Path(out + ".partial.json").read_text())
    ev = rec["evidence"]
    pos = [e for e in ev if e["outcome"] == "immunogenic"]
    neg = [e for e in ev if e["outcome"] == "not_immunogenic"]
    assert len(pos) == 21
    assert len(neg) == 11
    assert all(e["magnitude"]["raw"] == "1" and e["magnitude"].get("value") is None for e in pos)
    assert all(e["magnitude"]["raw"] == "0" for e in neg)


def test_patch_tool_updates_matching_rows_only(tmp_path):
    from vaxtract.reactivity_adapter import add_zero_evidence, patch_table_magnitudes
    from vaxtract.schema import ExtractedPaper

    paper = _book(tmp_path)
    out = _record(tmp_path, paper)
    # Build positives-only by taking the full load then stripping zeros + magnitudes.
    assert ac.build_reactivity_evidence(out, paper_dir=str(paper))[0]
    rec = json.loads(pathlib.Path(out + ".partial.json").read_text())
    rec["evidence"] = [e for e in rec["evidence"] if e["outcome"] == "immunogenic"]
    for e in rec["evidence"]:
        e["magnitude"] = None
    rec["evidence"].append({
        "patient_paper_id": "Pt2",
        "target_kind": "immunizing_peptide",
        "immunizing_peptide_paper_id": "IMP1",
        "epitope_paper_id": None,
        "assay": "ics",
        "outcome": "immunogenic",
        "section_ref": "Fig. 1a not a matrix",
        "quoted_text": "CD4 ICS not on the reactivity sheet",
        "assay_detail": "unrelated",
        "timepoint_phase": "post_vaccine",
        "timepoint_label": "week 16",
        "assay_stimulation": "ex_vivo",
        "assay_antigen_format": "peptide_pulsed",
        "t_cell_subset": "cd4",
        "vaccine_induced": True,
        "magnitude": None,
    })
    n_pos = sum(1 for e in rec["evidence"] if e["outcome"] == "immunogenic"
                and "4b" in (e.get("section_ref") or ""))
    assert n_pos == 21
    report = patch_table_magnitudes(rec, str(paper))
    assert report["n_patched"] == 21
    assert report["n_unmatched"] == 0
    matrix = [e for e in rec["evidence"] if "4b" in (e.get("section_ref") or "")]
    other = [e for e in rec["evidence"] if "4b" not in (e.get("section_ref") or "")]
    assert all(e["magnitude"]["raw"] == "1" and e["magnitude"].get("value") is None for e in matrix)
    assert other[0]["magnitude"] is None
    zreport = add_zero_evidence(rec, str(paper))
    assert zreport["n_zeros_added"] == 11
    assert sum(1 for e in rec["evidence"] if e["outcome"] == "immunogenic") == 22  # 21 + ICS
    assert sum(1 for e in rec["evidence"] if e["outcome"] == "not_immunogenic") == 11
    ExtractedPaper.model_validate(rec)
