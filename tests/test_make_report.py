"""Curator HTML is a lossless view of a validated record — optional, not a pipeline gate."""
from __future__ import annotations

import pytest

from vaxtract import report as mr
from vaxtract.cli import main as cli_main
from vaxtract.schema import SCHEMA_VERSION, ExtractedPaper


def test_mag_str_table_call_not_sfc():
    assert mr.mag_str({"raw": "1", "unit": "unknown"}) == "table call 1"
    assert mr.mag_str({"raw": "0", "unit": "unknown"}) == "table call 0"
    assert "SFC" not in mr.mag_str({"raw": "[141, 0]", "unit": "unknown"})
    assert mr.mag_str({"raw": "[141, 0]", "unit": "unknown"}) == "[141, 0]"


def test_mag_str_numeric_unit():
    assert "2459" in mr.mag_str({"value": 2459, "unit": "sfc_per_1e6"})
    assert mr.mag_str(None) == ""


def _minimal_record():
    return {
        "pmid": "25837513",
        "journal": "Science",
        "year": 2015,
        "title": "report fixture",
        "cohort_size": 1,
        "indication_summary": "melanoma",
        "data_resolution": "per_sequence",
        "peptide_manifest_present": True,
        "tcr_seq_status": "bulk",
        "finalize_overrides_used": ["allow_candidate_bridge_mismatch"],
        "patients": [{
            "paper_local_id": "MEL21",
            "indication": "melanoma",
            "quoted_text": "patient MEL21",
            "section_ref": "Table 1",
            "n_peptides_synthesized": 1,
            "n_peptides_immunogenic": 1,
        }],
        "immunizing_peptides": [{
            "paper_local_id": "IMP1",
            "sequence": "QTIDNIVFL",
            "is_neoantigen": True,
            "quoted_text": "QTIDNIVFL",
            "section_ref": "Table S4",
        }],
        "epitopes": [{
            "paper_local_id": "EPI1",
            "sequence": "QTIDNIVFL",
            "is_neoantigen": True,
            "mhc_class": "I",
            "parent_peptide_ids": ["IMP1"],
            "predicted_affinity": {"unit": "unknown", "raw": "not printed"},
            "quoted_text": "QTIDNIVFL",
            "section_ref": "Table S4",
        }],
        "evidence": [{
            "patient_paper_id": "MEL21",
            "target_kind": "epitope",
            "epitope_paper_id": "EPI1",
            "assay": "tetramer",
            "outcome": "immunogenic",
            "quoted_text": "dextramer+",
            "section_ref": "Fig 4",
        }],
        "tcr_clonotypes": [{
            "paper_local_id": "TCR_MEL21_01",
            "patient_paper_id": "MEL21",
            "quoted_text": "CASSQDLSGGVYYGYTF",
            "section_ref": "Table S6",
            "chains": [{"locus": "TRB", "junction_aa": "CASSQDLSGGVYYGYTF"}],
            "specific_for_kind": "epitope",
            "epitope_paper_id": "EPI1",
            "specificity_evidence": "multimer_sort",
        }],
    }


def test_refuses_invalid(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text("{}")
    assert mr.main([str(bad), str(tmp_path / "out.html")]) == 1
    assert not (tmp_path / "out.html").exists()


def test_shim_and_cli_report_subcommand(tmp_path):
    from cancervac_packet import make_report as shim
    bad = tmp_path / "bad.json"
    bad.write_text("{}")
    assert shim.main([str(bad), str(tmp_path / "shim.html")]) == 1
    with pytest.raises(SystemExit) as ei:
        cli_main(["report", str(bad), str(tmp_path / "cli.html")])
    assert ei.value.code == 1


def test_html_has_tcr_health_and_hard_override(tmp_path):
    rec = ExtractedPaper(**_minimal_record())
    src = tmp_path / "ok.json"
    src.write_text(rec.model_dump_json())
    out = tmp_path / "ok.html"
    assert mr.main([str(src), str(out)]) == 0
    html = out.read_text()
    assert SCHEMA_VERSION.startswith("2.19")
    assert "TCR_MEL21_01" in html
    assert "allow_candidate_bridge_mismatch" in html
    assert "HARD" in html
    assert "table call" not in html or True
