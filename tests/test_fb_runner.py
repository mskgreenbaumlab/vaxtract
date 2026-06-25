import json
import pathlib

import pytest

import figure_benchmark.model as m
import figure_benchmark.run_benchmark as rb

PDF = pathlib.Path("agentBranch/cancerVacExtrac_claudeAi/data/raw/rojas_37165196/main.pdf")
SRC = pathlib.Path("agentBranch/cancerVacExtrac_claudeAi/data/raw/rojas_37165196/sourceData")
pytestmark = pytest.mark.skipif(not PDF.exists(), reason="Rojas corpus not present")


class FakeProvider:
    name = "fake-opus"

    def extract(self, image_path, prompt, effort):
        # echo Fig 4B truth back as a near-perfect read
        payload = {"series": [{"series_label": "blood_T_pct", "unit": "%", "points": [
            {"key": "Blood (30 weeks)", "value": 1.66},
            {"key": "Liver mass (35 weeks)", "value": 2.48},
            {"key": "Blood (41 weeks)", "value": 0.68},
        ]}]}
        return json.dumps(payload), {"wall_time_s": 0.01, "cost_usd": 0.0}


def test_run_one_panel_scores_and_reports(tmp_path):
    panel = m.Panel(
        panel_id="fig4b", figure_label="Figure 4B", chart_type="bar",
        page=5, bbox_frac=(0.1, 0.1, 0.45, 0.45), truth_file="41586_2023_6063_MOESM7_ESM.xlsx",
        truth_sheet="Figure 4B", truth_adapter="fig4b_bars", align_mode="by_key",
        unit="% all blood T cells",
    )
    results = rb.run_benchmark(
        [panel], provider=FakeProvider(), pdf_path=PDF, source_dir=SRC,
        out_dir=tmp_path, dpi=120,
    )
    assert len(results) == 1
    r = results[0]
    assert r["panel_id"] == "fig4b"
    assert r["metrics"]["per_point"]["n_matched"] == 3
    assert r["metrics"]["per_point"]["median_abs_log10"] < 0.05  # near-perfect
    assert r["metrics_vector"] is None  # fig4b has no axis_calib → vector arm skipped
    assert (tmp_path / "report.json").exists()
    assert (tmp_path / "report.html").exists()
    html = (tmp_path / "report.html").read_text()
    assert "Figure 4B" in html
    assert "vision" in html


def test_run_skips_panels_without_bbox(tmp_path):
    panel = m.Panel(
        panel_id="nobbox", figure_label="x", chart_type="bar", page=0, bbox_frac=None,
        truth_file="41586_2023_6063_MOESM7_ESM.xlsx", truth_sheet="Figure 4B",
        truth_adapter="fig4b_bars",
    )
    results = rb.run_benchmark([panel], provider=FakeProvider(), pdf_path=PDF,
                               source_dir=SRC, out_dir=tmp_path, dpi=120)
    assert results[0]["status"] == "skipped_no_bbox"
