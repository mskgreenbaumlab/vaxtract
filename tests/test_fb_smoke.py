import os
import pathlib
import shutil

import pytest

import figure_benchmark.panels as P
import figure_benchmark.run_benchmark as rb

PDF = pathlib.Path("agentBranch/cancerVacExtrac_claudeAi/data/raw/rojas_37165196/main.pdf")
SRC = pathlib.Path("agentBranch/cancerVacExtrac_claudeAi/data/raw/rojas_37165196/sourceData")


@pytest.mark.live_vision
@pytest.mark.skipif(os.environ.get("RUN_LIVE_VISION") != "1",
                    reason="set RUN_LIVE_VISION=1 to hit the real claude CLI")
@pytest.mark.skipif(shutil.which("claude") is None, reason="claude CLI not on PATH")
def test_live_reads_easy_bar_panel(tmp_path):
    panel = P.get("fig4b")
    assert panel.bbox_frac is not None, "Task 9 must set bbox_frac first"
    results = rb.run_benchmark(
        [panel], provider=rb._default_provider(), pdf_path=PDF, source_dir=SRC,
        out_dir=tmp_path, dpi=200,
    )
    r = results[0]
    assert r["status"] == "scored"
    # easy 3-bar panel: vision should match all 3 within ~0.3 log10 (~2x)
    assert r["metrics"]["per_point"]["n_matched"] >= 2
    assert r["metrics"]["per_point"]["median_abs_log10"] < 0.3
