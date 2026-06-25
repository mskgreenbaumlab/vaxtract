import pathlib
import pytest

import figure_benchmark.render_panel as rp

PDF = pathlib.Path(
    "agentBranch/cancerVacExtrac_claudeAi/data/raw/rojas_37165196/main.pdf"
)
pytestmark = pytest.mark.skipif(not PDF.exists(), reason="Rojas main.pdf not present")


def test_render_full_page_returns_image(tmp_path):
    img = rp.render_page(PDF, page=1, dpi=150)
    assert img.width > 500 and img.height > 500


def test_crop_frac_returns_subregion():
    img = rp.render_page(PDF, page=1, dpi=150)
    crop = rp.crop_frac(img, (0.0, 0.0, 0.5, 0.5))
    assert crop.width == img.width // 2
    assert crop.height == img.height // 2


def test_render_panel_writes_cached_png(tmp_path):
    from figure_benchmark.model import Panel
    panel = Panel(
        panel_id="probe", figure_label="probe", chart_type="bar",
        page=1, bbox_frac=(0.05, 0.05, 0.45, 0.45),
        truth_file="x", truth_sheet="x", truth_adapter="x",
    )
    out1 = rp.render_panel(panel, PDF, out_dir=tmp_path, dpi=150)
    assert out1.exists() and out1.suffix == ".png"
    mtime1 = out1.stat().st_mtime_ns
    out2 = rp.render_panel(panel, PDF, out_dir=tmp_path, dpi=150)  # cache hit
    assert out2 == out1
    assert out2.stat().st_mtime_ns == mtime1  # not rewritten


def test_render_panel_requires_bbox():
    from figure_benchmark.model import Panel
    panel = Panel(
        panel_id="nobbox", figure_label="x", chart_type="bar", page=1,
        bbox_frac=None, truth_file="x", truth_sheet="x", truth_adapter="x",
    )
    with pytest.raises(ValueError):
        rp.render_panel(panel, PDF, out_dir="/tmp", dpi=150)


def test_render_panel_uses_per_panel_dpi_when_dpi_arg_omitted(tmp_path):
    from figure_benchmark.model import Panel
    # render_dpi=600 should produce a larger crop than the 150-dpi default panel
    big = Panel(panel_id="big", figure_label="x", chart_type="scatter_log", page=1,
                bbox_frac=(0.1, 0.1, 0.4, 0.4), truth_file="x", truth_sheet="x",
                truth_adapter="x", render_dpi=300)
    small = Panel(panel_id="small", figure_label="x", chart_type="bar", page=1,
                  bbox_frac=(0.1, 0.1, 0.4, 0.4), truth_file="x", truth_sheet="x",
                  truth_adapter="x", render_dpi=100)
    from PIL import Image
    big_png = rp.render_panel(big, PDF, out_dir=tmp_path)   # no dpi arg → uses 300
    small_png = rp.render_panel(small, PDF, out_dir=tmp_path)  # no dpi arg → uses 100
    assert Image.open(big_png).width > Image.open(small_png).width
