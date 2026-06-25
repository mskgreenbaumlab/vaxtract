import math

import pytest

import figure_benchmark.model as m
import figure_benchmark.vector_extract as ve


def test_interp_axis_linear():
    ax = {"p0_pt": 0.0, "v0": 0.0, "p1_pt": 100.0, "v1": 50.0, "log": False}
    assert ve.interp_axis(0.0, ax) == pytest.approx(0.0)
    assert ve.interp_axis(50.0, ax) == pytest.approx(25.0)
    assert ve.interp_axis(100.0, ax) == pytest.approx(50.0)


def test_interp_axis_log():
    # pixels 0..100 map to values 1..1000 on a log axis → midpoint = 10^1.5 ≈ 31.6
    ax = {"p0_pt": 0.0, "v0": 1.0, "p1_pt": 100.0, "v1": 1000.0, "log": True}
    assert ve.interp_axis(0.0, ax) == pytest.approx(1.0)
    assert ve.interp_axis(50.0, ax) == pytest.approx(10 ** 1.5, rel=1e-6)
    assert ve.interp_axis(100.0, ax) == pytest.approx(1000.0)


def test_interp_axis_log_rejects_nonpositive_anchor():
    # a log axis with a zero/negative anchor is a calibration error → clear ValueError,
    # not a deep math-domain traceback
    ax = {"p0_pt": 0.0, "v0": 0.0, "p1_pt": 100.0, "v1": 1000.0, "log": True}
    with pytest.raises(ValueError):
        ve.interp_axis(50.0, ax)


def test_assign_series_by_x_bucket():
    buckets = [{"label": "early", "x_lo": 0.0, "x_hi": 10.0},
               {"label": "peak", "x_lo": 10.0, "x_hi": 20.0}]
    assert ve.assign_series(5.0, buckets) == "early"
    assert ve.assign_series(15.0, buckets) == "peak"
    assert ve.assign_series(99.0, buckets) is None  # outside all buckets
    assert ve.assign_series(20.0, buckets) == "peak"  # inclusive LAST upper edge
    # a gap between buckets returns None even at a non-last upper edge
    gapped = [{"label": "a", "x_lo": 0.0, "x_hi": 10.0},
              {"label": "b", "x_lo": 15.0, "x_hi": 20.0}]
    assert ve.assign_series(10.0, gapped) is None


def test_extract_filled_marks_offline_synthetic(tmp_path):
    # Build a tiny PDF with 3 filled circles at known page-point centers.
    import fitz
    doc = fitz.open()
    page = doc.new_page(width=200, height=200)
    centers = [(40, 160), (100, 100), (160, 40)]  # (x_pt, y_pt), PDF origin top-left in fitz
    for cx, cy in centers:
        page.draw_circle((cx, cy), 3, color=(0, 0, 0), fill=(0, 0, 0))
    pdf = tmp_path / "synthetic.pdf"
    doc.save(pdf)
    doc.close()

    marks = ve.extract_filled_marks(
        pdf, page=0, plot_box_pt=(0, 0, 200, 200),
        min_area_pt2=1.0, max_area_pt2=200.0,
    )
    # 3 marks found, centroids near the drawn centers (order-independent)
    assert len(marks) == 3
    for cx, cy in centers:
        assert any(math.hypot(mx - cx, my - cy) < 4 for mx, my in marks)


def test_extract_panel_maps_marks_to_points(tmp_path):
    import fitz
    doc = fitz.open()
    page = doc.new_page(width=200, height=200)
    # one mark in the "early" x-bucket, one in "peak"; y maps linearly 200px->0, 0px->100
    page.draw_circle((5, 100), 3, color=(0, 0, 0), fill=(0, 0, 0))    # x=5(early), y=100px
    page.draw_circle((15, 50), 3, color=(0, 0, 0), fill=(0, 0, 0))    # x=15(peak), y=50px
    pdf = tmp_path / "s2.pdf"
    doc.save(pdf)
    doc.close()

    panel = m.Panel(
        panel_id="syn", figure_label="syn", chart_type="scatter_log", page=0,
        bbox_frac=(0, 0, 1, 1), truth_file="x", truth_sheet="x", truth_adapter="x",
        align_mode="by_nearest",
        axis_calib={
            "plot_box_pt": (0, 0, 200, 200),
            "x": {"p0_pt": 0.0, "v0": 0.0, "p1_pt": 20.0, "v1": 20.0, "log": False},
            "y": {"p0_pt": 200.0, "v0": 0.0, "p1_pt": 0.0, "v1": 100.0, "log": False},
            "mark": {"min_area_pt2": 1.0, "max_area_pt2": 200.0},
            "series_x_buckets": [
                {"label": "early", "x_lo": 0.0, "x_hi": 10.0},
                {"label": "peak", "x_lo": 10.0, "x_hi": 20.0},
            ],
        },
    )
    pts = ve.extract_panel(panel, pdf)
    g = m.group_by_series(pts)
    assert set(g) == {"early", "peak"}
    # early mark at y=100px → value 50; peak mark at y=50px → value 75
    assert g["early"][0].value == pytest.approx(50.0, abs=2.0)
    assert g["peak"][0].value == pytest.approx(75.0, abs=2.0)


def test_extract_panel_without_calib_returns_empty():
    panel = m.Panel(panel_id="nc", figure_label="x", chart_type="bar", page=0,
                    bbox_frac=None, truth_file="x", truth_sheet="x", truth_adapter="x",
                    axis_calib=None)
    assert ve.extract_panel(panel, "/nonexistent.pdf") == []
