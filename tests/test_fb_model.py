import figure_benchmark.model as m


def test_point_roundtrips_through_dict():
    p = m.Point(series_label="peak", key="pt1", value=2459.0, unit="SFC/1e6")
    d = p.to_dict()
    assert d == {"series_label": "peak", "key": "pt1", "value": 2459.0, "unit": "SFC/1e6"}
    assert m.Point.from_dict(d) == p


def test_point_allows_null_value_for_unreadable():
    p = m.Point(series_label="peak", key="pt9", value=None, unit=None)
    assert p.value is None
    assert m.Point.from_dict(p.to_dict()) == p


def test_group_by_series_buckets_points():
    pts = [
        m.Point("early", "a", 1.0, "u"),
        m.Point("peak", "a", 9.0, "u"),
        m.Point("early", "b", 2.0, "u"),
    ]
    grouped = m.group_by_series(pts)
    assert set(grouped) == {"early", "peak"}
    assert [p.key for p in grouped["early"]] == ["a", "b"]


def test_panel_lookup_fields_present():
    panel = m.Panel(
        panel_id="fig1g", figure_label="Figure 1G", chart_type="scatter_log",
        page=1, bbox_frac=(0.1, 0.2, 0.5, 0.6),
        truth_file="x.xlsx", truth_sheet="Figure 1G",
        truth_adapter="fig1g_three_block", align_mode="by_nearest",
        positivity_threshold=None, unit="SFC/1e6", notes="",
    )
    assert panel.panel_id == "fig1g"
    assert panel.align_mode == "by_nearest"


def test_panel_render_dpi_defaults_and_axis_calib_optional():
    # render_dpi defaults to 300; dense panels override higher.
    default = m.Panel(
        panel_id="x", figure_label="x", chart_type="bar", page=0, bbox_frac=None,
        truth_file="x", truth_sheet="x", truth_adapter="x",
    )
    assert default.render_dpi == 300
    assert default.axis_calib is None
    dense = m.Panel(
        panel_id="fig1g", figure_label="Figure 1G", chart_type="scatter_log",
        page=1, bbox_frac=None, truth_file="x", truth_sheet="x",
        truth_adapter="x", render_dpi=600,
        axis_calib={"y": {"p0_pt": 0.0, "v0": 1.0, "p1_pt": 100.0, "v1": 1000.0, "log": True}},
    )
    assert dense.render_dpi == 600
    assert dense.axis_calib["y"]["log"] is True
