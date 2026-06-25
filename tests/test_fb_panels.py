import figure_benchmark.panels as P
import figure_benchmark.truth_loader as tl


def test_registry_has_fig1g_locked_as_scatter_log():
    g = P.get("fig1g")
    assert g.figure_label == "Figure 1G"
    assert g.chart_type == "scatter_log"
    assert g.align_mode == "by_nearest"
    assert g.truth_adapter == "fig1g_peak"  # published panel plots the peak block only
    assert g.bbox_frac is not None  # bbox located during discovery
    assert g.render_dpi == 600  # dense crop renders at higher DPI


def test_every_panel_adapter_is_registered():
    for panel in P.PANEL_REGISTRY:
        assert panel.truth_adapter in tl.ADAPTERS, panel.panel_id


def test_panel_ids_unique():
    ids = [p.panel_id for p in P.PANEL_REGISTRY]
    assert len(ids) == len(set(ids))
