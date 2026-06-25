import pathlib
import pytest

import figure_benchmark.model as m
import figure_benchmark.truth_loader as tl

SRC = pathlib.Path(
    "agentBranch/cancerVacExtrac_claudeAi/data/raw/rojas_37165196/sourceData"
)
FIG1G = SRC / "41586_2023_6063_MOESM5_ESM.xlsx"
FIG2A = SRC / "41586_2023_6063_MOESM6_ESM.xlsx"
FIG4B = SRC / "41586_2023_6063_MOESM7_ESM.xlsx"

pytestmark = pytest.mark.skipif(
    not FIG1G.exists(), reason="Rojas source data symlink not present"
)


def test_fig1g_three_block_anchors_and_invariants():
    pts = tl.load_fig1g_three_block(FIG1G, "Figure 1G")
    g = m.group_by_series(pts)
    # three timepoint series
    assert set(g) == {"early", "peak", "late"}
    # patient 1's first row: early=2, peak=2459, late=757 (from row2 of the sheet)
    early1 = g["early"][0]
    assert early1.key.startswith("1#")
    assert early1.value == 2.0
    assert g["peak"][0].value == 2459.0
    assert g["late"][0].value == 757.0
    # peak series max is the patient-1 priming spike (2459); all values >= 0
    peak_vals = [p.value for p in g["peak"] if p.value is not None]
    assert max(peak_vals) == 2459.0
    assert min(p.value for p in pts if p.value is not None) >= 0.0
    # forward-fill: every early point has a non-empty patient prefix
    assert all("#" in p.key and p.key.split("#")[0] for p in g["early"])
    # unit carried through
    assert early1.unit and "IFN" in early1.unit or early1.unit == "SFC/1e6"


def test_fig1g_peak_only_one_series_matches_plotted_panel():
    pts = tl.load_fig1g_peak(FIG1G, "Figure 1G")
    g = m.group_by_series(pts)
    assert set(g) == {"peak"}  # the published panel plots only the peak block
    by_patient = {}
    for p in pts:
        by_patient.setdefault(p.key.split("#")[0], []).append(p.value)
    # anchors verified against the plotted dots
    assert by_patient["1"] == [2459.0]   # patient 1: single peak dot
    assert 484.0 in by_patient["5"]
    assert 442.0 in by_patient["14"]
    # peak-only is a strict subset of the full 3-block load
    assert len(pts) < len(tl.load_fig1g_three_block(FIG1G, "Figure 1G"))


def test_fig2a_bar_contingency():
    pts = tl.load_fig2a_bar(FIG2A, "Figure 2A bar graph")
    # 2x2 contingency: ELISPOT+/- x TCRVb+/- as percentages
    by = {(p.series_label, p.key): p.value for p in pts}
    assert by[("TCRVβ +", "ELISPOT +")] == 100.0
    assert by[("TCRVβ -", "ELISPOT -")] == 88.0


def test_fig4b_three_bars():
    pts = tl.load_fig4b_bars(FIG4B, "Figure 4B")
    vals = {p.key: p.value for p in pts}
    assert vals["Blood (30 weeks)"] == pytest.approx(1.663, abs=1e-3)
    assert vals["Liver mass (35 weeks)"] == pytest.approx(2.476, abs=1e-3)
    assert len(pts) == 3


def test_dispatch_by_adapter_name():
    panel = m.Panel(
        panel_id="fig4b", figure_label="Figure 4B", chart_type="bar",
        page=0, bbox_frac=None, truth_file=FIG4B.name, truth_sheet="Figure 4B",
        truth_adapter="fig4b_bars", unit="% all blood T cells",
    )
    pts = tl.load_truth(panel, source_dir=SRC)
    assert len(pts) == 3
