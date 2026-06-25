"""Curated panel registry — the single source of truth for WHAT is benchmarked.

bbox_frac is None until set in Task 9 (page-render bbox discovery). The runner
skips bbox-less panels. page is the 0-based PDF page index in main.pdf.

Source-data filenames (Rojas, Nature 2023):
  MOESM5 = Figure 1 panels; MOESM6 = Figure 2 panels; MOESM7 = Figure 4 panels.
"""
from __future__ import annotations

from .model import Panel

_M5 = "41586_2023_6063_MOESM5_ESM.xlsx"
_M6 = "41586_2023_6063_MOESM6_ESM.xlsx"
_M7 = "41586_2023_6063_MOESM7_ESM.xlsx"

PANEL_REGISTRY: list[Panel] = [
    Panel(
        panel_id="fig1g", figure_label="Figure 1G", chart_type="scatter_log",
        page=2, bbox_frac=(0.806, 0.362, 0.982, 0.492), truth_file=_M5, truth_sheet="Figure 1G",
        truth_adapter="fig1g_peak", align_mode="by_nearest",
        positivity_threshold=None, unit="SFC/1e6", render_dpi=600, axis_calib=None,
        notes="HARD case / run5 gap: per-patient IFNγ ELISpot scatter (page 2 = Fig 1, "
              "0-based). The PUBLISHED panel plots only the PEAK block (truth_adapter="
              "fig1g_peak), x-axis = Patient. y-axis is SEGMENTED (linear 0-500 then "
              "compressed 500/1500/3000) — not clean log/linear, so the vector arm "
              "(axis_calib) is left None; vision arm only.",
    ),
    Panel(
        panel_id="fig4b", figure_label="Figure 4B", chart_type="bar",
        page=5, bbox_frac=(0.487, 0.200, 0.648, 0.326), truth_file=_M7, truth_sheet="Figure 4B",
        truth_adapter="fig4b_bars", align_mode="by_nearest",
        positivity_threshold=None, unit="% all blood T cells", render_dpi=400, axis_calib=None,
        notes="EASY: 3 linear bars (% all blood T cells), page 5 = Fig 4 (0-based). "
              "Bars are not circular marks, so the vector arm does not apply (axis_calib None).",
    ),
    # NOTE: 'Figure 2A bar graph' (source sheet MOESM6) was dropped from the curated
    # set — the published Fig 2 panel 'a' is per-patient line charts, and the 2x2
    # contingency bar graph the sheet contains is not locatable at that panel letter.
    # The load_fig2a_bar adapter remains in truth_loader for if/when it is relocated.
]

_BY_ID = {p.panel_id: p for p in PANEL_REGISTRY}


def get(panel_id: str) -> Panel:
    return _BY_ID[panel_id]
