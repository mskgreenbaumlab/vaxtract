"""Per-panel ground-truth adapters: Source Data xlsx → list[Point].

Each adapter knows the idiosyncratic layout of one panel's sheet. Adapters are
registered in ADAPTERS and dispatched by Panel.truth_adapter via load_truth().
"""
from __future__ import annotations

import pathlib

import openpyxl

from .model import Panel, Point

_SFC_UNIT = "SFC/1e6"


def _rows(path: pathlib.Path, sheet: str) -> list[list]:
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        ws = wb[sheet]
        return [list(r) for r in ws.iter_rows(values_only=True)]
    finally:
        wb.close()


def _num(v) -> float | None:
    if isinstance(v, (int, float)):
        return float(v)
    return None


def load_fig1g_three_block(path, sheet="Figure 1G") -> list[Point]:
    """Three side-by-side timepoint blocks; col A patient is forward-filled.

    Cols (0-based): A=0 patient; early B=1/C=2; peak E=4/F=5; late H=7/I=8.
    Data rows are sheet rows 2.. (after the title row 0 and header row 1).
    """
    rows = _rows(pathlib.Path(path), sheet)
    blocks = [("early", 2), ("peak", 5), ("late", 8)]
    out: list[Point] = []
    patient = None
    for idx, row in enumerate(rows[2:], start=2):
        a = row[0] if len(row) > 0 else None
        if isinstance(a, (int, float, str)) and str(a).strip():
            patient = str(int(a)) if isinstance(a, float) and a.is_integer() else str(a).strip()
        if patient is None:
            continue
        key = f"{patient}#{idx}"
        for label, col in blocks:
            val = _num(row[col]) if len(row) > col else None
            # emit only when the SFC cell is a real number (blank tails differ per block)
            if val is not None:
                out.append(Point(series_label=label, key=key, value=val, unit=_SFC_UNIT))
    return out


def load_fig1g_peak(path, sheet="Figure 1G") -> list[Point]:
    """The PUBLISHED Fig 1G panel plots only the PEAK timepoint block (cols E/F) —
    one dot per patient/pool, x-axis = Patient. The sheet also carries early/late
    blocks that are NOT drawn in the panel (use load_fig1g_three_block for all
    three). Verified by anchors: patient 1 = 2459, patient 5 = 484, patient 14 =
    442, patient 6 = 85 — all peak-block singletons matching the plotted dots."""
    return [p for p in load_fig1g_three_block(path, sheet) if p.series_label == "peak"]


def load_fig2a_bar(path, sheet="Figure 2A bar graph") -> list[Point]:
    """2x2 contingency percentages. Row1 = column headers (TCRVβ +/-),
    col0 of rows 2.. = row headers (ELISPOT +/-)."""
    rows = _rows(pathlib.Path(path), sheet)
    col_headers = [c for c in rows[1][1:] if c]  # ['TCRVβ +', 'TCRVβ -']
    out: list[Point] = []
    for row in rows[2:]:
        if not row or not row[0]:
            continue
        rh = str(row[0]).strip()
        for j, ch in enumerate(col_headers, start=1):
            val = _num(row[j]) if len(row) > j else None
            if val is not None:
                out.append(Point(series_label=str(ch).strip(), key=rh, value=val, unit="%"))
    return out


def load_fig4b_bars(path, sheet="Figure 4B") -> list[Point]:
    """Single-series bars: col0 = sample label, col1 = % all blood T cells.
    Data starts after the title row and a header row."""
    rows = _rows(pathlib.Path(path), sheet)
    out: list[Point] = []
    for row in rows:
        if not row or not row[0]:
            continue
        label = str(row[0]).strip()
        val = _num(row[1]) if len(row) > 1 else None
        if val is not None and label.lower() != "sample":
            out.append(Point(series_label="blood_T_pct", key=label, value=val, unit="%"))
    return out


ADAPTERS = {
    "fig1g_three_block": load_fig1g_three_block,
    "fig1g_peak": load_fig1g_peak,
    "fig2a_bar": load_fig2a_bar,
    "fig4b_bars": load_fig4b_bars,
}


def load_truth(panel: Panel, source_dir) -> list[Point]:
    fn = ADAPTERS[panel.truth_adapter]
    path = pathlib.Path(source_dir) / panel.truth_file
    return fn(path, panel.truth_sheet)
