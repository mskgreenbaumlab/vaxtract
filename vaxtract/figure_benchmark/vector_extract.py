"""Vector-extraction COMPARISON arm (experimental).

Reads plotted marker coordinates from the PDF's vector drawing operators and
maps them to data values via a per-panel axis calibration. This is a "ceiling"
reference for the vision arm on vector figures — NOT a graded gate, and NOT a
general capability (it cannot read raster figures, which is exactly why the
vision arm is the thing being validated). Runs only when Panel.axis_calib is set.
"""
from __future__ import annotations

import math

import fitz

from .model import Panel, Point


def interp_axis(pt: float, axis: dict) -> float:
    """Map a coordinate (PDF points) to a data value via two anchors.

    Linear interpolation in value space, or in log10 space when axis['log'].
    """
    p0, v0, p1, v1 = axis["p0_pt"], axis["v0"], axis["p1_pt"], axis["v1"]
    frac = (pt - p0) / (p1 - p0) if p1 != p0 else 0.0
    if axis.get("log"):
        if v0 <= 0 or v1 <= 0:
            raise ValueError(
                f"log axis requires positive anchors; got v0={v0}, v1={v1}"
            )
        lv = math.log10(v0) + frac * (math.log10(v1) - math.log10(v0))
        return 10 ** lv
    return v0 + frac * (v1 - v0)


def assign_series(x_value: float, buckets: list[dict]) -> str | None:
    for b in buckets:
        if b["x_lo"] <= x_value < b["x_hi"]:
            return b["label"]
    # inclusive upper edge for the LAST bucket only (so a mark exactly on the
    # overall right edge is kept). Interior boundaries are already caught by the
    # half-open pass above; gaps between buckets correctly return None.
    if buckets and x_value == buckets[-1]["x_hi"]:
        return buckets[-1]["label"]
    return None


def _centroid_and_area(rect: "fitz.Rect") -> tuple[float, float, float]:
    cx = (rect.x0 + rect.x1) / 2.0
    cy = (rect.y0 + rect.y1) / 2.0
    area = abs(rect.width * rect.height)
    return cx, cy, area


def extract_filled_marks(pdf_path, page: int, plot_box_pt, min_area_pt2, max_area_pt2):
    """Centroids (x_pt, y_pt) of small filled drawings inside the plot box.

    Heuristic: a data marker is a *filled* path whose bbox area is within
    [min_area_pt2, max_area_pt2] and whose centroid is inside plot_box_pt.
    """
    x0, y0, x1, y1 = plot_box_pt
    doc = fitz.open(pdf_path)
    try:
        pg = doc[page]
        out = []
        for d in pg.get_drawings():
            if not d.get("fill"):
                continue
            rect = d.get("rect")
            if rect is None:
                continue
            cx, cy, area = _centroid_and_area(rect)
            if not (min_area_pt2 <= area <= max_area_pt2):
                continue
            if not (x0 <= cx <= x1 and y0 <= cy <= y1):
                continue
            out.append((cx, cy))
        return out
    finally:
        doc.close()


def extract_panel(panel: Panel, pdf_path) -> list[Point]:
    calib = panel.axis_calib
    if not calib:
        return []
    marks = extract_filled_marks(
        pdf_path, panel.page, calib["plot_box_pt"],
        calib["mark"]["min_area_pt2"], calib["mark"]["max_area_pt2"],
    )
    pts: list[Point] = []
    for i, (cx, cy) in enumerate(marks):
        x_val = interp_axis(cx, calib["x"])
        y_val = interp_axis(cy, calib["y"])
        series = assign_series(x_val, calib["series_x_buckets"])
        if series is None:
            continue
        pts.append(Point(series_label=series, key=f"vec#{i}", value=y_val, unit=panel.unit))
    return pts
