"""Shared data model for figure_benchmark.

Both ground-truth (truth_loader) and prediction (numeric_reader) are expressed
as ``list[Point]`` so the scorer is symmetric. ``value=None`` means "present but
unreadable" (vision) or "blank cell" (truth) — never fabricated.
"""
from __future__ import annotations

from dataclasses import dataclass
from collections import OrderedDict


@dataclass(frozen=True)
class Point:
    series_label: str
    key: str
    value: float | None
    unit: str | None = None

    def to_dict(self) -> dict:
        return {
            "series_label": self.series_label,
            "key": self.key,
            "value": self.value,
            "unit": self.unit,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "Point":
        return cls(
            series_label=d["series_label"],
            key=d["key"],
            value=d["value"],
            unit=d.get("unit"),
        )


@dataclass(frozen=True)
class Panel:
    panel_id: str
    figure_label: str
    chart_type: str                 # 'bar' | 'line' | 'scatter_log' | 'heatmap'
    page: int                       # 0-based PDF page index
    bbox_frac: tuple | None         # (x0,y0,x1,y1) as fractions of page; None until set
    truth_file: str                 # filename under data/raw/.../sourceData/
    truth_sheet: str
    truth_adapter: str              # key into truth_loader.ADAPTERS
    align_mode: str = "by_key"      # 'by_key' | 'by_nearest'
    positivity_threshold: float | None = None
    unit: str | None = None
    render_dpi: int = 300           # per-panel raster DPI; dense panels override higher
    axis_calib: dict | None = None  # vector-arm calibration (see vector_extract); None = vector arm skipped
    notes: str = ""


def group_by_series(points: list[Point]) -> "OrderedDict[str, list[Point]]":
    out: "OrderedDict[str, list[Point]]" = OrderedDict()
    for p in points:
        out.setdefault(p.series_label, []).append(p)
    return out
