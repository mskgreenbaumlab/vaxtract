"""Numeric-readout vision reader.

Builds a numeric-extraction prompt (distinct from the categorical SKILL.md in
docs/bootstrap/figure_vision), calls a Provider (default: the bootstrap
ClaudeCodeProvider), parses the response into list[Point], and caches by
(image bytes, prompt, provider name).
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import re
import time
from dataclasses import dataclass

from .model import Point

_FENCE_RE = re.compile(r"```(?:json)?\s*(.*?)\s*```", re.DOTALL)


@dataclass
class ReadResult:
    points: list[Point]
    raw_response: str
    cost_usd: float | None = None
    wall_time_s: float = 0.0
    from_cache: bool = False
    error: str | None = None


def build_numeric_prompt(panel, image_path) -> str:
    return "\n".join([
        f"You are reading numeric data points off a single chart panel: {panel.figure_label}.",
        f"The panel image is at the absolute path: {image_path}",
        f"Chart type: {panel.chart_type}. Expected value unit: {panel.unit!r}.",
        "",
        "Read the data series and their numeric values directly off the plotted marks.",
        "RULES:",
        "- Inspect the y-axis scale FIRST. If it is logarithmic, read values on the log scale "
        "(do not linearly interpolate pixel positions).",
        "- Return the actual plotted value for every visible data point.",
        "- If a point is present but you cannot read its value, set value to null. NEVER guess "
        "or fabricate a number.",
        "- Do not invent points that are not plotted.",
        "",
        "Return ONLY JSON in this exact shape (no prose):",
        '{"series": [{"series_label": "<name>", "unit": "<unit>", '
        '"points": [{"key": "<x or category>", "value": <number|null>}]}]}',
    ])


def _extract_json(text: str) -> str:
    m = _FENCE_RE.search(text)
    return (m.group(1) if m else text).strip()


def _coerce_value(v) -> float | None:
    """Coerce a JSON value to float | None without fabricating.

    bool is excluded (it is an int subclass, so a JSON true/false must NOT become
    1.0/0.0 — that would turn a clearly wrong model output into a plausible
    number). Numeric strings ("1.7") are parsed; non-numeric strings and anything
    else become None (present-but-unreadable), never a guessed value.
    """
    if isinstance(v, bool):
        return None
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, str):
        try:
            return float(v.strip())
        except ValueError:
            return None
    return None


def parse_numeric_response(text: str) -> list[Point]:
    try:
        parsed = json.loads(_extract_json(text))
    except json.JSONDecodeError as e:
        raise ValueError(f"could not parse vision response as JSON: {e}; first 200: {text[:200]!r}")
    if not isinstance(parsed, dict) or "series" not in parsed:
        raise ValueError(f"expected object with 'series'; got keys={list(parsed) if isinstance(parsed, dict) else type(parsed)}")
    out: list[Point] = []
    for s in parsed["series"]:
        label = str(s.get("series_label", ""))
        unit = s.get("unit")
        for p in s.get("points", []):
            v = _coerce_value(p.get("value"))
            out.append(Point(series_label=label, key=str(p.get("key", "")), value=v, unit=unit))
    return out


def _cache_key(image_path, prompt, provider_name) -> str:
    h = hashlib.sha256()
    h.update(pathlib.Path(image_path).read_bytes())
    h.update(b"\x00")
    h.update(prompt.encode())
    h.update(b"\x00")
    h.update(provider_name.encode())
    return h.hexdigest()


def read_panel(panel, image_path, provider, *, cache_dir=None, effort="medium") -> ReadResult:
    prompt = build_numeric_prompt(panel, image_path)
    cpath = None
    if cache_dir is not None:
        cache_dir = pathlib.Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        cpath = cache_dir / f"{_cache_key(image_path, prompt, provider.name)}.json"
        if cpath.exists():
            d = json.loads(cpath.read_text())
            return ReadResult(
                points=[Point.from_dict(x) for x in d["points"]],
                raw_response=d["raw_response"], cost_usd=d.get("cost_usd"),
                wall_time_s=d.get("wall_time_s", 0.0), from_cache=True,
            )
    t0 = time.monotonic()
    raw, meta = provider.extract(image_path=pathlib.Path(image_path), prompt=prompt, effort=effort)
    wall = meta.get("wall_time_s", time.monotonic() - t0)
    try:
        points = parse_numeric_response(raw)
        error = None
    except ValueError as e:
        points, error = [], str(e)
    result = ReadResult(points=points, raw_response=raw, cost_usd=meta.get("cost_usd"),
                        wall_time_s=wall, error=error)
    if cpath is not None and error is None:
        cpath.write_text(json.dumps({
            "points": [p.to_dict() for p in points], "raw_response": raw,
            "cost_usd": result.cost_usd, "wall_time_s": wall,
        }, indent=2))
    return result
