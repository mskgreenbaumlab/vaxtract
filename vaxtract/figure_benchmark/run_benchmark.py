"""Benchmark runner: registry → render → truth → [vision read + vector read] → score → report.

The agent's read_figure stub is NOT touched. Vision/vector numbers never enter a
real record — this writes only report.json / report.html under out_dir. The
vector arm (Task 7) runs only when a panel has axis_calib; it is an experimental
comparison ceiling, never a gate.
"""
from __future__ import annotations

import argparse
import html
import json
import pathlib
import sys

from . import panels as panels_mod
from . import render_panel, score_numeric, truth_loader, vector_extract
from .numeric_reader import read_panel


def run_benchmark(panels, *, provider, pdf_path, source_dir, out_dir, dpi=None,
                  cache_dir=None) -> list[dict]:
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fixtures = out_dir / "fixtures"
    cache_dir = pathlib.Path(cache_dir) if cache_dir else out_dir / ".cache"
    results: list[dict] = []
    for panel in panels:
        if panel.bbox_frac is None:
            results.append({"panel_id": panel.panel_id, "status": "skipped_no_bbox"})
            continue
        # dpi=None → render_panel falls back to panel.render_dpi (Fig 1G → 600)
        img = render_panel.render_panel(panel, pdf_path, out_dir=fixtures, dpi=dpi)
        truth = truth_loader.load_truth(panel, source_dir=source_dir)

        # The vision model invents its own series labels/keys; when the truth is a
        # single series we collapse both sides so matching is by value, not by a
        # label the model could never reproduce.
        collapse = len({t.series_label for t in truth}) <= 1

        # --- vision arm ---
        read = read_panel(panel, img, provider, cache_dir=cache_dir)
        metrics = score_numeric.score(
            read.points, truth, align_mode=panel.align_mode,
            positivity_threshold=panel.positivity_threshold, collapse_series=collapse,
        )

        # --- vector arm (optional comparison ceiling) ---
        vec_points = vector_extract.extract_panel(panel, pdf_path)  # [] if no axis_calib
        metrics_vector = None
        if vec_points:
            metrics_vector = score_numeric.score(
                vec_points, truth, align_mode=panel.align_mode,
                positivity_threshold=panel.positivity_threshold, collapse_series=collapse,
            )

        results.append({
            "panel_id": panel.panel_id,
            "figure_label": panel.figure_label,
            "chart_type": panel.chart_type,
            "status": "scored",
            "image": str(img),
            "metrics": metrics,
            "metrics_vector": metrics_vector,
            "read_error": read.error,
            "cost_usd": read.cost_usd,
            "wall_time_s": read.wall_time_s,
            "from_cache": read.from_cache,
            "truth": [p.to_dict() for p in truth],
            "pred": [p.to_dict() for p in read.points],
            "pred_vector": [p.to_dict() for p in vec_points],
        })
    (out_dir / "report.json").write_text(json.dumps(results, indent=2))
    (out_dir / "report.html").write_text(_render_html(results))
    return results


def _render_html(results) -> str:
    cards = []
    for r in results:
        if r.get("status") != "scored":
            cards.append(f"<section><h2>{html.escape(r['panel_id'])}</h2>"
                         f"<p><em>{html.escape(r.get('status',''))}</em></p></section>")
            continue
        cards.append(
            "<section style='border:1px solid #ccc;margin:1em 0;padding:1em'>"
            f"<h2>{html.escape(r['figure_label'])} "
            f"<small>({html.escape(r['chart_type'])})</small></h2>"
            f"<img src='{html.escape(r['image'])}' style='max-width:520px;border:1px solid #eee'><br>"
            + _arm_html("vision", r["metrics"])
            + (_arm_html("vector (experimental ceiling)", r["metrics_vector"])
               if r.get("metrics_vector") else "")
            + f"<b>cost</b>: {_fmt_cost(r['cost_usd'])} "
            f"wall={r['wall_time_s']:.1f}s cache={r['from_cache']}"
            "</section>"
        )
    return ("<!doctype html><meta charset='utf-8'>"
            "<title>figure_benchmark report</title>"
            "<h1>Figure-vision numeric-readout benchmark</h1>" + "".join(cards))


def _fmt_cost(cost) -> str:
    # cost_usd may be None when the provider reports no telemetry (e.g. cache hit)
    return f"${cost:.4f}" if isinstance(cost, (int, float)) else "$N/A"


def _arm_html(label: str, metrics: dict) -> str:
    pp = metrics["per_point"]
    ps = metrics["per_series"]
    cat = metrics["categorical"]
    return (
        f"<div style='margin:.4em 0;padding:.4em;background:#f7f7f7'>"
        f"<b>{html.escape(label)}</b><br>"
        f"per-point: matched={pp['n_matched']} missed={pp['n_missed']} "
        f"extra={pp['n_extra']} | median|log10|={pp['median_abs_log10']:.3f} "
        f"p90={pp['p90_abs_log10']:.3f} (n_log_pairs={pp['n_log_pairs']})<br>"
        f"per-series: peak_log_rmse={ps['peak_log_rmse']:.3f} "
        f"mean_log_rmse={ps['mean_log_rmse']:.3f}<br>"
        f"categorical: pos_truth={cat['n_positive_truth']} "
        f"pos_pred={cat['n_positive_pred']} agreement={cat['call_agreement']}"
        f"</div>"
    )


# Reuse the bootstrap ClaudeCodeProvider TRANSPORT but NOT its categorical
# SKILL.md (that would instruct evidence-row output and fight our numeric ask).
# All numeric instructions live in build_numeric_prompt; the system prompt is a
# one-line framing only.
_NUMERIC_SYSTEM = (
    "You are a precise chart-reading assistant. You read numeric values off a "
    "single plotted panel and return them as JSON exactly as the user instructs. "
    "You never fabricate values; unreadable points are null."
)


def _default_provider():
    """CLI-only: import the bootstrap ClaudeCodeProvider TRANSPORT lazily."""
    boot = pathlib.Path("docs/bootstrap")
    if str(boot) not in sys.path:
        sys.path.insert(0, str(boot))
    from figure_vision.worker import ClaudeCodeProvider  # noqa: E402
    # Construct directly (NOT from_skill_md_file) so the categorical SKILL.md is
    # not appended as the system prompt.
    return ClaudeCodeProvider(skill_md=_NUMERIC_SYSTEM, model="claude-opus-4-8")


def main(argv=None):
    ap = argparse.ArgumentParser(description="figure_benchmark runner")
    ap.add_argument("--pdf", default="agentBranch/cancerVacExtrac_claudeAi/data/raw/rojas_37165196/main.pdf")
    ap.add_argument("--source-dir", default="agentBranch/cancerVacExtrac_claudeAi/data/raw/rojas_37165196/sourceData")
    ap.add_argument("--out-dir", default="agentBranch/cancerVacExtrac_claudeAi/outputs/figure_benchmark")
    ap.add_argument("--dpi", type=int, default=None,
                    help="override render DPI for all panels; omit to use each panel's render_dpi")
    ap.add_argument("--only", nargs="*", help="restrict to these panel_ids")
    args = ap.parse_args(argv)

    sel = panels_mod.PANEL_REGISTRY
    if args.only:
        sel = [p for p in sel if p.panel_id in set(args.only)]
    results = run_benchmark(
        sel, provider=_default_provider(), pdf_path=args.pdf,
        source_dir=args.source_dir, out_dir=args.out_dir, dpi=args.dpi,
    )
    scored = [r for r in results if r.get("status") == "scored"]
    print(f"scored {len(scored)}/{len(results)} panels → {args.out_dir}/report.html")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
