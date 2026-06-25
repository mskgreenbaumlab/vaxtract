"""Symmetric scorer: prediction list[Point] vs truth list[Point].

Three metric blocks, all reported (no committed pass threshold):
  per_point   — |log10(pred/truth)| distribution over matched, both-positive pairs;
                match/miss/extra counts.
  per_series  — log-RMSE of each series' peak and mean summary (what the DB stores).
  categorical — responder/non calls vs a positivity threshold; n-positive; agreement.

Alignment:
  by_key     — match points sharing the same key within a series.
  by_nearest — greedy nearest-value matching within a series (for unordered scatter).
"""
from __future__ import annotations

import math
import statistics

from .model import Point, group_by_series


def _abs_log10_ratio(a: float, b: float) -> float:
    return abs(math.log10(a) - math.log10(b))


def _first_by_key(points):
    """Index points by key, keeping the FIRST per key; return (index, duplicates).

    Duplicate keys are never silently dropped — they are returned so the caller
    can route them to missed/extra (a duplicated row is a grading signal, not a
    no-op)."""
    by, dups = {}, []
    for p in points:
        if p.key in by:
            dups.append(p)
        else:
            by[p.key] = p
    return by, dups


def _align_by_key(pred, truth):
    pred_by, pred_dups = _first_by_key(pred)
    truth_by, truth_dups = _first_by_key(truth)
    matched = [(pred_by[k], truth_by[k]) for k in truth_by if k in pred_by]
    missed = [t for k, t in truth_by.items() if k not in pred_by] + truth_dups
    extra = [p for k, p in pred_by.items() if k not in truth_by] + pred_dups
    return matched, missed, extra


def _align_by_nearest(pred, truth):
    # value=None predictions are NOT match candidates — an unreadable point must
    # not consume a truth slot (that would inflate the match rate). They go to extra.
    remaining = [p for p in pred if p.value is not None]
    extra_null = [p for p in pred if p.value is None]
    matched, missed = [], []
    for t in sorted(truth, key=lambda p: (p.value is None, p.value or 0.0)):
        if t.value is None or not remaining:
            missed.append(t)
            continue
        best = min(remaining, key=lambda p: abs(p.value - t.value))
        matched.append((best, t))
        remaining.remove(best)
    extra = remaining + extra_null
    return matched, missed, extra


def _align(pred, truth, align_mode):
    pg, tg = group_by_series(pred), group_by_series(truth)
    matched, missed, extra = [], [], []
    for label in tg:
        p_series = pg.get(label, [])
        if align_mode == "by_nearest":
            mm, ms, ex = _align_by_nearest(p_series, tg[label])
        else:
            mm, ms, ex = _align_by_key(p_series, tg[label])
        matched += mm
        missed += ms
        extra += ex
    # series present in pred but not truth → all extra
    for label in pg:
        if label not in tg:
            extra += pg[label]
    return matched, missed, extra


def _per_point(matched, missed, extra):
    # NB: p90 is a rank-based percentile. For n_log_pairs < 10 it collapses to the
    # max (ceil(0.9*n)==n), so on tiny panels read p90 as "worst point", not a tail.
    logs = [
        _abs_log10_ratio(p.value, t.value)
        for p, t in matched
        if p.value is not None and t.value is not None and p.value > 0 and t.value > 0
    ]
    return {
        "n_matched": len(matched),
        "n_missed": len(missed),
        "n_extra": len(extra),
        "n_log_pairs": len(logs),
        "median_abs_log10": statistics.median(logs) if logs else 0.0,
        "p90_abs_log10": (sorted(logs)[max(0, math.ceil(0.9 * len(logs)) - 1)] if logs else 0.0),
    }


def _summ_log_rmse(pred_vals, truth_vals, reducer):
    errs = []
    for label in truth_vals:
        if label not in pred_vals or not pred_vals[label] or not truth_vals[label]:
            continue
        pv, tv = reducer(pred_vals[label]), reducer(truth_vals[label])
        if pv > 0 and tv > 0:
            errs.append((math.log10(pv) - math.log10(tv)) ** 2)
    return math.sqrt(sum(errs) / len(errs)) if errs else 0.0


def _per_series(pred, truth):
    pg, tg = group_by_series(pred), group_by_series(truth)
    pv = {k: [p.value for p in v if p.value is not None and p.value > 0] for k, v in pg.items()}
    tv = {k: [p.value for p in v if p.value is not None and p.value > 0] for k, v in tg.items()}
    return {
        "peak_log_rmse": _summ_log_rmse(pv, tv, max),
        "mean_log_rmse": _summ_log_rmse(pv, tv, lambda xs: sum(xs) / len(xs)),
    }


def _categorical(matched, threshold):
    if threshold is None:
        return {"n_positive_truth": None, "n_positive_pred": None, "call_agreement": None}
    n_pos_t = n_pos_p = agree = total = 0
    for p, t in matched:
        if p.value is None or t.value is None:
            continue
        tp, pp = t.value >= threshold, p.value >= threshold
        n_pos_t += int(tp)
        n_pos_p += int(pp)
        agree += int(tp == pp)
        total += 1
    return {
        "n_positive_truth": n_pos_t,
        "n_positive_pred": n_pos_p,
        "call_agreement": (agree / total) if total else None,
    }


def _collapse_to_one_series(points):
    # The vision model invents its own series_label/key strings, which will not
    # match the truth adapter's. For single-series panels we pool every point into
    # one bucket so alignment is by value (by_nearest) rather than by a label the
    # model could never reproduce.
    return [Point("_all", p.key, p.value, p.unit) for p in points]


def score(pred, truth, *, align_mode="by_key", positivity_threshold=None,
          collapse_series=False) -> dict:
    if collapse_series:
        pred = _collapse_to_one_series(pred)
        truth = _collapse_to_one_series(truth)
    matched, missed, extra = _align(pred, truth, align_mode)
    return {
        "per_point": _per_point(matched, missed, extra),
        "per_series": _per_series(pred, truth),
        "categorical": _categorical(matched, positivity_threshold),
        "n_truth": len(truth),
        "n_pred": len(pred),
    }
