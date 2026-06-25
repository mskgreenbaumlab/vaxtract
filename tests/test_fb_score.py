import math

import figure_benchmark.model as m
import figure_benchmark.score_numeric as sc


def _pts(series, pairs, unit="u"):
    return [m.Point(series, k, v, unit) for k, v in pairs]


def test_perfect_match_zero_log_error():
    truth = _pts("s", [("a", 10.0), ("b", 100.0), ("c", 1000.0)])
    pred = _pts("s", [("a", 10.0), ("b", 100.0), ("c", 1000.0)])
    r = sc.score(pred, truth, align_mode="by_key")
    assert r["per_point"]["median_abs_log10"] == 0.0
    assert r["per_point"]["n_matched"] == 3
    assert r["per_point"]["n_missed"] == 0
    assert r["per_point"]["n_extra"] == 0


def test_known_offset_known_log_error():
    truth = _pts("s", [("a", 100.0)])
    pred = _pts("s", [("a", 1000.0)])  # 10x high → log10 error = 1.0
    r = sc.score(pred, truth, align_mode="by_key")
    assert r["per_point"]["median_abs_log10"] == 1.0


def test_missing_and_extra_counted():
    truth = _pts("s", [("a", 10.0), ("b", 20.0)])
    pred = _pts("s", [("a", 10.0), ("z", 99.0)])  # b missed, z extra
    r = sc.score(pred, truth, align_mode="by_key")
    assert r["per_point"]["n_matched"] == 1
    assert r["per_point"]["n_missed"] == 1
    assert r["per_point"]["n_extra"] == 1


def test_per_series_peak_and_mean_log_rmse():
    truth = _pts("s", [("a", 10.0), ("b", 1000.0)])  # peak 1000, mean 505
    pred = _pts("s", [("a", 10.0), ("b", 1000.0)])
    r = sc.score(pred, truth, align_mode="by_key")
    assert r["per_series"]["peak_log_rmse"] == 0.0
    assert r["per_series"]["mean_log_rmse"] == 0.0


def test_categorical_responder_calls_with_threshold():
    # threshold 50: truth a=neg(10) b=pos(900); pred a=neg(20) b=pos(800)
    truth = _pts("s", [("a", 10.0), ("b", 900.0)])
    pred = _pts("s", [("a", 20.0), ("b", 800.0)])
    r = sc.score(pred, truth, align_mode="by_key", positivity_threshold=50.0)
    assert r["categorical"]["n_positive_truth"] == 1
    assert r["categorical"]["n_positive_pred"] == 1
    assert r["categorical"]["call_agreement"] == 1.0  # both points agree


def test_nearest_alignment_for_unordered_scatter():
    truth = _pts("s", [("p1", 5.0), ("p2", 500.0)])
    pred = _pts("s", [("x", 480.0), ("y", 6.0)])  # unordered; nearest pairs 500<->480, 5<->6
    r = sc.score(pred, truth, align_mode="by_nearest")
    assert r["per_point"]["n_matched"] == 2
    assert r["per_point"]["median_abs_log10"] < 0.1


def test_zero_values_routed_not_logged():
    truth = _pts("s", [("a", 0.0), ("b", 100.0)])
    pred = _pts("s", [("a", 0.0), ("b", 100.0)])
    r = sc.score(pred, truth, align_mode="by_key", positivity_threshold=50.0)
    # the zero pair is excluded from log stats but matched + categorical-correct
    assert r["per_point"]["n_log_pairs"] == 1
    assert r["categorical"]["call_agreement"] == 1.0


def test_collapse_series_matches_across_mismatched_labels():
    # vision invents its own series_label; collapse must let value-matching work
    truth = _pts("peak", [("p1", 0.7), ("p2", 1.7), ("p3", 2.5)])
    pred = _pts("blood_T_pct", [("x", 0.7), ("y", 1.6), ("z", 2.5)])  # different label
    # without collapse, the label mismatch zeroes out all matches
    bad = sc.score(pred, truth, align_mode="by_nearest")
    assert bad["per_point"]["n_matched"] == 0
    # with collapse, all three align by value
    good = sc.score(pred, truth, align_mode="by_nearest", collapse_series=True)
    assert good["per_point"]["n_matched"] == 3
    assert good["per_point"]["median_abs_log10"] < 0.05


def test_by_key_duplicate_keys_not_silently_dropped():
    truth = _pts("s", [("a", 10.0)])
    pred = _pts("s", [("a", 10.0), ("a", 11.0)])  # duplicate pred key
    r = sc.score(pred, truth, align_mode="by_key")
    assert r["per_point"]["n_matched"] == 1
    assert r["per_point"]["n_extra"] == 1  # the duplicate is extra, not dropped


def test_by_nearest_null_pred_is_not_a_match_candidate():
    truth = _pts("s", [("p1", 5.0), ("p2", 500.0)])
    pred = [m.Point("s", "x", None, "u"), m.Point("s", "y", 6.0, "u")]
    r = sc.score(pred, truth, align_mode="by_nearest")
    # only the real value 6 matches (to 5); 500 missed; the None pred is extra
    assert r["per_point"]["n_matched"] == 1
    assert r["per_point"]["n_missed"] == 1
    assert r["per_point"]["n_extra"] == 1
