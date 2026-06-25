"""Unit tests for the HPC batch wrapper's pure logic (scale/extract_one.py).

Covers the exit-code/marker contract that makes the subscription batch resumable + quota-safe:
parse_summary (ledger fields from the runner line), usage_limit_hit (retry-later signal), and
classify (the ok/invalid/usage_limit/error -> retryable decision). No network/subprocess."""
import pathlib
import sys

PKT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PKT / "scale"))
import extract_one as E  # noqa: E402


# ---- parse_summary ----

def test_parse_summary_success_line():
    s = E.parse_summary("[info] agent loop ok: subtype=success turns=27 cost=$4.08 errors=None")
    assert s == {"subtype": "success", "turns": 27, "cost_usd": 4.08}

def test_parse_summary_error_line():
    s = E.parse_summary("[warn] agent loop ENDED IN ERROR: subtype=error_max_turns turns=120 cost=$13.34 errors=None")
    assert s["subtype"] == "error_max_turns" and s["turns"] == 120 and s["cost_usd"] == 13.34

def test_parse_summary_absent():
    assert E.parse_summary("nothing useful here") == {}


# ---- usage_limit_hit ----

def test_usage_limit_patterns():
    for t in ["Claude usage limit reached", "HTTP 429 Too Many Requests",
              "rate limit exceeded", "your limit resets at 5pm", "overloaded_error", "quota exceeded"]:
        assert E.usage_limit_hit(t), t

def test_usage_limit_negative():
    assert not E.usage_limit_hit("subtype=success turns=27 cost=$4.08")
    assert not E.usage_limit_hit("")


# ---- classify (the marker/exit contract) ----

def test_classify_ok_is_terminal():
    assert E.classify(out_exists=True, valid=True, log_text="success") == ("ok", False)

def test_classify_invalid_is_terminal_not_retried():
    # extractor produced a record but it fails the schema -> human review, do NOT auto-retry
    assert E.classify(out_exists=True, valid=False, log_text="success") == ("invalid", False)

def test_classify_usage_limit_is_retryable():
    assert E.classify(out_exists=False, valid=False, log_text="usage limit reached") == ("usage_limit", True)

def test_classify_transient_error_is_retryable():
    assert E.classify(out_exists=False, valid=False, log_text="ECONNRESET") == ("error", True)

def test_classify_quota_signal_only_matters_when_no_output():
    # a valid output wins even if the log also mentions a limit (e.g. a recovered 429 mid-run)
    assert E.classify(out_exists=True, valid=True, log_text="429 then recovered") == ("ok", False)
