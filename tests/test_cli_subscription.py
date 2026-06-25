import os

import pytest

import extraction_agent as ea


def test_parse_cli_plain():
    paper_dir, out, sub = ea._parse_cli(["data/raw/foo"])
    assert paper_dir == "data/raw/foo"
    assert out == "newpaper_extracted.json"  # default
    assert sub is False


def test_parse_cli_with_out_and_subscription_anywhere():
    # --subscription is position-independent and stripped from positionals
    paper_dir, out, sub = ea._parse_cli(["--subscription", "data/raw/foo", "out.json"])
    assert (paper_dir, out, sub) == ("data/raw/foo", "out.json", True)
    paper_dir, out, sub = ea._parse_cli(["data/raw/foo", "out.json", "--subscription"])
    assert (paper_dir, out, sub) == ("data/raw/foo", "out.json", True)


def test_parse_cli_requires_paper_dir():
    with pytest.raises(SystemExit):
        ea._parse_cli(["--subscription"])


def test_apply_auth_mode_subscription_drops_api_key(monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-ant-test")
    mode = ea._apply_auth_mode(True)
    assert "ANTHROPIC_API_KEY" not in os.environ   # popped → claude CLI uses the subscription
    assert "subscription" in mode.lower()


def test_apply_auth_mode_default_keeps_api_key(monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-ant-test")
    mode = ea._apply_auth_mode(False)
    assert os.environ.get("ANTHROPIC_API_KEY") == "sk-ant-test"  # untouched
    assert "api" in mode.lower()


def test_host_file_tools_are_hard_denied():
    # The Grep/Read leak fix: read-only host tools bypass can_use_tool in default mode,
    # so they must be on the disallowed_tools denylist to actually be blocked.
    for t in ("Grep", "Read", "Glob", "Bash", "Write"):
        assert t in ea._DENY_HOST_TOOLS


def test_denylist_does_not_block_antvac_tools():
    assert not any(t.startswith("mcp__antvac__") for t in ea._DENY_HOST_TOOLS)
