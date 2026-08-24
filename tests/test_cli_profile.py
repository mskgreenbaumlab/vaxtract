import os

import pytest

import extraction_agent as ea
from vaxtract.backends.base import load_profiles, resolve_profile


def test_profiles_yaml_has_expected_runtimes():
    profiles = load_profiles()
    assert profiles["claude"].runtime == "claude_sdk"
    assert profiles["claude"].experimental is False
    for name in ("grok", "gpt", "inkling", "nemotron", "qwen"):
        assert profiles[name].runtime == "openai_compat"
        assert profiles[name].experimental is True
    assert profiles["grok"].base_url == "https://api.x.ai/v1"
    assert profiles["grok"].api_key_env == "XAI_API_KEY"
    assert profiles["nemotron"].vision is False
    assert profiles["qwen"].vision is False


def test_parse_profile_and_backend_flags():
    cfg = ea.parse_argv(["--profile", "grok", "--model", "grok-4.6", "data/raw/foo", "out.json"])
    assert cfg.profile == "grok"
    assert cfg.model == "grok-4.6"
    assert cfg.paper_dir == "data/raw/foo"
    assert cfg.subscription is False
    cfg = ea.parse_argv(["--backend", "openai_compat", "--base-url", "http://127.0.0.1:8000/v1",
                         "data/raw/foo"])
    assert cfg.backend == "openai_compat"
    assert cfg.base_url == "http://127.0.0.1:8000/v1"


def test_parse_cli_four_tuple_still_works():
    paper_dir, out, sub, model = ea._parse_cli(["data/raw/foo"])
    assert (paper_dir, out, sub, model) == ("data/raw/foo", "newpaper_extracted.json", False, ea.MODEL)


def test_help_lists_profiles():
    with pytest.raises(SystemExit) as ei:
        ea.parse_argv(["--help"])
    msg = str(ei.value)
    assert "--profile" in msg
    assert "grok" in msg
    assert "experimental" in msg


def test_subscription_rejected_on_openai_compat():
    with pytest.raises(SystemExit) as ei:
        ea._apply_auth_mode(True, runtime="openai_compat")
    assert "Claude-only" in str(ei.value)


def test_subscription_still_drops_anthropic_key(monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-ant-test")
    mode = ea._apply_auth_mode(True, runtime="claude_sdk")
    assert "ANTHROPIC_API_KEY" not in os.environ
    assert "subscription" in mode.lower()


def test_resolve_profile_unknown():
    with pytest.raises(ValueError, match="unknown profile"):
        resolve_profile(profile_name="claude-but-typo")


def test_resolve_profile_overlay_model():
    p = resolve_profile(profile_name="grok", model="grok-4")
    assert p.model == "grok-4"
    assert p.runtime == "openai_compat"
