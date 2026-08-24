"""Canned Chat Completions loop — no network, no real model."""
from __future__ import annotations

import asyncio
import json
import pathlib
from types import SimpleNamespace

from vaxtract.backends.base import ExtractionRequest
from vaxtract.backends.openai_compat import OpenAICompatBackend, tools_for_openai
from vaxtract.tool_registry import allowed_tools

_PKT = pathlib.Path(__file__).resolve().parents[1]
_REF = json.loads((_PKT / "reference_records" / "rojas_extracted.json").read_text())
_META = {k: v for k, v in _REF.items() if not isinstance(v, list)}


class _Fn(SimpleNamespace):
    pass


class _FakeCompletions:
    def __init__(self, script):
        self.script = list(script)
        self.requests = []

    def create(self, **kwargs):
        self.requests.append(kwargs)
        msg = self.script.pop(0)
        return SimpleNamespace(
            choices=[SimpleNamespace(message=msg)],
            usage=SimpleNamespace(prompt_tokens=3, completion_tokens=5),
        )


class _FakeClient:
    def __init__(self, script):
        completions = _FakeCompletions(script)
        self.chat = SimpleNamespace(completions=completions)
        self.completions = completions


def _tool_call(name, arguments, call_id="c1"):
    args = arguments if isinstance(arguments, str) else json.dumps(arguments)
    return SimpleNamespace(
        id=call_id,
        function=_Fn(name=name, arguments=args),
    )


def _req(out_path, vision=True):
    return ExtractionRequest(
        paper_dir="/tmp/paper",
        out_path=str(out_path),
        system_prompt="sys",
        initial_prompt="extract",
        model="grok-4.6",
        max_turns=8,
        vision=vision,
        api_key_env="XAI_API_KEY",
        profile_name="grok",
    )


def test_loop_init_then_stop(tmp_path):
    asyncio.run(_test_loop_init_then_stop(tmp_path))


async def _test_loop_init_then_stop(tmp_path):
    out = tmp_path / "r.json"
    script = [
        SimpleNamespace(
            content="",
            tool_calls=[_tool_call("init_record", {
                "out_path": str(out),
                "paper_meta_json": json.dumps(_META),
            }, "c-init")],
        ),
        SimpleNamespace(content="done", tool_calls=[]),
    ]
    client = _FakeClient(script)
    backend = OpenAICompatBackend(client=client)
    result = await backend.run(_req(out), allowed_tools())
    # outer_guard fails (empty record) — that's fine; the loop itself ran
    assert result.turns == 2
    assert any(t["function"]["name"] == "init_record"
               for t in client.completions.requests[0]["tools"])
    tool_names = {t["function"]["name"] for t in client.completions.requests[0]["tools"]}
    assert "Bash" not in tool_names
    assert "Read" not in tool_names
    assert (tmp_path / "r.json.partial.json").exists()


def test_unknown_tool_is_refused_not_executed(tmp_path):
    asyncio.run(_test_unknown_tool_is_refused_not_executed(tmp_path))


async def _test_unknown_tool_is_refused_not_executed(tmp_path):
    out = tmp_path / "r.json"
    script = [
        SimpleNamespace(
            content="",
            tool_calls=[_tool_call("Bash", {"command": "rm -rf /"}, "c-bash")],
        ),
        SimpleNamespace(content="ok", tool_calls=[]),
    ]
    client = _FakeClient(script)
    backend = OpenAICompatBackend(client=client)
    await backend.run(_req(out), allowed_tools())
    tool_msgs = [m for m in client.completions.requests[1]["messages"] if m.get("role") == "tool"]
    assert tool_msgs
    assert "not available" in tool_msgs[0]["content"]


def test_malformed_tool_json_returns_error_not_crash(tmp_path):
    asyncio.run(_test_malformed_tool_json_returns_error_not_crash(tmp_path))


async def _test_malformed_tool_json_returns_error_not_crash(tmp_path):
    out = tmp_path / "r.json"
    script = [
        SimpleNamespace(
            content="",
            tool_calls=[_tool_call("init_record", "{not json", "c1")],
        ),
        SimpleNamespace(content="ok", tool_calls=[]),
    ]
    client = _FakeClient(script)
    backend = OpenAICompatBackend(client=client)
    await backend.run(_req(out), allowed_tools())
    tool_msgs = [m for m in client.completions.requests[1]["messages"] if m.get("role") == "tool"]
    assert "not valid JSON" in tool_msgs[0]["content"]


def test_openai_tool_schema_is_functions():
    specs = allowed_tools()
    oai = tools_for_openai(specs)
    assert all(t["type"] == "function" for t in oai)
    names = {t["function"]["name"] for t in oai}
    assert "finalize" in names and "Bash" not in names
    fin = next(t for t in oai if t["function"]["name"] == "finalize")
    assert "allow_missing_class_ii" in fin["function"]["parameters"]["properties"]
