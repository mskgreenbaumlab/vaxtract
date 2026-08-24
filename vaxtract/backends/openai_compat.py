"""OpenAI-compatible Chat Completions + tool-call loop.

Covers xAI Grok, OpenAI GPT, Together/Inkling, NIM/vLLM Nemotron and Qwen
via base_url + api_key. The app owns confinement: only registry tools are
callable; unknown names are refused; host tools are never offered.

`read_figure` images are attached as an `image_url` content part on the
following turn — Chat Completions tool-role messages are text-only.
"""
from __future__ import annotations

import json
import os
from typing import Any

from .. import agent_core
from ..tool_registry import allowed_tools, dispatch
from .base import ExtractionRequest, RunResult, ToolResult, ToolSpec


def tools_for_openai(specs: list[ToolSpec]) -> list[dict]:
    return [{
        "type": "function",
        "function": {
            "name": spec.name,
            "description": spec.description,
            "parameters": spec.input_schema,
        },
    } for spec in specs]


def _parse_args(raw: Any) -> dict:
    if raw is None or raw == "":
        return {}
    if isinstance(raw, dict):
        return raw
    if not isinstance(raw, str):
        return {}
    try:
        val = json.loads(raw)
    except json.JSONDecodeError:
        return {"__parse_error__": raw}
    return val if isinstance(val, dict) else {}


def _text_from_result(result: ToolResult) -> str:
    if result.image_png_b64 and result.text:
        return result.text + "\n[image attached as the next user message]"
    return result.text or ""


class OpenAICompatBackend:
    name = "openai_compat"

    def __init__(self, client=None):
        self.client = client  # injected in tests; live path builds OpenAI() in run()

    def _live_client(self, req: ExtractionRequest):
        try:
            from openai import OpenAI
        except ModuleNotFoundError as exc:
            raise ImportError(
                "vaxtract's OpenAI-compatible runtime needs the 'compat' extra:\n"
                "    pip install 'vaxtract[compat]'"
            ) from exc
        env_name = req.api_key_env or "OPENAI_API_KEY"
        if env_name == "EMPTY":
            key = "EMPTY"
        else:
            key = os.environ.get(env_name)
            if not key:
                raise SystemExit(
                    f"[auth] missing {env_name}. Set that env var for this profile "
                    f"(or pass --backend claude_sdk / --profile claude)."
                )
        kwargs: dict[str, Any] = {"api_key": key}
        if req.base_url:
            kwargs["base_url"] = req.base_url
        return OpenAI(**kwargs)

    async def run(self, req: ExtractionRequest, tools: list[ToolSpec]) -> RunResult:
        specs = tools or allowed_tools()
        by_name = {s.name for s in specs if s.allowed}
        client = self.client or self._live_client(req)
        messages: list[dict] = [
            {"role": "system", "content": req.system_prompt},
            {"role": "user", "content": req.initial_prompt},
        ]
        oai_tools = tools_for_openai([s for s in specs if s.allowed])
        turns = 0
        prompt_tokens = 0
        completion_tokens = 0

        while turns < req.max_turns:
            turns += 1
            resp = client.chat.completions.create(
                model=req.model,
                messages=messages,
                tools=oai_tools,
                tool_choice="auto",
            )
            usage = getattr(resp, "usage", None)
            if usage is not None:
                prompt_tokens += int(getattr(usage, "prompt_tokens", 0) or 0)
                completion_tokens += int(getattr(usage, "completion_tokens", 0) or 0)
            choice = resp.choices[0]
            msg = choice.message
            tool_calls = list(getattr(msg, "tool_calls", None) or [])
            assistant: dict[str, Any] = {"role": "assistant", "content": getattr(msg, "content", None) or ""}
            if tool_calls:
                assistant["tool_calls"] = [
                    {
                        "id": getattr(tc, "id", f"call_{i}"),
                        "type": "function",
                        "function": {
                            "name": tc.function.name,
                            "arguments": tc.function.arguments
                            if isinstance(tc.function.arguments, str)
                            else json.dumps(tc.function.arguments or {}),
                        },
                    }
                    for i, tc in enumerate(tool_calls)
                ]
            messages.append(assistant)
            print(msg)

            if not tool_calls:
                break

            pending_images: list[str] = []
            for tc in tool_calls:
                name = tc.function.name
                call_id = getattr(tc, "id", name)
                if name not in by_name:
                    result = ToolResult(
                        text=f"ERROR: tool {name!r} is not available. Use ONLY the antvac tools.",
                        is_error=True,
                    )
                else:
                    parsed = _parse_args(tc.function.arguments)
                    if "__parse_error__" in parsed:
                        result = ToolResult(
                            text=("ERROR: tool arguments were not valid JSON. "
                                  f"Resend {name} with a JSON object. Got: {parsed['__parse_error__'][:200]!r}"),
                            is_error=True,
                        )
                    else:
                        result = await dispatch(name, parsed, vision=req.vision)
                if result.image_png_b64:
                    pending_images.append(result.image_png_b64)
                messages.append({
                    "role": "tool",
                    "tool_call_id": call_id,
                    "content": _text_from_result(result),
                })
            if pending_images:
                content: list[dict] = [{"type": "text", "text": "Figure image(s) from read_figure:"}]
                for b64 in pending_images:
                    content.append({
                        "type": "image_url",
                        "image_url": {"url": f"data:image/png;base64,{b64}"},
                    })
                messages.append({"role": "user", "content": content})

        ok, gate_msg = agent_core.outer_guard(req.out_path)
        print(gate_msg)
        if not ok:
            return RunResult(ok=False, turns=turns, error=gate_msg,
                             extra={"prompt_tokens": prompt_tokens,
                                    "completion_tokens": completion_tokens})
        if turns >= req.max_turns:
            return RunResult(ok=False, turns=turns, error=f"max_turns={req.max_turns} reached",
                             extra={"prompt_tokens": prompt_tokens,
                                    "completion_tokens": completion_tokens})
        return RunResult(ok=True, turns=turns,
                         extra={"prompt_tokens": prompt_tokens,
                                "completion_tokens": completion_tokens})
