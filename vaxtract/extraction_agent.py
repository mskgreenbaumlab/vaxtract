#!/usr/bin/env python3
"""extraction_agent.py — paper dir -> validated JSON via a tool-using agent loop.

WHAT MAKES THIS AN AGENT: the model is given TOOLS and a LOOP. It reads the
tables/text itself, drafts candidate JSON, calls validate(), reads the schema's
errors, and FIXES its own output — repeating until the record passes the packet
schema or it gives up (-> outer-guard quarantine).

Pure logic lives in agent_core.py / prompt_render.py (SDK-free, unit-tested).
Tools are defined once in tool_registry.py. This file is the CLI/runtime shell:
it picks a backend (Claude Agent SDK, or OpenAI-compatible Chat Completions)
and runs the loop.

Setup:
    pip install 'vaxtract[agent]'            # Claude runtime
    pip install 'vaxtract[compat]'           # Grok / GPT / Inkling / vLLM
    export ANTHROPIC_API_KEY=...             # Claude BYOK (or --subscription)
    export XAI_API_KEY=...                   # --profile grok
Usage:
    vaxtract <paper_dir> [out.json]
    vaxtract --profile grok <paper_dir> [out.json]

Author: Samuel Ahuno (ekwame001@gmail.com)
"""
from __future__ import annotations

import os
import pathlib
import sys
from dataclasses import dataclass

from . import agent_core
from .backends.base import ExtractionRequest, resolve_profile
from .prompt_render import build_system_prompt, field_guidance_only
from .schema_digest import build_schema_digest
from .tool_registry import (
    DENY_HOST_TOOLS,
    allowed_mcp_names,
    handle_read_figure,
)

# 1M-context variant (B, 2026-06-04): Keskin overflowed the 200K window -> 2 compactions
# -> re-orientation churn (extra turns/thinking). The [1m] model removes compaction entirely.
# Trade-off: 1M-tier pricing applies only ABOVE 200K, so small papers (Rojas, ~64K/turn) are
# unaffected; only a paper big enough to have compacted pays the premium -- exactly when it
# helps. (The SDK `betas=['context-1m-2025-08-07']` flag is Sonnet-only and does NOT enable
# 1M on Opus; the context window is selected via this model id.)
MODEL = "claude-opus-4-8[1m]"
MAX_TURNS = 120
MAX_BUDGET_USD = 18.0

_DENY_HOST_TOOLS = DENY_HOST_TOOLS
_ANTVAC_TOOLS = allowed_mcp_names()

SYSTEM_PROMPT = build_system_prompt(
    field_guidance_only(agent_core.DELTAS_TEXT),
    build_schema_digest(agent_core.schema, agent_core.vocab),
    agent_core.vocab,
)

USAGE = (
    "usage: vaxtract [--subscription] [--model NAME] [--profile NAME] "
    "[--backend claude_sdk|openai_compat] [--base-url URL] "
    "<paper_dir> [out.json]"
)


async def read_figure(args):
    """MCP-shaped wrapper kept for tests that call the figure tool directly."""
    result = await handle_read_figure(args)
    return result.to_mcp_content()


@dataclass
class CliArgs:
    paper_dir: str
    out_path: str
    subscription: bool = False
    model: str | None = None
    profile: str | None = None
    backend: str | None = None
    base_url: str | None = None


def parse_argv(argv: list[str]) -> CliArgs:
    """Parse argv (excluding the program name)."""
    if "--help" in argv or "-h" in argv:
        profiles = "claude, grok, gpt, inkling, nemotron, qwen"
        raise SystemExit(
            USAGE + "\n"
            "  --profile NAME     named endpoint (" + profiles + "); default claude\n"
            "  --backend NAME     claude_sdk | openai_compat (overrides the profile runtime)\n"
            "  --model NAME       model id (overrides the profile)\n"
            "  --base-url URL     OpenAI-compatible endpoint (Grok, Together, vLLM, ...)\n"
            "  --subscription     Claude plan quota (Claude runtime only)\n"
            "Non-Claude profiles are experimental until the parity suite passes.\n"
            "Keep-files under outputs/classii_* are frozen; do not overwrite them."
        )
    subscription = "--subscription" in argv
    args = [a for a in argv if a != "--subscription"]
    model: str | None = None
    profile: str | None = None
    backend: str | None = None
    base_url: str | None = None
    positional: list[str] = []
    i = 0
    while i < len(args):
        if args[i] == "--model":
            if i + 1 >= len(args):
                raise SystemExit(USAGE)
            model = args[i + 1]
            i += 2
            continue
        if args[i] == "--profile":
            if i + 1 >= len(args):
                raise SystemExit(USAGE)
            profile = args[i + 1]
            i += 2
            continue
        if args[i] == "--backend":
            if i + 1 >= len(args):
                raise SystemExit(USAGE)
            backend = args[i + 1]
            i += 2
            continue
        if args[i] == "--base-url":
            if i + 1 >= len(args):
                raise SystemExit(USAGE)
            base_url = args[i + 1]
            i += 2
            continue
        positional.append(args[i])
        i += 1
    if not positional:
        raise SystemExit(USAGE)
    paper_dir = positional[0]
    out = positional[1] if len(positional) > 1 else "newpaper_extracted.json"
    return CliArgs(
        paper_dir=paper_dir, out_path=out, subscription=subscription,
        model=model, profile=profile, backend=backend, base_url=base_url,
    )


def _parse_cli(argv):
    """Parse argv → (paper_dir, out_path, subscription, model).

    Kept as a 4-tuple for existing tests. New flags live on parse_argv().
    """
    cfg = parse_argv(argv)
    return cfg.paper_dir, cfg.out_path, cfg.subscription, cfg.model or MODEL


def _apply_auth_mode(subscription, runtime: str = "claude_sdk"):
    """Decide how authentication works and return a label.

    `--subscription` is Claude-only: drop ANTHROPIC_API_KEY so the spawned
    `claude` CLI uses plan quota. Combining it with openai_compat is a hard error.
    """
    if subscription and runtime != "claude_sdk":
        raise SystemExit(
            "[auth] --subscription is Claude-only; do not combine it with "
            f"--backend {runtime} / a non-Claude --profile"
        )
    if runtime != "claude_sdk":
        return f"api-key ({runtime}; see the profile's api_key_env)"
    if subscription:
        os.environ.pop("ANTHROPIC_API_KEY", None)
        return "subscription (ANTHROPIC_API_KEY unset; claude CLI uses your plan quota)"
    return ("api-key (ANTHROPIC_API_KEY present)"
            if os.environ.get("ANTHROPIC_API_KEY")
            else "subscription (no ANTHROPIC_API_KEY in env)")


def _list_paper_files(paper_dir: str) -> list[str]:
    files = [str(p) for p in pathlib.Path(paper_dir).rglob("*")
             if p.suffix.lower() in (".pdf", ".xlsx", ".docx")]
    if not files:
        print(f"[abort] no .pdf/.xlsx/.docx files found in {paper_dir}")
        sys.exit(1)
    return files


async def extract_paper(
    paper_dir: str,
    out_path: str,
    model: str | None = None,
    *,
    profile: str | None = None,
    backend: str | None = None,
    base_url: str | None = None,
) -> None:
    files = _list_paper_files(paper_dir)
    prof = resolve_profile(
        profile_name=profile, runtime=backend, model=model, base_url=base_url)
    if model:
        prof.model = model
    print(f"[model] {prof.model}")
    print(f"[profile] {prof.name}  [runtime] {prof.runtime}  "
          f"[vision] {str(prof.vision).lower()}"
          + ("  [experimental]" if prof.experimental else ""))
    prompt = (f"Extract this paper into {out_path}. Files available:\n"
              + "\n".join(files) +
              "\nStart by reading the neoantigen/ELISpot supplementary tables.")
    req = ExtractionRequest(
        paper_dir=paper_dir,
        out_path=out_path,
        system_prompt=SYSTEM_PROMPT,
        initial_prompt=prompt,
        model=prof.model,
        max_turns=MAX_TURNS,
        max_budget_usd=MAX_BUDGET_USD,
        vision=prof.vision,
        max_context=prof.max_context,
        base_url=prof.base_url,
        api_key_env=prof.api_key_env,
        profile_name=prof.name,
    )
    from .backends import select_backend
    from .tool_registry import allowed_tools
    runner = select_backend(prof.runtime)
    result = await runner.run(req, allowed_tools())
    if not result.ok:
        sys.exit(1)
