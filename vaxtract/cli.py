"""Console entry point for ``vaxtract``.

    vaxtract [--subscription] [--model NAME] [--profile NAME]
             [--backend claude_sdk|openai_compat] [--base-url URL]
             <paper_dir> [out.json]

BYOK Claude: set ``ANTHROPIC_API_KEY`` or pass ``--subscription`` for a logged-in
Claude plan. OpenAI-compatible profiles (grok, gpt, inkling, nemotron, qwen)
read the profile's ``api_key_env``. The agent writes a schema-validated *silver*
extraction to ``out.json`` for human sign-off. Non-Claude profiles are
experimental until the parity suite passes.

Author: Samuel Ahuno (ekwame001@gmail.com)
"""
from __future__ import annotations

import asyncio
import sys

from .backends.base import resolve_profile
from .extraction_agent import _apply_auth_mode, extract_paper, parse_argv


def main(argv: list[str] | None = None) -> None:
    argv = sys.argv[1:] if argv is None else argv
    cfg = parse_argv(argv)
    prof = resolve_profile(
        profile_name=cfg.profile, runtime=cfg.backend,
        model=cfg.model, base_url=cfg.base_url)
    print(f"[auth] {_apply_auth_mode(cfg.subscription, runtime=prof.runtime)}")
    asyncio.run(extract_paper(
        cfg.paper_dir, cfg.out_path, model=cfg.model,
        profile=cfg.profile, backend=cfg.backend, base_url=cfg.base_url))


if __name__ == "__main__":
    main()
