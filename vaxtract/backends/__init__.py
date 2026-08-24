"""Extraction runtimes. Claude SDK and OpenAI-compatible HTTP are the two.

Import backends lazily so `import vaxtract.schema` and the CLI parser do not
require `claude_agent_sdk` or `openai`.
"""
from __future__ import annotations

from .base import (
    ExtractionRequest,
    Profile,
    RunResult,
    ToolResult,
    ToolSpec,
    load_profiles,
    resolve_profile,
)


def select_backend(runtime: str):
    """Return a backend instance for `claude_sdk` or `openai_compat`."""
    if runtime == "claude_sdk":
        from .claude_sdk import ClaudeSdkBackend
        return ClaudeSdkBackend()
    if runtime == "openai_compat":
        from .openai_compat import OpenAICompatBackend
        return OpenAICompatBackend()
    raise ValueError(
        f"unknown runtime {runtime!r}; expected 'claude_sdk' or 'openai_compat'"
    )


__all__ = [
    "ExtractionRequest",
    "Profile",
    "RunResult",
    "ToolResult",
    "ToolSpec",
    "load_profiles",
    "resolve_profile",
    "select_backend",
]
