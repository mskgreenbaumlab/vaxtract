"""Provider-neutral types for the extraction loop.

`agent_core` is the logic boundary; this module is the *provider* boundary.
Tools are defined once as ToolSpec; each runtime adapts them to its API.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Awaitable, Callable, Protocol

_PKG = Path(__file__).resolve().parent.parent
PROFILES_PATH = _PKG / "profiles.yaml"


def _coerce(raw: str) -> Any:
    s = raw.strip().strip('"').strip("'")
    if s.lower() in ("true", "yes"):
        return True
    if s.lower() in ("false", "no"):
        return False
    if s.lower() in ("null", "none", "~", ""):
        return None
    try:
        if s.isdigit() or (s.startswith("-") and s[1:].isdigit()):
            return int(s)
        return float(s)
    except ValueError:
        return s


def load_yaml_map(path: Path) -> dict[str, dict[str, Any]]:
    """Load a restricted two-level YAML map (profile/paper name -> key: value).

    No extra dependency. Nested lists and multiline values are not supported.
    """
    text = path.read_text(encoding="utf-8")
    out: dict[str, dict[str, Any]] = {}
    current: str | None = None
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if line.startswith(" ") or line.startswith("\t"):
            if current is None:
                raise ValueError(f"{path}: indented line with no parent key: {line!r}")
            if ":" not in line:
                raise ValueError(f"{path}: expected key: value, got {line!r}")
            key, _, val = line.strip().partition(":")
            out[current][key.strip()] = _coerce(val)
            continue
        if line.endswith(":"):
            current = line[:-1].strip()
            out[current] = {}
            continue
        raise ValueError(f"{path}: expected a section header or indented key, got {line!r}")
    return out


@dataclass
class ToolResult:
    """Provider-neutral tool return. Backends render this to MCP / Chat Completions."""
    text: str | None = None
    image_png_b64: str | None = None
    is_error: bool = False

    def to_mcp_content(self) -> dict:
        content: list[dict] = []
        if self.image_png_b64:
            content.append({"type": "image", "data": self.image_png_b64, "mimeType": "image/png"})
        if self.text is not None:
            content.append({"type": "text", "text": self.text})
        if not content:
            content.append({"type": "text", "text": ""})
        return {"content": content}


@dataclass
class ToolSpec:
    name: str
    description: str
    input_schema: dict
    handler: Callable[[dict], Awaitable[ToolResult]]
    allowed: bool = True  # False = registered but not on the model allowlist (read_docx)


@dataclass
class Profile:
    name: str
    runtime: str  # claude_sdk | openai_compat
    model: str
    vision: bool = True
    max_context: int | None = None
    base_url: str | None = None
    api_key_env: str | None = None
    experimental: bool = True

    @classmethod
    def from_mapping(cls, name: str, raw: dict[str, Any]) -> "Profile":
        runtime = raw.get("runtime")
        model = raw.get("model")
        if not runtime or not model:
            raise ValueError(f"profile {name!r} needs runtime and model")
        return cls(
            name=name,
            runtime=str(runtime),
            model=str(model),
            vision=bool(raw.get("vision", True)),
            max_context=raw.get("max_context"),
            base_url=raw.get("base_url"),
            api_key_env=raw.get("api_key_env"),
            experimental=bool(raw.get("experimental", True)),
        )


def load_profiles(path: Path | None = None) -> dict[str, Profile]:
    raw = load_yaml_map(path or PROFILES_PATH)
    return {name: Profile.from_mapping(name, body) for name, body in raw.items()}


def resolve_profile(
    *,
    profile_name: str | None = None,
    runtime: str | None = None,
    model: str | None = None,
    base_url: str | None = None,
) -> Profile:
    """CLI resolution: --profile NAME, with optional --backend/--model/--base-url overlays."""
    profiles = load_profiles()
    name = profile_name or "claude"
    if name not in profiles:
        known = ", ".join(sorted(profiles))
        raise ValueError(f"unknown profile {name!r}; known: {known}")
    p = profiles[name]
    if runtime:
        p.runtime = runtime
    if model:
        p.model = model
    if base_url:
        p.base_url = base_url
    return p


@dataclass
class ExtractionRequest:
    paper_dir: str
    out_path: str
    system_prompt: str
    initial_prompt: str
    model: str
    max_turns: int
    max_budget_usd: float | None = None
    vision: bool = True
    max_context: int | None = None
    base_url: str | None = None
    api_key_env: str | None = None
    profile_name: str = "claude"


@dataclass
class RunResult:
    ok: bool
    turns: int = 0
    cost_usd: float | None = None
    error: str | None = None
    extra: dict[str, Any] = field(default_factory=dict)


class ExtractionBackend(Protocol):
    name: str

    async def run(self, req: ExtractionRequest, tools: list[ToolSpec]) -> RunResult: ...
