"""Claude Agent SDK runtime. Behaviour matches the pre-refactor extraction_agent loop."""
from __future__ import annotations

from .. import agent_core
from ..tool_registry import DENY_HOST_TOOLS, TOOLS, allowed_mcp_names
from .base import ExtractionRequest, RunResult, ToolSpec

try:
    from claude_agent_sdk import (
        ClaudeAgentOptions,
        ClaudeSDKClient,
        PermissionResultAllow,
        PermissionResultDeny,
        create_sdk_mcp_server,
        tool,
    )
except ModuleNotFoundError as exc:
    raise ImportError(
        "vaxtract's Claude runtime requires the optional 'agent' dependencies "
        "(the Claude Agent SDK + file readers). Install them with:\n"
        "    pip install 'vaxtract[agent]'\n"
        "The schema/vocab data contract (`import vaxtract.schema`) works without them."
    ) from exc

def _wrap(spec: ToolSpec):
    @tool(spec.name, spec.description, spec.input_schema)
    async def _fn(args):
        result = await spec.handler(args)
        return result.to_mcp_content()
    _fn.__name__ = spec.name
    return _fn


def build_server():
    return create_sdk_mcp_server(
        name="antvac",
        version="1.0.0",
        tools=[_wrap(t) for t in TOOLS],
    )


async def _only_antvac(tool_name, tool_input, context):
    if tool_name.startswith("mcp__antvac__"):
        return PermissionResultAllow()
    return PermissionResultDeny(
        message=(f"Tool '{tool_name}' is not available. Use ONLY the antvac tools: "
                 "read_table, read_docx, read_pdf_text, read_figure, init_record, add_entities, "
                 "add_table, clear_entities, partial_status, finalize."))


class ClaudeSdkBackend:
    name = "claude_sdk"

    async def run(self, req: ExtractionRequest, tools: list[ToolSpec]) -> RunResult:
        del tools  # catalog is the in-process MCP server built from TOOLS
        options = ClaudeAgentOptions(
            model=req.model,
            system_prompt=req.system_prompt,
            mcp_servers={"antvac": build_server()},
            tools=[],
            strict_mcp_config=True,
            allowed_tools=allowed_mcp_names(),
            disallowed_tools=list(DENY_HOST_TOOLS),
            can_use_tool=_only_antvac,
            setting_sources=[],
            permission_mode="default",
            max_turns=req.max_turns,
            max_budget_usd=req.max_budget_usd,
            max_buffer_size=32 * 1024 * 1024,
        )
        result = None
        async with ClaudeSDKClient(options=options) as client:
            await client.query(req.initial_prompt)
            async for msg in client.receive_response():
                print(msg)
                if type(msg).__name__ == "ResultMessage":
                    result = msg

        turns = int(getattr(result, "num_turns", 0) or 0) if result is not None else 0
        cost = getattr(result, "total_cost_usd", None) if result is not None else None
        is_error = bool(getattr(result, "is_error", False)) if result is not None else False
        if result is not None:
            tag = "[warn] agent loop ENDED IN ERROR" if is_error else "[info] agent loop ok"
            print(f"{tag}: subtype={getattr(result, 'subtype', '?')} "
                  f"turns={getattr(result, 'num_turns', '?')} "
                  f"cost=${(cost or 0):.2f} "
                  f"errors={getattr(result, 'errors', None)}")

        ok, gate_msg = agent_core.outer_guard(req.out_path)
        print(gate_msg)
        if not ok:
            return RunResult(ok=False, turns=turns, cost_usd=cost, error=gate_msg)
        return RunResult(ok=not is_error, turns=turns, cost_usd=cost)
