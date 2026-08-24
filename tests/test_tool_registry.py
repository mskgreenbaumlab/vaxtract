import asyncio
import inspect

import agent_core
from vaxtract.tool_registry import (
    DENY_HOST_TOOLS,
    allowed_mcp_names,
    allowed_tools,
    dispatch,
    finalize_override_names,
    get_tool,
)


def test_catalog_has_expected_tools():
    names = {t.name for t in allowed_tools()}
    for required in ("read_table", "read_pdf_text", "survey_sources", "read_figure",
                     "init_record", "add_entities", "add_table", "finalize",
                     "build_pools", "build_pool_evidence", "build_crossreactivity_evidence"):
        assert required in names
    assert "read_docx" not in names
    assert get_tool("read_docx").allowed is False


def test_mcp_allowlist_prefix():
    names = allowed_mcp_names()
    assert all(n.startswith("mcp__antvac__") for n in names)
    assert "mcp__antvac__finalize" in names
    assert "mcp__antvac__read_docx" not in names


def test_finalize_schema_tracks_finalize_partial_signature():
    params = [p for p in inspect.signature(agent_core.finalize_partial).parameters
              if p.startswith("allow_")]
    assert set(finalize_override_names()) == set(params)
    props = get_tool("finalize").input_schema["properties"]
    for p in params:
        assert p in props


def test_host_tools_denied():
    for t in ("Grep", "Read", "Glob", "Bash", "Write"):
        assert t in DENY_HOST_TOOLS
    assert not any(t.startswith("mcp__antvac__") for t in DENY_HOST_TOOLS)


def test_dispatch_refuses_unknown_tool():
    result = asyncio.run(dispatch("Bash", {}))
    assert result.is_error
    assert "not available" in result.text


def test_dispatch_blocks_read_docx_allowlist():
    result = asyncio.run(dispatch("read_docx", {"path": "/x"}))
    assert result.is_error
    assert "allowlist" in result.text


def test_vision_false_does_not_render_a_figure():
    result = asyncio.run(
        dispatch("read_figure", {"path": "/no.pdf", "page": 0, "what": "x"}, vision=False))
    assert result.is_error
    assert "vision=false" in result.text
    assert result.image_png_b64 is None
