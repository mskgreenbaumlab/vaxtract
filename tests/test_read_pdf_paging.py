"""read_pdf_text_from paging window.

The default window keeps each result under the CLI's per-tool-result cap so the text is
never spilled to a file -- the agent has no host Read to re-open a spill (host tools are
denied). A truncated window reports the next offset so the agent can page through.
"""
import pathlib

import agent_core


def test_default_window_caps_and_reports_next_offset(monkeypatch):
    big = "X" * 100_000
    monkeypatch.setattr(agent_core, "read_pdf_text_from",
                        agent_core.read_pdf_text_from)  # ensure real fn

    # stub PdfReader so no real PDF is needed
    import types
    fake = types.SimpleNamespace(pages=[types.SimpleNamespace(extract_text=lambda: big)])
    import sys
    pypdf = types.ModuleType("pypdf")
    pypdf.PdfReader = lambda _p: fake
    monkeypatch.setitem(sys.modules, "pypdf", pypdf)

    out = agent_core.read_pdf_text_from("whatever.pdf")
    assert len(out) < 41_000  # 40k window + short suffix, well under the spill cap
    assert "offset=40000" in out
    assert "of 100000" in out


def test_offset_reads_next_chunk_and_marks_end(monkeypatch):
    body = "".join(chr(65 + (i % 26)) for i in range(50_000))
    import types, sys
    fake = types.SimpleNamespace(pages=[types.SimpleNamespace(extract_text=lambda: body)])
    pypdf = types.ModuleType("pypdf")
    pypdf.PdfReader = lambda _p: fake
    monkeypatch.setitem(sys.modules, "pypdf", pypdf)

    out = agent_core.read_pdf_text_from("x.pdf", max_chars=40_000, offset=40_000)
    assert out.startswith(body[40_000:40_050])
    assert "end of text" in out and "of 50000" in out
