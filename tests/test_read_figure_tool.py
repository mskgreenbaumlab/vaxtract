import asyncio
import pathlib

import pytest

import extraction_agent as ea

PDF = pathlib.Path(
    "agentBranch/cancerVacExtrac_claudeAi/data/raw/rojas_37165196/main.pdf"
)
needs_pdf = pytest.mark.skipif(not PDF.exists(), reason="Rojas main.pdf not present")


def _call(**args):
    fn = ea.read_figure.handler if hasattr(ea.read_figure, "handler") else ea.read_figure
    return asyncio.get_event_loop().run_until_complete(fn(args))


@needs_pdf
def test_read_figure_returns_image_block_for_a_page():
    out = _call(path=str(PDF), page=2, what="Figure 1G ELISpot")
    blocks = out["content"]
    kinds = [b.get("type") for b in blocks]
    assert "image" in kinds          # an image block is present
    assert "text" in kinds           # plus the recording reminder
    img = next(b for b in blocks if b.get("type") == "image")
    assert img.get("data")           # base64 PNG payload present
    assert img.get("mimeType") == "image/png"


@needs_pdf
def test_read_figure_region_path_also_returns_image():
    out = _call(path=str(PDF), page=2, what="zoom", region=[0.8, 0.36, 0.98, 0.49])
    assert any(b.get("type") == "image" for b in out["content"])


def test_read_figure_bad_path_returns_text_error_not_crash():
    out = _call(path="/no/such.pdf", page=0, what="x")
    assert out["content"][0]["type"] == "text"
    assert "ERROR" in out["content"][0]["text"]
