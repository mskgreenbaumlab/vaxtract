import base64
import io
import pathlib

import pytest
from PIL import Image

import agent_core

PDF = pathlib.Path(
    "agentBranch/cancerVacExtrac_claudeAi/data/raw/rojas_37165196/main.pdf"
)
needs_pdf = pytest.mark.skipif(not PDF.exists(), reason="Rojas main.pdf not present")


def _decode(b64):
    return Image.open(io.BytesIO(base64.b64decode(b64)))


@needs_pdf
def test_render_full_page_caps_long_side():
    b64, w, h = agent_core.render_figure_image(str(PDF), page=2, max_side=1200)
    img = _decode(b64)
    assert img.format == "PNG"
    assert max(img.width, img.height) <= 1200
    assert (w, h) == (img.width, img.height)


@needs_pdf
def test_render_region_is_a_subcrop():
    full, fw, fh = agent_core.render_figure_image(str(PDF), page=2, max_side=4000)
    crop, cw, ch = agent_core.render_figure_image(
        str(PDF), page=2, region=(0.8, 0.36, 0.98, 0.49), max_side=4000
    )
    # the region crop covers a fraction of the page → fewer pixels than the full page
    assert cw * ch < fw * fh


def test_bad_region_raises_valueerror():
    with pytest.raises(ValueError):
        agent_core.render_figure_image(str(PDF), page=2, region=(0.5, 0.5, 0.4, 0.6))  # x1<x0
