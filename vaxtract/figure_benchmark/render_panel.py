"""Render a main.pdf page to a raster and crop a panel by fractional bbox.

bbox_frac is (x0,y0,x1,y1) as fractions of page width/height so it is
resolution-independent. Output PNG is content-addressed: filename embeds a hash
of (page, dpi, bbox) so re-runs are cache hits. The benchmark uses exactly one
PDF, so pdf identity is intentionally NOT in the key; add it before reusing this
across multiple PDFs writing to the same out_dir.
"""
from __future__ import annotations

import hashlib
import pathlib

import fitz  # PyMuPDF
from PIL import Image

from .model import Panel


def render_page(pdf_path, page: int, dpi: int = 300) -> Image.Image:
    doc = fitz.open(pdf_path)
    try:
        pg = doc[page]
        zoom = dpi / 72.0
        pix = pg.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)
        return Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    finally:
        doc.close()


def crop_frac(img: Image.Image, bbox_frac: tuple) -> Image.Image:
    x0, y0, x1, y1 = bbox_frac
    W, H = img.width, img.height
    box = (int(x0 * W), int(y0 * H), int(x1 * W), int(y1 * H))
    return img.crop(box)


def _panel_hash(panel: Panel, dpi: int) -> str:
    h = hashlib.sha256()
    h.update(f"{panel.page}|{dpi}|{panel.bbox_frac}".encode())
    return h.hexdigest()[:16]


def render_panel(panel: Panel, pdf_path, out_dir, dpi: int | None = None) -> pathlib.Path:
    if panel.bbox_frac is None:
        raise ValueError(f"panel {panel.panel_id!r} has no bbox_frac set")
    dpi = dpi if dpi is not None else panel.render_dpi
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{panel.panel_id}_{_panel_hash(panel, dpi)}.png"
    if out.exists():
        return out
    img = render_page(pdf_path, panel.page, dpi=dpi)
    crop_frac(img, panel.bbox_frac).save(out)
    return out
