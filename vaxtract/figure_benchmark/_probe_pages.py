"""One-off: render full pages of main.pdf to PNGs so panel bboxes can be read.

Usage:
  python -m figure_benchmark._probe_pages --pages 1 2 5 --out /tmp/probe
Then open the PNGs, measure each panel rectangle as fractions of the page
(x0,y0,x1,y1 from top-left), and paste into panels.py bbox_frac.
"""
from __future__ import annotations

import argparse
import pathlib

from .render_panel import render_page

DEF_PDF = "agentBranch/cancerVacExtrac_claudeAi/data/raw/rojas_37165196/main.pdf"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", default=DEF_PDF)
    ap.add_argument("--pages", nargs="+", type=int, required=True)
    ap.add_argument("--out", default="/tmp/fb_probe")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args(argv)
    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    for pg in args.pages:
        img = render_page(args.pdf, pg, dpi=args.dpi)
        # 10% gridlines to read fractions off the image
        from PIL import ImageDraw
        d = ImageDraw.Draw(img)
        for f in range(1, 10):
            d.line([(img.width * f / 10, 0), (img.width * f / 10, img.height)], fill=(255, 0, 0), width=1)
            d.line([(0, img.height * f / 10), (img.width, img.height * f / 10)], fill=(255, 0, 0), width=1)
        p = out / f"page_{pg}_grid.png"
        img.save(p)
        print(f"wrote {p} ({img.width}x{img.height})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
