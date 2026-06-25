"""Console entry point for ``vaxtract``.

    vaxtract [--subscription] <paper_dir> [out.json]

BYOK: set ``ANTHROPIC_API_KEY`` (pay-per-token) or pass ``--subscription`` to use a
logged-in Claude plan via the ``claude`` CLI. The agent reads the PDF/XLSX/DOCX files
in ``<paper_dir>`` and writes a schema-validated *silver* extraction to ``out.json``
for human sign-off.

Author: Samuel Ahuno (ekwame001@gmail.com)
"""
from __future__ import annotations

import asyncio
import sys

from .extraction_agent import _apply_auth_mode, _parse_cli, extract_paper


def main(argv: list[str] | None = None) -> None:
    argv = sys.argv[1:] if argv is None else argv
    paper_dir, out, subscription = _parse_cli(argv)
    print(f"[auth] {_apply_auth_mode(subscription)}")
    asyncio.run(extract_paper(paper_dir, out))


if __name__ == "__main__":
    main()
