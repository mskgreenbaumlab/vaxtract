"""vaxtract — schema-validated extraction of neoantigen cancer-vaccine
immunogenicity data from primary papers.

The agent loop is model-agnostic: Claude Agent SDK (default) or any
OpenAI-compatible Chat Completions endpoint (Grok, GPT, Inkling, vLLM).

Public API:
    from vaxtract import extract_paper      # async (paper_dir, out_path) -> None

``extract_paper`` is imported lazily so that ``import vaxtract.schema`` (the
data contract) and the other SDK-free modules work without a model SDK
installed.

Author: Samuel Ahuno (ekwame001@gmail.com)
"""
from __future__ import annotations

__version__ = "0.3.0"
__all__ = ["extract_paper", "__version__"]


def __getattr__(name: str):  # PEP 562 lazy attribute
    if name == "extract_paper":
        from .extraction_agent import extract_paper
        return extract_paper
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
