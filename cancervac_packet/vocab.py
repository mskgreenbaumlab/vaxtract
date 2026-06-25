"""Back-compat shim — vocab.py moved into the vaxtract package (2026-06-17).

Re-exported here so ``cancervac_packet.vocab`` keeps resolving to the canonical module
regardless of the current working directory (see schema.py shim for rationale).

Author: Samuel Ahuno (ekwame001@gmail.com)
"""
import pathlib
import sys

_PKT_ROOT = str(pathlib.Path(__file__).resolve().parents[1])
if _PKT_ROOT not in sys.path:
    sys.path.insert(0, _PKT_ROOT)

from vaxtract import vocab as _vocab  # noqa: E402

sys.modules[__name__] = _vocab
