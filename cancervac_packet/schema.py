"""Back-compat shim — schema.py moved into the vaxtract package (2026-06-17).

Re-exported here so ``cancervac_packet.schema`` (and dev tools like predeploy_gate.py /
make_report.py / schema_overview.py that put cancervac_packet on sys.path) keep
resolving to the single canonical module — regardless of the current working directory.

Author: Samuel Ahuno (ekwame001@gmail.com)
"""
import pathlib
import sys

# Ensure the package root (parent of cancervac_packet/, which holds vaxtract/)
# is importable even when only cancervac_packet was added to sys.path.
_PKT_ROOT = str(pathlib.Path(__file__).resolve().parents[1])
if _PKT_ROOT not in sys.path:
    sys.path.insert(0, _PKT_ROOT)

from vaxtract import schema as _schema  # noqa: E402

sys.modules[__name__] = _schema
