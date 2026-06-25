import importlib
import pathlib
import sys

# The agent modules now live in the vaxtract package (2026-06-17 restructure).
# PKT holds that package + cancervac_packet/.
PKT = pathlib.Path(__file__).resolve().parents[1]
if str(PKT) not in sys.path:
    sys.path.insert(0, str(PKT))

# Existing tests import the modules by their old flat names (agent_core, schema, ...).
# Alias those to the package submodules so the suite keeps working without per-test
# edits. The SDK-free core modules are aliased UNCONDITIONALLY — a real import error in
# them must surface loudly, not be swallowed into a confusing downstream
# "ModuleNotFoundError: schema". Only the optional-dependency modules are guarded
# (extraction_agent needs claude_agent_sdk; figure_benchmark.* needs pymupdf/Pillow),
# and only ImportError is caught so a genuine bug there still raises.
_CORE = ["agent_core", "prompt_render", "schema_digest", "table_map", "schema", "vocab"]
_OPTIONAL = ["extraction_agent", "figure_benchmark"]
for _name in _CORE:
    sys.modules.setdefault(_name, importlib.import_module(f"vaxtract.{_name}"))
for _name in _OPTIONAL:
    if _name in sys.modules:
        continue
    try:
        sys.modules[_name] = importlib.import_module(f"vaxtract.{_name}")
    except ImportError:
        pass  # optional deps absent; only tests using these modules are affected


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "live_vision: hits the real `claude` CLI for a vision read; "
        "skipped unless RUN_LIVE_VISION=1 is set in the environment.",
    )
