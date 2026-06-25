import json
import pathlib

import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]
REF = json.loads((PKT / "reference_records" / "rojas_extracted.json").read_text())
META = {k: v for k, v in REF.items() if not isinstance(v, list)}
SECTIONS = [s for s in agent_core.SECTION_MODEL if REF.get(s)]


def _seed(out, mutate=None):
    agent_core.init_partial(str(out), json.dumps(META))
    for s in SECTIONS:
        items = REF[s]
        if mutate:
            items = mutate(s, json.loads(json.dumps(items)))
        ok, msg = agent_core.append_section(str(out), s, json.dumps(items))
        assert ok, (s, msg)


def test_reference_with_raw_magnitudes_finalizes(tmp_path):
    # run-5 reference: immunogenic rows carry a raw → magnitude present → passes the nudge
    out = tmp_path / "r.json"
    _seed(out)
    ok, msg = agent_core.finalize_partial(str(out))
    assert ok, msg


def _null_all_magnitudes(section, items):
    if section == "evidence":
        for e in items:
            e["magnitude"] = None
    return items


def test_null_magnitudes_block_finalize_once(tmp_path):
    out = tmp_path / "r.json"
    _seed(out, mutate=_null_all_magnitudes)
    ok, msg = agent_core.finalize_partial(str(out))
    assert not ok
    assert "magnitude" in msg.lower()
    assert (tmp_path / "r.json.partial.json").exists()   # partial kept so the agent can fix


def test_override_lets_null_magnitudes_through(tmp_path):
    out = tmp_path / "r.json"
    _seed(out, mutate=_null_all_magnitudes)
    ok, msg = agent_core.finalize_partial(str(out), allow_missing_magnitudes=True)
    assert ok, msg
    assert out.exists()
    assert not (tmp_path / "r.json.partial.json").exists()  # cleaned up on success


def test_a_single_raw_magnitude_satisfies_a_row(tmp_path):
    # a deliberate "not reported" raw counts as present (only fully-null trips the nudge)
    def one_raw(section, items):
        if section == "evidence":
            for e in items:
                e["magnitude"] = {"unit": "unknown", "raw": "not reported", "tier": "reported"}
        return items
    out = tmp_path / "r.json"
    _seed(out, mutate=one_raw)
    ok, msg = agent_core.finalize_partial(str(out))
    assert ok, msg
