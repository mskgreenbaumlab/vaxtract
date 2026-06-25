import json
import pathlib

import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]
REF = json.loads((PKT / "reference_records" / "rojas_extracted.json").read_text())


def test_linked_record_passes_guard(tmp_path):
    good = tmp_path / "good.json"
    good.write_text(json.dumps(REF))
    ok, msg = agent_core.outer_guard(str(good))
    assert ok, msg


def test_orphaned_epitopes_fail_guard_without_quarantine(tmp_path):
    rec = json.loads(json.dumps(REF))
    for e in rec["epitopes"]:
        e["parent_peptide_ids"] = []           # break the link, as run 6 did
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(rec))
    ok, msg = agent_core.outer_guard(str(bad))
    assert not ok
    assert "parent_peptide_ids" in msg
    assert bad.exists()                          # schema-valid → not quarantined, agent can fix


def test_dangling_epitope_reference_fails_guard(tmp_path):
    rec = json.loads(json.dumps(REF))
    rec["epitopes"][0]["parent_peptide_ids"] = ["NO_SUCH_PEPTIDE"]
    bad = tmp_path / "bad2.json"
    bad.write_text(json.dumps(rec))
    ok, msg = agent_core.outer_guard(str(bad))
    assert not ok
    assert "unknown" in msg.lower() or "parent_peptide_ids" in msg
