import json
import pathlib

import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]
REF = json.loads((PKT / "reference_records" / "rojas_extracted.json").read_text())


def test_matching_counts_pass(tmp_path):
    # reference: sum(n_peptides_synthesized) == len(immunizing_peptides)
    good = tmp_path / "good.json"
    good.write_text(json.dumps(REF))
    ok, msg = agent_core.outer_guard(str(good))
    assert ok, msg


def test_declared_more_than_present_fails(tmp_path):
    # a silent drop looks like: patients declare more synthesized than peptides exist
    rec = json.loads(json.dumps(REF))
    rec["patients"][0]["n_peptides_synthesized"] += 5   # declared 5 the record doesn't hold
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(rec))
    ok, msg = agent_core.outer_guard(str(bad))
    assert not ok
    assert "peptide count" in msg.lower() or "synthesized" in msg.lower()
    assert bad.exists()                                 # schema-valid → not quarantined
