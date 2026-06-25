"""Pool-entity reconciliation (root-cause fix for ExtractedPeptidePool run-to-run variance).

The schema already enforces referential integrity for evidence that LINKS to a pool
(target_kind="pool" + pool_paper_id). But the agent's in-pool responses are usually tagged
target_kind="immunizing_peptide" with pool_paper_id=None, so nothing references the pool and
the pool ENTITY gets created ad hoc — present in Rojas run13, dropped in run14 (same 30
evidence rows, no pool). This backstop keys on the word "pool" in evidence quoted_text and
nudges finalize ONCE (overridable) when a patient has pooled evidence but no pool entity.
It fires only when NO pool exists for that patient, so explicit-pool papers (Keskin A-D)
never trip it.
"""
import json
import pathlib
from types import SimpleNamespace as NS

import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]
REF = json.loads((PKT / "reference_records" / "rojas_extracted.json").read_text())  # has pool_P25


def _ev(pid, q):
    return NS(patient_paper_id=pid, quoted_text=q)


# ---- unit: the detection logic (duck-typed; no schema construction needed) ----

def test_pooled_evidence_without_pool_entity_is_flagged():
    rec = NS(evidence=[_ev("P25", "DOCK11 neoantigen: De novo response in pool")], pools=[])
    assert agent_core._patients_needing_pool(rec) == ["P25"]


def test_main_text_neoantigen_pools_wording_also_flagged():  # the run14 wording
    rec = NS(evidence=[_ev("P25", "2 of 5 ELISpot responses were against 2 neoantigen pools")], pools=[])
    assert agent_core._patients_needing_pool(rec) == ["P25"]


def test_not_flagged_when_a_pool_entity_exists():
    rec = NS(evidence=[_ev("P25", "response in pool")], pools=[NS(patient_paper_id="P25")])
    assert agent_core._patients_needing_pool(rec) == []


def test_explicit_pool_paper_never_trips():  # Keskin-style: pool mentioned, pools present
    rec = NS(evidence=[_ev("Pt7", "responded primarily to pool C peptides")],
             pools=[NS(patient_paper_id="Pt7")])
    assert agent_core._patients_needing_pool(rec) == []


def test_no_flag_without_a_pool_mention():
    rec = NS(evidence=[_ev("P25", "De novo response")], pools=[])
    assert agent_core._patients_needing_pool(rec) == []


def test_intact_reference_reconciles():  # pool_P25 present + pool-target evidence
    assert agent_core._patients_needing_pool(agent_core.ExtractedPaper(**REF)) == []


# ---- integration: finalize blocks once, then the override proceeds ----

def _run14_style(rec):
    """Mimic the run14 failure: no pool entity, no pool-target evidence row, but a
    surviving per-peptide P25 row still quotes the pool. Schema-valid (no dangling ref)."""
    rec = json.loads(json.dumps(rec))
    rec["evidence"] = [e for e in rec["evidence"] if e.get("target_kind") != "pool"]
    rec["pools"] = []
    for e in rec["evidence"]:
        if e.get("patient_paper_id") == "P25":
            e["quoted_text"] = (e.get("quoted_text") or "") + " (detected in pool)"
            break
    return rec


def test_finalize_blocks_once_then_overrides(tmp_path):
    rec = _run14_style(REF)
    out = tmp_path / "r.json"
    part = tmp_path / "r.json.partial.json"

    part.write_text(json.dumps(rec))
    ok, msg = agent_core.finalize_partial(str(out), allow_missing_magnitudes=True)
    assert not ok
    assert "pool" in msg.lower() and "P25" in msg

    part.write_text(json.dumps(rec))                    # re-stage (block left it; success deletes it)
    # _run14_style strips the pool-target row, so this degraded record also drops one
    # immunogenic row (24 vs Rojas anchor 25) -> override the evidence anchor too.
    ok2, msg2 = agent_core.finalize_partial(
        str(out), allow_missing_magnitudes=True, allow_missing_pools=True,
        allow_evidence_count_mismatch=True)
    assert ok2, msg2
    assert not part.exists()                            # success cleans up the partial


def test_finalize_passes_when_pool_present(tmp_path):
    out = tmp_path / "r.json"
    (tmp_path / "r.json.partial.json").write_text(json.dumps(REF))
    ok, msg = agent_core.finalize_partial(str(out), allow_missing_magnitudes=True)
    assert ok, msg


# ---- canonical pool-EVIDENCE rule (root-cause fix for the Rojas P25 evidence-count variance) ----
# Distinct from the entity check above: here the pool ENTITY exists, but the pooled response was
# encoded as per-member rows instead of ONE pool-target row.

def _evk(pid, tk, q):  # duck evidence row with the fields the collapse-guard reads
    return NS(patient_paper_id=pid, target_kind=tk, quoted_text=q)

def test_pooled_member_rows_without_pool_row_are_flagged():
    rec = NS(pools=[NS(patient_paper_id="P25")], evidence=[
        _evk("P25", "immunizing_peptide", "De novo response in pool"),
        _evk("P25", "immunizing_peptide", "De novo response in pool"),
    ])
    assert agent_core._pool_evidence_not_collapsed(rec) == ["P25"]

def test_not_flagged_when_a_pool_target_row_consolidates():
    rec = NS(pools=[NS(patient_paper_id="P25")], evidence=[
        _evk("P25", "pool", "De novo response in pool"),
        _evk("P25", "immunizing_peptide", "De novo response in pool"),  # one deconvoluted member, OK
    ])
    assert agent_core._pool_evidence_not_collapsed(rec) == []

def test_single_pooled_member_row_is_not_flagged():  # needs >=2 to be the per-member anti-pattern
    rec = NS(pools=[NS(patient_paper_id="P25")], evidence=[
        _evk("P25", "immunizing_peptide", "De novo response in pool")])
    assert agent_core._pool_evidence_not_collapsed(rec) == []

def test_deconvoluted_rows_without_pool_wording_are_not_flagged():
    rec = NS(pools=[NS(patient_paper_id="P25")], evidence=[
        _evk("P25", "immunizing_peptide", "De novo response"),
        _evk("P25", "immunizing_peptide", "De novo response"),
    ])
    assert agent_core._pool_evidence_not_collapsed(rec) == []

def test_collapse_guard_silent_without_pool_entity():  # that case is _patients_needing_pool's job
    rec = NS(pools=[], evidence=[
        _evk("P25", "immunizing_peptide", "De novo response in pool"),
        _evk("P25", "immunizing_peptide", "De novo response in pool"),
    ])
    assert agent_core._pool_evidence_not_collapsed(rec) == []

def test_intact_reference_does_not_trip_collapse_guard():  # gold P25 = 1 pool row + deconvoluted
    assert agent_core._pool_evidence_not_collapsed(agent_core.ExtractedPaper(**REF)) == []


def _per_member_pool_style(rec):
    """Pool ENTITY kept, but its consolidating pool-target evidence row removed and >=2 of the
    patient's member rows tagged 'in pool' — the P25 anti-pattern. Schema-valid."""
    rec = json.loads(json.dumps(rec))
    rec["evidence"] = [e for e in rec["evidence"] if e.get("target_kind") != "pool"]
    n = 0
    for e in rec["evidence"]:
        if e.get("patient_paper_id") == "P25" and e.get("target_kind") == "immunizing_peptide":
            e["quoted_text"] = (e.get("quoted_text") or "") + " in pool"
            n += 1
            if n >= 2:
                break
    return rec

def test_collapse_finalize_blocks_once_then_overrides(tmp_path):
    rec = _per_member_pool_style(REF)
    out = tmp_path / "r.json"
    part = tmp_path / "r.json.partial.json"
    part.write_text(json.dumps(rec))
    ok, msg = agent_core.finalize_partial(str(out), allow_missing_magnitudes=True)
    assert not ok and "pool" in msg.lower() and "P25" in msg
    assert "allow_member_level_pool_evidence=true" in msg
    part.write_text(json.dumps(rec))
    # this collapsed variant also diverges from the Rojas immunogenic anchor (25) -> override it too.
    ok2, msg2 = agent_core.finalize_partial(
        str(out), allow_missing_magnitudes=True, allow_member_level_pool_evidence=True,
        allow_evidence_count_mismatch=True)
    assert ok2, msg2
