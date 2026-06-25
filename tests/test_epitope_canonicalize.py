"""Deterministic epitope canonicalization (agent_core.canonicalize_epitopes).

Root-cause fix for the run-to-run epitope-count variance: the agent merged same-identity
duplicate epitopes inconsistently (Keskin run9=100 vs gold/other runs=121, same 99 distinct
sequences). canonicalize_epitopes pins this down — one MinimalEpitope per scientific identity
with many-to-many parent_peptide_ids — WITHOUT collapsing genuinely-distinct records (e.g.
Rojas same-sequence/different-HLA). Runs deterministically at finalize.
"""
import json

import agent_core


def _epi(pid, seq, parents, hla="HLA-A*02:01", **extra):
    e = {"paper_local_id": pid, "sequence": seq, "mhc_class": "I", "hla_allele": hla,
         "is_neoantigen": True, "parent_peptide_ids": list(parents),
         "quoted_text": "q", "section_ref": "S5"}
    e.update(extra)
    return e


def test_merges_identical_except_parents_and_unions_parents():
    # Keskin case: same sequence + same HLA/gene/affinity, different parent IMP -> ONE record.
    rec = {"epitopes": [
        _epi("EPI_IMP22_X", "TGKWKLSDL", ["IMP22"], gene_symbol="DOK7"),
        _epi("EPI_IMP23_X", "TGKWKLSDL", ["IMP23"], gene_symbol="DOK7"),
    ]}
    merged = agent_core.canonicalize_epitopes(rec)
    assert merged == 1
    assert len(rec["epitopes"]) == 1
    assert rec["epitopes"][0]["paper_local_id"] == "EPI_IMP22_X"      # smallest id survives
    assert rec["epitopes"][0]["parent_peptide_ids"] == ["IMP22", "IMP23"]  # sorted union


def test_keeps_same_sequence_different_hla_separate():
    # Rojas case: same sequence, DIFFERENT HLA restriction -> genuinely distinct, NEVER merge.
    rec = {"epitopes": [
        _epi("EPI_A", "TEYKLVVVGAV", ["IMPa"], hla="HLA-B*41:01"),
        _epi("EPI_B", "TEYKLVVVGAV", ["IMPb"], hla="HLA-C*12:03"),
    ]}
    merged = agent_core.canonicalize_epitopes(rec)
    assert merged == 0
    assert {e["paper_local_id"] for e in rec["epitopes"]} == {"EPI_A", "EPI_B"}


def test_different_affinity_stays_separate():
    rec = {"epitopes": [
        _epi("EPI_A", "QTQKHLDLY", ["IMP1"], predicted_affinity={"value": 29.9, "unit": "nM", "raw": "29.9"}),
        _epi("EPI_B", "QTQKHLDLY", ["IMP2"], predicted_affinity={"value": 50.0, "unit": "nM", "raw": "50"}),
    ]}
    assert agent_core.canonicalize_epitopes(rec) == 0
    assert len(rec["epitopes"]) == 2


def test_remaps_evidence_and_curator_refs_to_survivor():
    rec = {
        "epitopes": [
            _epi("EPI_IMP22_X", "TGKWKLSDL", ["IMP22"]),
            _epi("EPI_IMP23_X", "TGKWKLSDL", ["IMP23"]),
        ],
        "evidence": [{"epitope_paper_id": "EPI_IMP23_X", "patient_paper_id": "P1"},
                     {"epitope_paper_id": "EPI_IMP22_X", "patient_paper_id": "P1"}],
        "curator_notes": [{"refs": ["EPI_IMP23_X", "P1"]}],
    }
    agent_core.canonicalize_epitopes(rec)
    assert [e["epitope_paper_id"] for e in rec["evidence"]] == ["EPI_IMP22_X", "EPI_IMP22_X"]
    assert rec["curator_notes"][0]["refs"] == ["EPI_IMP22_X", "P1"]  # non-epitope ref untouched


def test_needs_review_propagates_on_merge():
    rec = {"epitopes": [
        _epi("EPI_A", "TGKWKLSDL", ["IMP22"], needs_review=False),
        _epi("EPI_B", "TGKWKLSDL", ["IMP23"], needs_review=True),
    ]}
    agent_core.canonicalize_epitopes(rec)
    assert rec["epitopes"][0]["needs_review"] is True


def test_idempotent_and_noop_on_unique():
    rec = {"epitopes": [
        _epi("EPI_A", "TGKWKLSDL", ["IMP22"]),
        _epi("EPI_B", "TGKWKLSDL", ["IMP23"]),
        _epi("EPI_C", "VSYQGRIPY", ["IMP16"]),
    ]}
    first = agent_core.canonicalize_epitopes(rec)
    assert first == 1
    snapshot = json.dumps(rec, sort_keys=True)
    second = agent_core.canonicalize_epitopes(rec)
    assert second == 0
    assert json.dumps(rec, sort_keys=True) == snapshot   # idempotent


def test_empty_is_safe():
    rec = {"epitopes": []}
    assert agent_core.canonicalize_epitopes(rec) == 0
    assert agent_core.canonicalize_epitopes({}) == 0


def test_reference_gold_is_canonical_and_zero_loss():
    # The shipped gold references are already canonical (re-canonicalized 2026-06-05), so
    # canonicalizing them again is a NO-OP and loses nothing. Proves the invariant on real data
    # AND that the gold stays at a fixed point. (Pre-canon real-dup merging is covered by the
    # archived fixtures + the synthetic tests above.)
    import pathlib, glob
    PKT = pathlib.Path(__file__).resolve().parents[1]
    for name in ("keskin", "rojas", "li"):
        g = json.loads((PKT / "reference_records" / f"{name}_extracted.json").read_text())
        before_seqs = {e["sequence"] for e in g["epitopes"]}
        before_parents = {(e["sequence"], p) for e in g["epitopes"] for p in (e.get("parent_peptide_ids") or [])}
        rec = json.loads(json.dumps(g))
        assert agent_core.canonicalize_epitopes(rec) == 0, f"{name} gold is not canonical"
        assert {e["sequence"] for e in rec["epitopes"]} == before_seqs
        assert {(e["sequence"], p) for e in rec["epitopes"] for p in (e.get("parent_peptide_ids") or [])} == before_parents

    # And on the ARCHIVED pre-canonicalization gold (real duplicates), merging removes records
    # without losing any distinct sequence or parent link.
    arch = sorted(glob.glob(str(PKT / "reference_records" / "archived" / "*preCanon*.json")))
    for path in arch:
        g = json.loads(pathlib.Path(path).read_text())
        before_seqs = {e["sequence"] for e in g["epitopes"]}
        before_parents = {(e["sequence"], p) for e in g["epitopes"] for p in (e.get("parent_peptide_ids") or [])}
        n = agent_core.canonicalize_epitopes(g)
        assert n > 0
        assert {e["sequence"] for e in g["epitopes"]} == before_seqs
        assert {(e["sequence"], p) for e in g["epitopes"] for p in (e.get("parent_peptide_ids") or [])} == before_parents
