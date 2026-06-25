"""v2.8 P19 — NeoantigenCandidate funnel (prediction -> selection -> outcome).

Covers the ported candidate wiring and the THREE design-review corrections:
  (1) PrioritizationScore is a SIBLING of Measurement, NOT `ranking_scores:
      list[Measurement]` — a non-affinity funnel score keeps `value` as a
      queryable FLOAT instead of degrading to value=None+raw (a lossless string).
  (2) SOFT bridge sequence-consistency: a selected candidate whose
      selected_peptide_id resolves to an IMP of a DIFFERENT sequence is flagged
      (block-once / overridable), never hard-rejected.
  (3) SOFT funnel-completeness signal: candidates present with
      n_predicted_reported unset is nudged (overridable).

HARD structural guards (candidate id disjointness, bridge resolves to an IMP,
candidate evidence target, unique-id + curator-ref participation) RAISE in the
schema model_validators. Conventions mirror tests/test_class2_paired_allele.py
(direct schema construction) and tests/test_pool_reconciliation.py (driving
agent_core.finalize_partial through the .partial.json staging file).
"""
import json
import pathlib
import sys

import pytest
from pydantic import ValidationError

PKT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PKT / "cancervac_packet"))
import schema  # noqa: E402

import agent_core  # noqa: E402

# Reference fixtures (no candidates) — used for back-compat (E).
REF_NAMES = ("rojas", "keskin", "li")


# ---------------------------------------------------------------------------
# Minimal-record builders (one patient + one IMP is a valid ExtractedPaper).
# ---------------------------------------------------------------------------
def _patient(pid="P1"):
    return dict(paper_local_id=pid, indication="melanoma", n_peptides_synthesized=1,
                n_peptides_immunogenic=0, quoted_text="q", section_ref="s")


def _imp(iid="IMP1", seq="NAQVRKCPPVITVNA", patient="P1"):
    return dict(paper_local_id=iid, sequence=seq, is_neoantigen=True,
                patient_paper_id=patient, quoted_text="q", section_ref="s")


def _candidate(cid="C1", seq="NAQVRKCPPVITVNA", status="selected",
               selected="IMP1", scores=None):
    c = dict(paper_local_id=cid, sequence=seq, is_neoantigen=True,
             candidate_status=status, quoted_text="q", section_ref="s",
             ranking_scores=scores if scores is not None else [])
    if selected is not None:
        c["selected_peptide_id"] = selected
    return c


def _paper(*, patients=None, imps=None, candidates=None, evidence=None,
           curator_notes=None, n_predicted_reported=None, n_selected_reported=None):
    rec = dict(pmid="123", journal="J", year=2023, title="T", cohort_size=1,
               indication_summary="mel",
               patients=patients if patients is not None else [_patient()],
               immunizing_peptides=imps if imps is not None else [_imp()])
    if candidates is not None:
        rec["candidates"] = candidates
    if evidence is not None:
        rec["evidence"] = evidence
    if curator_notes is not None:
        rec["curator_notes"] = curator_notes
    if n_predicted_reported is not None:
        rec["n_predicted_reported"] = n_predicted_reported
    if n_selected_reported is not None:
        rec["n_selected_reported"] = n_selected_reported
    return rec


# ===========================================================================
# A. PrioritizationScore — the sibling that keeps non-affinity scores queryable
# ===========================================================================
def test_prioritization_score_keeps_nonaffinity_value_as_float():
    """The whole point of the sibling: a non-affinity score (expression TPM) is a
    queryable FLOAT, not value=None + a lossless string as Measurement would force."""
    s = schema.PrioritizationScore(score_kind="expression_tpm", value=12.5)
    assert s.value == 12.5
    assert isinstance(s.value, float)


def test_prioritization_score_requires_value_or_raw():
    """_lossless: neither a parsed value nor a raw token -> never an empty score."""
    with pytest.raises(ValidationError):
        schema.PrioritizationScore(score_kind="rank")


def test_prioritization_score_raw_only_is_lossless_ok():
    """A figure-only / unparsed score survives losslessly via raw alone."""
    s = schema.PrioritizationScore(score_kind="quality_score", raw="top-decile")
    assert s.value is None and s.raw == "top-decile"


def test_prioritization_score_kind_is_enforced():
    """score_kind is a controlled vocab (lockstep with vocab.SCORE_KINDS)."""
    with pytest.raises(ValidationError):
        schema.PrioritizationScore(score_kind="not_a_real_kind", value=1.0)


# ===========================================================================
# B. NeoantigenCandidate — status<->bridge consistency (model_validator)
# ===========================================================================
def test_predicted_candidate_with_bridge_raises():
    """_selected_link_requires_selected_status: a 'predicted' candidate may not
    carry a selected_peptide_id (only selected/administered may bridge)."""
    with pytest.raises(ValidationError):
        schema.NeoantigenCandidate(
            paper_local_id="C1", sequence="NAQVRKCPPVITVNA", is_neoantigen=True,
            candidate_status="predicted", selected_peptide_id="IMP1",
            quoted_text="q", section_ref="s")


def test_selected_candidate_with_bridge_is_valid():
    """A 'selected' candidate carrying a bridge id is structurally fine at the
    model level (paper-level resolution is checked separately, in C)."""
    c = schema.NeoantigenCandidate(
        paper_local_id="C1", sequence="NAQVRKCPPVITVNA", is_neoantigen=True,
        candidate_status="selected", selected_peptide_id="IMP1",
        ranking_scores=[schema.PrioritizationScore(score_kind="rank", value=1.0)],
        quoted_text="q", section_ref="s")
    assert c.selected_peptide_id == "IMP1"
    assert c.ranking_scores[0].value == 1.0


def test_selected_candidate_resolving_bridge_paper_valid():
    """End-to-end: a selected candidate bridged to a real IMP validates as a paper."""
    rec = _paper(candidates=[_candidate(
        scores=[dict(score_kind="expression_tpm", value=12.5)])])
    p = schema.ExtractedPaper(**rec)
    assert p.candidates[0].selected_peptide_id == "IMP1"


# ===========================================================================
# C. ExtractedPaper candidate guards (HARD — raise)
# ===========================================================================
def test_candidate_id_colliding_with_patient_raises():
    """_cross_reference_check: candidate ids must be disjoint from all other
    entities (here a patient id)."""
    rec = _paper(candidates=[_candidate(cid="P1")])  # collides with the patient
    with pytest.raises(ValidationError, match="collide"):
        schema.ExtractedPaper(**rec)


def test_candidate_id_colliding_with_imp_raises():
    """Same disjointness rule against an immunizing-peptide id."""
    rec = _paper(candidates=[_candidate(cid="IMP1")])  # collides with the IMP
    with pytest.raises(ValidationError, match="collide"):
        schema.ExtractedPaper(**rec)


def test_selected_bridge_must_resolve_to_an_imp():
    """A selected candidate whose selected_peptide_id is not a known IMP raises."""
    rec = _paper(candidates=[_candidate(selected="NO_SUCH_IMP")])
    with pytest.raises(ValidationError, match="does not resolve"):
        schema.ExtractedPaper(**rec)


def test_candidate_target_evidence_resolves_when_candidate_exists():
    """A candidate-target evidence row resolves when the candidate is present."""
    ev = dict(patient_paper_id="P1", target_kind="candidate", candidate_paper_id="C1",
              assay="elispot", outcome="not_immunogenic", quoted_text="q", section_ref="s")
    rec = _paper(candidates=[_candidate(status="predicted", selected=None)], evidence=[ev])
    p = schema.ExtractedPaper(**rec)
    assert p.evidence[0].candidate_paper_id == "C1"


def test_candidate_target_evidence_raises_when_candidate_missing():
    """The same evidence row raises when no such candidate exists."""
    ev = dict(patient_paper_id="P1", target_kind="candidate", candidate_paper_id="GHOST",
              assay="elispot", outcome="not_immunogenic", quoted_text="q", section_ref="s")
    rec = _paper(candidates=[_candidate(status="predicted", selected=None)], evidence=[ev])
    with pytest.raises(ValidationError, match="unknown candidate_paper_id"):
        schema.ExtractedPaper(**rec)


def test_duplicate_candidate_ids_raise():
    """_unique_local_ids: candidate ids participate in within-type uniqueness."""
    rec = _paper(candidates=[
        _candidate(cid="C1", status="predicted", selected=None),
        _candidate(cid="C1", status="predicted", selected=None)])
    with pytest.raises(ValidationError, match="duplicate paper_local_id"):
        schema.ExtractedPaper(**rec)


def test_curator_note_ref_to_candidate_resolves():
    """_curator_refs_resolve: candidate ids are part of the resolvable union."""
    note = dict(kind="highlight", text="strong candidate", refs=["C1"])
    rec = _paper(candidates=[_candidate(status="predicted", selected=None)],
                 curator_notes=[note])
    p = schema.ExtractedPaper(**rec)
    assert p.curator_notes[0].refs == ["C1"]


def test_curator_note_ref_to_unknown_candidate_raises():
    """A curator note referencing a non-existent candidate id raises."""
    note = dict(kind="highlight", text="strong candidate", refs=["GHOST"])
    rec = _paper(candidates=[_candidate(status="predicted", selected=None)],
                 curator_notes=[note])
    with pytest.raises(ValidationError, match="unknown paper_local_id"):
        schema.ExtractedPaper(**rec)


# ===========================================================================
# D. SOFT guards via agent_core.finalize_partial (block-once, overridable)
# ===========================================================================
def _finalize(tmp_path, rec, **kw):
    """Stage the partial and call finalize_partial, mirroring test_pool_reconciliation.
    allow_missing_magnitudes defaults True so unrelated nudges don't shadow ours."""
    out = tmp_path / "r.json"
    (tmp_path / "r.json.partial.json").write_text(json.dumps(rec))
    return agent_core.finalize_partial(str(out), allow_missing_magnitudes=True, **kw)


def test_bridge_seq_mismatch_blocks_then_overrides(tmp_path):
    """#2 SOFT: a selected candidate bridged to an IMP of a DIFFERENT sequence is
    flagged; allow_candidate_bridge_mismatch=True proceeds. n_predicted_reported is
    set so the funnel-size nudge doesn't fire first."""
    rec = _paper(
        candidates=[_candidate(seq="QTQKHLDLYAAAA", selected="IMP1")],  # IMP1 seq differs
        n_predicted_reported=10, n_selected_reported=1)
    out = tmp_path / "r.json"
    part = tmp_path / "r.json.partial.json"

    part.write_text(json.dumps(rec))
    ok, msg = agent_core.finalize_partial(str(out), allow_missing_magnitudes=True)
    assert not ok
    assert "C1" in msg and "sequence" in msg.lower()
    assert part.exists()  # block leaves the partial in place to fix

    part.write_text(json.dumps(rec))  # re-stage
    ok2, msg2 = agent_core.finalize_partial(
        str(out), allow_missing_magnitudes=True, allow_candidate_bridge_mismatch=True)
    assert ok2, msg2
    assert not part.exists()  # success cleans up the partial


def test_bridge_seq_match_does_not_block(tmp_path):
    """A candidate whose sequence MATCHES its bridged IMP passes the #2 nudge."""
    rec = _paper(candidates=[_candidate(seq="NAQVRKCPPVITVNA", selected="IMP1")],
                 n_predicted_reported=10)
    ok, msg = _finalize(tmp_path, rec)
    assert ok, msg


def test_unknown_funnel_size_blocks_then_overrides(tmp_path):
    """#3 SOFT: candidates present but n_predicted_reported unset is nudged;
    allow_unknown_funnel_size=True proceeds. Sequence matches so #2 doesn't fire."""
    rec = _paper(candidates=[_candidate(seq="NAQVRKCPPVITVNA", selected="IMP1")])
    out = tmp_path / "r.json"
    part = tmp_path / "r.json.partial.json"

    part.write_text(json.dumps(rec))
    ok, msg = agent_core.finalize_partial(str(out), allow_missing_magnitudes=True)
    assert not ok
    assert "n_predicted_reported" in msg
    assert part.exists()

    part.write_text(json.dumps(rec))  # re-stage
    ok2, msg2 = agent_core.finalize_partial(
        str(out), allow_missing_magnitudes=True, allow_unknown_funnel_size=True)
    assert ok2, msg2
    assert not part.exists()


def test_funnel_size_known_does_not_block(tmp_path):
    """With n_predicted_reported set, the #3 nudge stays silent."""
    rec = _paper(candidates=[_candidate(seq="NAQVRKCPPVITVNA", selected="IMP1")],
                 n_predicted_reported=322, n_selected_reported=1)
    ok, msg = _finalize(tmp_path, rec)
    assert ok, msg


def test_no_candidates_no_funnel_nudge(tmp_path):
    """A record with no candidates trips neither v2.8 nudge (back-compat path)."""
    rec = _paper()  # no candidates
    ok, msg = _finalize(tmp_path, rec)
    assert ok, msg


# ===========================================================================
# E. BACK-COMPAT — the three reference fixtures (no candidates) still validate
# ===========================================================================
@pytest.mark.parametrize("name", REF_NAMES)
def test_reference_fixture_still_validates(name):
    ref = json.loads((PKT / "reference_records" / f"{name}_extracted.json").read_text())
    p = schema.ExtractedPaper(**ref)
    assert p.candidates == []                 # candidates default to []
    assert p.n_predicted_reported is None
    assert p.n_selected_reported is None
