"""v2.15 safety-grounding guard — the SAE/severity vs CTCAE-grade vs relatedness conflation.

Root-caused on 39041242 (MVX-ONCO-1): the agent set any_grade3plus_related=true from a *serious*
adverse event (SAE = regulatory category, not a grade) and a "severe/life-threatening" descriptor the
paper attributed to DISEASE PROGRESSION, not treatment. A grade>=3 RELATED claim must be grounded in a
verbatim grade>=3 token in `raw`; otherwise the record routes to needs_review.

Modeled on tests/test_finalize_class_ii_gate.py (init -> append -> set_safety_summary -> finalize).
The record carries one patient and no evidence, so ONLY the safety guard can fire.
"""
import json
import pathlib

import agent_core

META = {"pmid": "12345678", "journal": "Test J", "year": 2024, "title": "safety grade guard test",
        "cohort_size": 1, "indication_summary": "melanoma"}

# the bug case: 'serious'/SAE wording + a disease-progression disclaimer, NO grade>=3 token.
UNGROUNDED = {"max_related_grade": 3, "any_grade3plus_related": True, "n_patients_with_related_ae": 3,
              "raw": "One serious adverse event was possibly related; two moderate related events. "
                     "Severe/life-threatening events were all related to disease progression, not treatment."}
# grounded: an explicit grade>=3 attributed to treatment (incl. spelled-out 'grade four').
GROUNDED_DIGIT = {"max_related_grade": 3, "any_grade3plus_related": True,
                  "raw": "Three patients had a grade 3 treatment-related adverse event."}
GROUNDED_WORD = {"max_related_grade": 4, "any_grade3plus_related": True,
                 "raw": "In one case, grade four toxicity with probable relation to treatment was observed."}
LOW_GRADE = {"max_related_grade": 2, "any_grade3plus_related": False,
             "raw": "Related TEAEs were of grade 1 or 2 severity."}


def _assemble(out, safety):
    agent_core.init_partial(str(out), json.dumps(META))
    agent_core.append_section(str(out), "patients", json.dumps(
        [{"quoted_text": "pt 1", "section_ref": "M", "paper_local_id": "P1",
          "indication": "melanoma", "n_peptides_synthesized": 0, "n_peptides_immunogenic": 0}]))
    ok, _ = agent_core.set_safety_summary(str(out), json.dumps(safety))
    assert ok


def test_blocks_on_ungrounded_grade(tmp_path):
    out = tmp_path / "r.json"
    _assemble(out, UNGROUNDED)
    ok, msg = agent_core.finalize_partial(str(out))
    assert not ok
    assert "serious adverse event" in msg.lower() or "grade>=3" in msg
    assert "allow_ungrounded_safety_grade=true" in msg


def test_override_routes_to_needs_review(tmp_path):
    out = tmp_path / "r.json"
    _assemble(out, UNGROUNDED)
    ok, msg = agent_core.finalize_partial(str(out), allow_ungrounded_safety_grade=True)
    assert ok and "OVERRIDES USED" in msg
    rec = json.loads(out.read_text())
    assert rec["finalize_overrides_used"] == ["allow_ungrounded_safety_grade"]
    # HARD override -> the scale lane must NOT treat it as clean
    assert not agent_core.overrides_are_soft_only(["allow_ungrounded_safety_grade"])


def test_grounded_grades_pass(tmp_path):
    for safety in (GROUNDED_DIGIT, GROUNDED_WORD, LOW_GRADE):
        out = tmp_path / f"r_{safety['max_related_grade']}_{safety['any_grade3plus_related']}.json"
        _assemble(out, safety)
        ok, msg = agent_core.finalize_partial(str(out))
        assert ok, f"grounded safety wrongly blocked: {safety['raw']!r} -> {msg}"


def test_unit_guard_direct():
    import types
    import schema

    def wrap(ss):
        return types.SimpleNamespace(safety_summary=schema.SafetySummary(**ss))

    assert agent_core._safety_grade_ungrounded(wrap(UNGROUNDED))          # flags the conflation
    assert agent_core._safety_grade_ungrounded(wrap(GROUNDED_DIGIT)) is None
    assert agent_core._safety_grade_ungrounded(wrap(GROUNDED_WORD)) is None   # spelled-out 'grade four'
    assert agent_core._safety_grade_ungrounded(wrap(LOW_GRADE)) is None
    # severe + treatment-related, no SAE/disclaimer cue -> conservative pass (avoid false positives)
    assert agent_core._safety_grade_ungrounded(
        wrap({"any_grade3plus_related": True,
              "raw": "Two patients had severe treatment-related events."})) is None
