import agent_core
from prompt_render import field_guidance_only

RAW = agent_core.DELTAS_TEXT
TRIM = field_guidance_only(RAW)


def test_keeps_the_field_guidance():
    assert "Two cross-cutting rules" in TRIM
    for d in ("Delta A", "Delta B", "Delta C", "Delta D", "Delta E"):
        assert d in TRIM, d


def test_drops_the_engineer_meta_prose():
    assert "Why this exists" not in TRIM            # version-history preamble
    assert "injection checklist" not in TRIM        # build_gold_outcomes engineer note
    assert "Done-checklist" not in TRIM             # prompt-update checklist
    assert "## v2.7" not in TRIM                     # trailing curator_notes section header


def test_meaningfully_smaller():
    assert len(TRIM) < len(RAW)
    assert len(RAW) - len(TRIM) > 3000              # ~3.7 KB of meta-prose removed


def test_failsafe_returns_original_when_markers_absent():
    # If the doc is restructured so the markers vanish, do NOT silently nuke the guidance.
    assert field_guidance_only("no markers here at all") == "no markers here at all"
