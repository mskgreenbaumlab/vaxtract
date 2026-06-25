import pathlib
import pytest
import vocab
from prompt_render import render_deltas, build_system_prompt

PKT = pathlib.Path(__file__).resolve().parents[1]
DELTAS_TEXT = (PKT / "vaxtract" / "layer2_prompt_deltas.md").read_text()

# The 10 tuples the deltas doc injects via {{join ...}} (its own checklist).
INJECTED = (
    "SPECIES", "VACCINE_PLATFORMS", "EFFICACY_READOUTS", "EFFICACY_RESULTS",
    "EFFICACY_SETTINGS", "COMBINATION_CLASSES", "SURVIVAL_ENDPOINTS",
    "SURVIVAL_TIME_UNITS", "RESPONSE_MAGNITUDE_UNITS", "RESPONSE_GRADES",
)


@pytest.mark.parametrize("name", INJECTED)
def test_each_injected_tuple_is_rendered_in_lockstep(name):
    out = render_deltas(DELTAS_TEXT, vocab)
    expected = ", ".join(getattr(vocab, name))
    assert expected in out                      # the joined vocab values appear
    assert "{{join " + name + "}}" not in out   # the placeholder is gone


def test_no_placeholder_survives():
    out = render_deltas(DELTAS_TEXT, vocab)
    assert "{{" not in out


def test_unknown_placeholder_raises():
    with pytest.raises((AttributeError, ValueError)):
        render_deltas("see {{join NOPE_NOT_A_TUPLE}} here", vocab)


def test_curator_note_kinds_is_never_injected():
    # Guard: curator_notes are human-authored; the deltas doc must not inject them.
    assert "{{join CURATOR_NOTE_KINDS}}" not in DELTAS_TEXT
    out = render_deltas(DELTAS_TEXT, vocab)
    assert ", ".join(vocab.CURATOR_NOTE_KINDS) not in out


def test_build_system_prompt_embeds_rendered_deltas_and_schema():
    schema_src = "SENTINEL_SCHEMA_SOURCE"
    prompt = build_system_prompt(DELTAS_TEXT, schema_src, vocab)
    assert ", ".join(vocab.SPECIES) in prompt   # deltas rendered
    assert "SENTINEL_SCHEMA_SOURCE" in prompt    # schema source embedded
    assert "Never fabricate" in prompt           # the non-negotiable rules
    assert "{{" not in prompt


def test_epitope_linkage_guidance_present():
    prompt = build_system_prompt(DELTAS_TEXT, "X", vocab)
    low = prompt.lower()
    assert "parent_peptide_ids" in low          # epitope→peptide link is called out
    assert "orphan" in low                       # and the failure mode is named


def test_indel_and_peptide_count_guidance_present():
    prompt = build_system_prompt(DELTAS_TEXT, "X", vocab)
    low = prompt.lower()
    assert "indel" in low and "frameshift" in low        # keep indel/frameshift neoantigens
    assert "n_peptides_synthesized" in low               # the count-reconciliation rule


def test_magnitude_finish_step_guidance_present():
    prompt = build_system_prompt(DELTAS_TEXT, "X", vocab)
    low = prompt.lower()
    assert "source-data" in low and "magnitude" in low          # the dig is called out
    assert "allow_missing_magnitudes" in low                    # the block-once override


def test_read_figure_guidance_has_zoom_and_figure_kind():
    prompt = build_system_prompt(DELTAS_TEXT, "X", vocab)
    low = prompt.lower()
    assert "region=" in low                       # two-step zoom
    assert "kind='figure'" in low                 # figure-kind provenance
    assert "needs_review=true" in low             # conservative recording


def test_companion_paper_guidance_present():
    # v2.11.4: secondary-analysis papers that defer their manifest to a prior paper.
    prompt = build_system_prompt(DELTAS_TEXT, "X", vocab)
    low = prompt.lower()
    assert "companion_paper_ref" in low                       # the field is named
    assert "secondary" in low and "earlier work" in low       # the tell-tale signs (defers to prior paper)
    assert "explicitly states" in low and "real miss" in low  # the anti-gaming guard


def test_multisheet_add_table_is_first_choice_guidance():
    # Item #1: sheets=[...] should be the FIRST-choice path, not a per-sheet read_table loop.
    prompt = build_system_prompt(DELTAS_TEXT, "X", vocab)
    low = prompt.lower()
    assert "sheets:[...]" in low or "sheets=[" in low          # the multi-sheet param
    assert "first choice" in low                               # framed as the default
    assert "__sheet__" in low                                  # the reserved patient-deriving column
