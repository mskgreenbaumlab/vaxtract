import schema
import vocab
from schema_digest import build_schema_digest

DIGEST = build_schema_digest(schema, vocab)
SCHEMA_SRC = (__import__("pathlib").Path(schema.__file__)).read_text()


def test_digest_lists_every_entity_model():
    for name in ("ExtractedPaper", "ExtractedPatient", "ImmunizingPeptide", "MinimalEpitope",
                 "ExtractedEvidence", "ExtractedPeptidePool", "SurvivalOutcome", "Measurement",
                 "ResponseMagnitude", "Provenance"):
        assert name in DIGEST, name


def test_digest_includes_key_field_names():
    for f in ("immunizing_peptides", "epitopes", "predicted_affinity", "mhc_class",
              "survival_outcomes", "vaccine_platform"):
        assert f in DIGEST, f


def test_digest_shows_literal_vocab_values():
    assert "synthetic_long_peptide" in DIGEST  # a VACCINE_PLATFORMS value
    assert "one of[" in DIGEST


def test_digest_is_far_smaller_than_the_source():
    assert len(DIGEST) < 15_000, len(DIGEST)
    assert len(DIGEST) < 0.25 * len(SCHEMA_SRC)


def test_digest_carries_no_validator_code():
    assert "@field_validator" not in DIGEST
    assert "def " not in DIGEST


def test_digest_renders_constraints_cleanly():
    assert "pattern=" in DIGEST          # a regex constraint renders
    assert ("maxlen=" in DIGEST) or ("minlen=" in DIGEST)
    assert ">==" not in DIGEST and "<==" not in DIGEST   # no doubled '='
    assert "{minlen=1, minlen=1}" not in DIGEST          # no duplicate constraint tokens
