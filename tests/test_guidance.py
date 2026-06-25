import agent_core
import prompt_render
from prompt_render import build_system_prompt
from schema_digest import build_schema_digest

PROMPT = build_system_prompt(
    prompt_render.field_guidance_only(agent_core.DELTAS_TEXT),
    build_schema_digest(agent_core.schema, agent_core.vocab),
    agent_core.vocab,
)


def test_evidence_not_one_per_peptide():
    p = PROMPT.lower()
    assert "evidence" in p and "not one per" in p


def test_no_clearing_patients_to_fix_counts():
    p = PROMPT.lower()
    assert "reconcil" in p
    assert "do not clear" in p or "never clear" in p


def test_omit_patient_id_on_bulk_peptides():
    assert "patient_paper_id" in PROMPT


def test_resume_with_partial_status_after_compaction():
    assert "partial_status" in PROMPT and "summari" in PROMPT.lower()


def test_finalize_promptly():
    assert "finalize" in PROMPT
    assert "{{" not in PROMPT
