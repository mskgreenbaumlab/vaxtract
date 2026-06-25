"""v2.12 P23: finalize wiring for the class-II minting nudge — block-once, override audit, tool exposure.

Modeled on tests/test_recall_guards.py (init_partial -> append_section -> finalize_partial). The record
carries a CD4 cue (in an evidence quoted_text) and ZERO class-II records, so _class_ii_minting_gap fires;
the evidence row carries a raw magnitude so the magnitude guard does NOT also fire (keeps the override
list to exactly the class-II one).
"""
import json
import pathlib

import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]
META = {"pmid": "12345678", "journal": "Test J", "year": 2024, "title": "class-II gate test",
        "cohort_size": 1, "indication_summary": "melanoma"}


def _patient(pid, nsyn):
    return {"quoted_text": f"pt {pid}", "section_ref": "M", "paper_local_id": pid,
            "indication": "melanoma", "n_peptides_synthesized": nsyn, "n_peptides_immunogenic": 0}

def _imp(pid):
    return {"paper_local_id": pid, "sequence": "SLLQHLIGL", "is_neoantigen": True,
            "quoted_text": "q", "section_ref": "s"}

def _cd4_evidence():
    # CD4 cue in quoted_text, no class-II epitope minted -> _class_ii_minting_gap fires.
    # raw magnitude present -> the magnitude guard does NOT fire (single override expected).
    return {"patient_paper_id": "P1", "target_kind": "immunizing_peptide",
            "immunizing_peptide_paper_id": "i0", "assay": "elispot", "outcome": "immunogenic",
            "quoted_text": "CD4+ helper response to the long peptide", "section_ref": "Results",
            "magnitude": {"raw": "120 SFC/1e6"}}

def _assemble(out):
    agent_core.init_partial(str(out), json.dumps(META))
    agent_core.append_section(str(out), "patients", json.dumps([_patient("P1", 2)]))
    agent_core.append_section(str(out), "immunizing_peptides", json.dumps([_imp("i0"), _imp("i1")]))
    agent_core.append_section(str(out), "evidence", json.dumps([_cd4_evidence()]))


def test_finalize_blocks_on_class_ii_gap(tmp_path):
    out = tmp_path / "r.json"
    _assemble(out)
    ok, msg = agent_core.finalize_partial(str(out))
    assert not ok and "class-II coverage" in msg and "allow_missing_class_ii=true" in msg

def test_finalize_override_recorded(tmp_path):
    out = tmp_path / "r.json"
    _assemble(out)
    ok, msg = agent_core.finalize_partial(str(out), allow_missing_class_ii=True)
    assert ok and "OVERRIDES USED" in msg
    rec = json.loads(out.read_text())
    assert rec["finalize_overrides_used"] == ["allow_missing_class_ii"]

def test_tool_schema_exposes_override():
    # the finalize tool must expose the override in its JSON schema AND forward it in the handler
    src = (PKT / "vaxtract" / "extraction_agent.py").read_text()
    body = src[src.index('@tool("finalize"'):src.index("return", src.index("async def finalize"))]
    assert '"allow_missing_class_ii"' in body                       # schema property
    assert 'allow_missing_class_ii=bool(args.get("allow_missing_class_ii"))' in body  # handler forward
