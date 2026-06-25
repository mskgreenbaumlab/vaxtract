import json
import pathlib
import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]
REF = json.loads((PKT / "reference_records" / "rojas_extracted.json").read_text())
META = {k: v for k, v in REF.items() if not isinstance(v, list)}
SECTIONS = [s for s in agent_core.SECTION_MODEL if REF.get(s)]  # sections present in the ref


def test_section_model_covers_entity_lists():
    # DERIVED from the schema (not a hardcoded list): every ExtractedPaper field that is a
    # list of an extracted entity (its item model carries a `paper_local_id`) MUST be wired
    # into SECTION_MODEL, or the agent can't write it via add_entities/add_table. This is the
    # invariant that the v2.8 `candidates` gap slipped through — a hardcoded list never
    # mentioned it. curator_notes is the one documented exception (human-authored, its item
    # model has no paper_local_id, so it is naturally excluded by the predicate below).
    import schema as _schema
    from pydantic import BaseModel
    for fname, fi in _schema.ExtractedPaper.model_fields.items():
        ann = fi.annotation
        args = getattr(ann, "__args__", ())
        item = args[0] if args else None
        is_entity_list = (
            getattr(ann, "__origin__", None) is list
            and isinstance(item, type) and issubclass(item, BaseModel)
            and "paper_local_id" in item.model_fields
        )
        if is_entity_list:
            assert fname in agent_core.SECTION_MODEL, (
                f"ExtractedPaper.{fname} is an extracted entity list but is not wired "
                f"into agent_core.SECTION_MODEL — the agent cannot write it"
            )


def test_finalize_tool_forwards_all_override_flags():
    # The MCP `finalize` tool must forward EVERY allow_* override that finalize_partial accepts.
    # When it doesn't, the agent's override is silently dropped and it resorts to clearing data
    # (the v2.8 gap: finalize_partial gained allow_unknown_funnel_size / allow_candidate_bridge_mismatch
    # but the tool forwarded only the older two, so a blocked agent deleted all 103 candidates).
    # Reads extraction_agent.py as text (no SDK import needed).
    import inspect
    params = [p for p in inspect.signature(agent_core.finalize_partial).parameters if p.startswith("allow_")]
    src = (PKT / "vaxtract" / "extraction_agent.py").read_text()
    start = src.index('@tool("finalize"')
    body = src[start:src.index("@tool(", start + 10)]
    for p in params:
        assert p in body, f"finalize_partial override {p!r} is not forwarded by the MCP finalize tool"


def test_can_append_candidates_with_scores_through_tool_path(tmp_path):
    # The funnel must be reachable through the ACTUAL tool the agent uses (add_entities ->
    # append_section), not just constructible in Python. Regression for the v2.8 gap where
    # `candidates` was a valid schema field but an "unknown section" to the tool.
    out = tmp_path / "r.json"
    agent_core.init_partial(str(out), json.dumps(META))
    cand = [{
        "quoted_text": "Tmem101.G96V\tQLASTYTAYIVGYVHYGDWLK", "section_ref": "Table S2",
        "paper_local_id": "CAND_Tmem101", "sequence": "QLASTYTAYIVGYVHYGDWLK",
        "gene_symbol": "Tmem101", "is_neoantigen": True, "candidate_status": "administered",
        "ranking_scores": [
            {"score_kind": "affinity", "value": 109.0, "raw": "MT IC50 109"},
            {"score_kind": "agretopicity_dai", "value": 109.49, "raw": "fold change 109.49"},
            {"score_kind": "expression_tpm", "name": "FPKM", "value": 7.83, "raw": "FPKM 7.83"},
        ],
    }]
    ok, msg = agent_core.append_section(str(out), "candidates", json.dumps(cand))
    assert ok, msg
    part = json.loads((tmp_path / "r.json.partial.json").read_text())
    assert len(part["candidates"]) == 1
    assert part["candidates"][0]["ranking_scores"][0]["value"] == 109.0


def test_init_accepts_good_meta_and_writes_empty_partial(tmp_path):
    out = tmp_path / "r.json"
    ok, msg = agent_core.init_partial(str(out), json.dumps(META))
    assert ok, msg
    part = json.loads((tmp_path / "r.json.partial.json").read_text())
    assert part["immunizing_peptides"] == [] and part["pmid"] == META["pmid"]


def test_init_rejects_incomplete_meta(tmp_path):
    ok, msg = agent_core.init_partial(str(tmp_path / "r.json"), json.dumps({"pmid": "x"}))
    assert not ok and "invalid" in msg.lower()


def test_append_validates_each_item_and_counts(tmp_path):
    out = tmp_path / "r.json"
    agent_core.init_partial(str(out), json.dumps(META))
    items = REF["immunizing_peptides"]
    ok, msg = agent_core.append_section(str(out), "immunizing_peptides", json.dumps(items))
    assert ok and str(len(items)) in msg


def test_append_rejects_a_batch_with_one_bad_item_and_leaves_partial_unchanged(tmp_path):
    out = tmp_path / "r.json"
    agent_core.init_partial(str(out), json.dumps(META))
    bad = [REF["immunizing_peptides"][0], {"nonsense": True}]
    ok, msg = agent_core.append_section(str(out), "immunizing_peptides", json.dumps(bad))
    assert not ok
    part = json.loads((tmp_path / "r.json.partial.json").read_text())
    assert part["immunizing_peptides"] == []  # nothing was appended


def test_clear_section_empties_it(tmp_path):
    out = tmp_path / "r.json"
    agent_core.init_partial(str(out), json.dumps(META))
    agent_core.append_section(str(out), "survival_outcomes", json.dumps(REF["survival_outcomes"]))
    ok, _ = agent_core.clear_section(str(out), "survival_outcomes")
    assert ok
    part = json.loads((tmp_path / "r.json.partial.json").read_text())
    assert part["survival_outcomes"] == []


def test_partial_status_reports_counts(tmp_path):
    out = tmp_path / "r.json"
    agent_core.init_partial(str(out), json.dumps(META))
    agent_core.append_section(str(out), "patients", json.dumps(REF["patients"]))
    ok, msg = agent_core.partial_status(str(out))
    assert ok and "patients" in msg and str(len(REF["patients"])) in msg


def test_finalize_assembles_the_full_record_and_validates(tmp_path):
    out = tmp_path / "r.json"
    agent_core.init_partial(str(out), json.dumps(META))
    for s in SECTIONS:
        ok, msg = agent_core.append_section(str(out), s, json.dumps(REF[s]))
        assert ok, (s, msg)
    ok, msg = agent_core.finalize_partial(str(out))
    assert ok, msg
    assert out.exists()
    ok2, _ = agent_core.validate_record(out.read_text())
    assert ok2
    assert not (tmp_path / "r.json.partial.json").exists()  # partial cleaned up on success


def test_finalize_surfaces_a_cross_entity_error_and_writes_nothing(tmp_path):
    out = tmp_path / "r.json"
    agent_core.init_partial(str(out), json.dumps(META))
    bad_ev = dict(REF["evidence"][0]); bad_ev["patient_paper_id"] = "NO_SUCH_PATIENT"
    agent_core.append_section(str(out), "evidence", json.dumps([bad_ev]))
    ok, msg = agent_core.finalize_partial(str(out))
    assert not ok
    assert not out.exists()  # nothing written on a failed finalize


def test_corrupt_partial_is_reported_not_raised(tmp_path):
    out = tmp_path / "r.json"
    (tmp_path / "r.json.partial.json").write_text("{ this is not valid json")
    for call in (lambda: agent_core.append_section(str(out), "patients", "[]"),
                 lambda: agent_core.clear_section(str(out), "patients"),
                 lambda: agent_core.partial_status(str(out)),
                 lambda: agent_core.finalize_partial(str(out))):
        ok, msg = call()
        assert ok is False and "corrupt" in msg.lower()
