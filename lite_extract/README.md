# lite_extract — prototype of the lighter extraction architecture

**Author:** Samuel Ahuno · **Date:** 2026-06-11 · **Status:** prototype (validated on 38538867; see below)

A response to the 2026-06-11 vanilla-vs-harness A/B (`docs/batch_extraction_report_2026-06-11.md` §7).
That experiment showed an uncoached native agent **out-recalled** the MCP harness (it recovered
image-locked peptide tables the harness's xlsx/table tools are blind to) but produced **schema-invalid**
records. The harness's irreducible value distilled to two things — **cross-field rule enforcement** and
**a deterministic validation/finalize layer** — both far lighter than the current ~5-tool, multi-page
MCP surface.

## The architecture

| Heavy harness (`extraction_agent.py`) | lite_extract |
|---|---|
| ~5 custom MCP tools (read_table, add_table, build_pools, survey_sources, …) | **native tools only** (Read, Bash, Grep, Write) |
| add_entities enforces rules at insert | **explicit `RULES.md`** in the prompt |
| read_table / add_table on xlsx (image-blind) | **Read renders PDF pages as images — first-class** (R-A1) |
| finalize_partial + soft/hard guards | **validate against the real `schema.py` via Bash python, repair-loop** |
| pages of guidance | schema digest + `RULES.md` |

Three pieces:
1. **`RULES.md`** — the cross-field `model_validator` rules a field-by-field digest omits (the evidence
   `target_kind` ownership rule, class-I `predicted_affinity`, NeoantigenMutation shape, referential
   integrity, peptide-count reconciliation, companion-deferred handling, …). This is most of what the
   heavy harness's tools silently carried.
2. **`build_prompt.py`** — assembles the full extraction prompt: image-first reading guidance + the
   schema digest + `RULES.md` + a validate-and-repair finalize block.
3. **The finalize loop** — the agent validates its own JSON against `agent_core.validate_record`
   (the *real* pydantic schema) via Bash python and repairs until it prints `PASS`. The schema IS the
   contract; no MCP tool needed to enforce it. (The A/B's Rojas agent did exactly this and produced a
   valid 16/232/451/225/2 record.)

## Run it

```bash
# 1) make sure the schema digest exists (regenerate if schema bumped):
python - <<'PY'
import sys; sys.path[:0]=['.','cancervac_packet']
import schema_digest, schema
open('outputs/vanilla_compare_2026-06-11/SCHEMA_DIGEST.txt','w').write(
    schema_digest.build_schema_digest(schema))
PY

# 2) build the prompt for a paper and hand it to a native-tool agent
python lite_extract/build_prompt.py --pmid 38538867 \
    --paper-dir /abs/path/data/raw/pubmed/38538867 \
    --out /abs/path/outputs/lite/38538867.json
# -> feed stdout to: claude -p "<prompt>"  (native tools), or the Agent tool, or any SDK shell with
#    default Read/Bash/Grep/Write and NO custom MCP config.
```

## Validation status
- **38538867** (harness got 19/**0/0/0** — image-locked manifest): lite_extract result recorded in
  `outputs/lite/` + the report §7 addendum. This is the headline test of R-A1 (image reading) + R-A3
  (rules + repair → valid).
- Open: run across the full corpus and compare clean-rate + recall + cost vs the heavy harness before
  deciding whether lite_extract replaces or augments it.
