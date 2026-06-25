# Reproducibility

`vaxtract` is an LLM agent, so "reproducible" has a specific, bounded meaning here. This document
states the guarantees and their limits honestly, for reviewers and downstream users.

## What is deterministic

- **The schema and validation.** `vaxtract.schema` (`SCHEMA_VERSION`) is a fixed, versioned Pydantic
  contract. Given the same JSON, validation is fully deterministic and reproducible across machines.
- **The pure logic.** `agent_core`, `prompt_render`, `schema_digest`, the deterministic record
  builders, and the `eval/` scorer are ordinary unit-tested Python with no model calls.
- **The test suite.** `pytest` runs offline and deterministically — it never calls the model.

## What is NOT byte-for-byte reproducible (and why)

- **Model outputs.** The agent calls a proprietary, versioned, non-deterministic model
  (`claude-opus-4-8[1m]`). Two runs on the same paper can differ in wording, ordering, and
  occasionally in borderline judgment calls; a future model version may differ more. We pin the
  model id in `extraction_agent.py` and record it, but we do **not** claim bit-identical re-runs.
- **The source corpus is copyright-restricted.** The papers the agent was validated on cannot be
  redistributed, so a third party cannot re-run the exact extractions from the original PDFs.

## The reproducibility anchor: `reference_records/`

Because of the two limits above, the unit of reproducibility is the **audited gold extraction**, not
a live re-run. `reference_records/` contains human-audited extractions; `eval/` scores any new
extraction's precision/recall against them. To reproduce a reported result:

1. Install a pinned version: `pip install "vaxtract[agent,figures]==<version>"`.
2. Validate the shipped gold records against the schema (no model needed):
   `python cancervac_packet/predeploy_gate.py reference_records`.
3. (Optional, needs a key + the source paper) re-run the agent and score against gold with `eval/`.

## Provenance

Every extracted fact carries `quoted_text` / source anchors, so a curator can verify each value
against the paper independently of how it was produced. Conservatively-gated records land in
`needs_review` rather than being emitted as clean — see [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) §6.

## Environment

- Python ≥ 3.10; dependencies pinned via the `[agent]` / `[figures]` extras in `pyproject.toml`.
- Running the agent additionally requires the Claude Code CLI on `PATH` (the Docker image bundles a
  fixed version).
