# antVacDB Extraction — Operations Runbook

How `schema.py`, `quarantine.py`, and `slurm_spool.py` fit together in
production, and the three loops you run to keep the pipeline healthy and
improving.

---

## Mental model

**`schema.py` is the contract.** It is the single definition of what a valid
extracted record is, and it refuses anything that doesn't fit. It is imported
anywhere data crosses a boundary, but the boundary that matters is exactly one:
where the LLM's output becomes database rows. Everywhere downstream you work
with validated `ExtractedPaper` objects and trust they are well-formed.

**`quarantine.py` is the operating loop wrapped around that one gate.** It runs
the gate, captures everything that fails (so failures become regression tests
and vocabulary proposals instead of vanishing), and keeps soft statistics on
what passes (so plausible-but-wrong values get flagged). A failure stops being a
dead end and becomes fuel.

**`slurm_spool.py` makes that loop safe under a SLURM array.** The store is a
single-writer resource; producers spool to per-task files, one reducer folds
them in.

---

## Files

| File | Role | Imports |
|---|---|---|
| `schema.py` | Hard contract (Pydantic). The single source of truth shared with the Layer-2 prompt builder. | pydantic, `vocab` |
| `vocab.py` | Controlled vocabularies (assays, outcomes). Lockstep-asserted against the schema at import. | — |
| `quarantine.py` | Quarantine → fixture → prior loop. `QuarantineLoop`, `QuarantineStore`. | `schema` |
| `slurm_spool.py` | SLURM-safe producer/reducer. `Spooler`, `reduce_spools`. | stdlib (producer); `quarantine` (reducer only) |

---

## The ingest boundary

The only place the contract and the loop meet. Validate the **whole paper** as
the primary gate — referential integrity and count reconciliation only run at
that level.

```python
import schema
from quarantine import QuarantineStore, QuarantineLoop

SCHEMA_VERSION = schema.SCHEMA_VERSION          # single source of truth, lives in schema.py
loop = QuarantineLoop(QuarantineStore(DB_PATH, FIXTURES_DIR))   # durable paths, NOT /tmp

def ingest_paper(candidate: dict, pmid: str) -> bool:
    src = {"pmid": pmid, "schema_version": SCHEMA_VERSION, "extractor": MODEL_ID, "run_id": RUN_ID}
    res = loop.ingest(schema.ExtractedPaper, candidate, source=src)
    if not res.ok:
        # optional: re-ingest each row to localize WHICH rows are malformed
        for row in candidate.get("epitopes", []):
            loop.ingest(schema.MinimalEpitope, row, source=src)
        return False                      # nothing partial reaches the DB
    for w in res.warnings:                # soft flags: surface, don't block
        log.warning("%s: %s", pmid, w)
    write_to_db(res.instance)             # only fully-validated data is persisted
    return True
```

Rule: **a partially-valid paper is never written to the DB.** It is quarantined
whole and salvaged on review.

---

## Loop 1 — Continuous ingest

Run `ingest_paper` once per extracted paper (directly, or via the SLURM path
below). Every success quietly sharpens the envelope priors; every failure lands
in the queue. Watch the **quarantine rate** as a health metric — a sudden spike
almost always means upstream breakage (a new journal table format, a model
regression), not a hundred genuinely-bad papers.

---

## Loop 2 — Review the quarantine queue

Periodically (daily, or on a quarantine-rate alert):

```python
for q in loop.store.list(status="open"):      # already sorted by seen_count (highest leverage first)
    ...                                        # inspect q["model"], q["errors"], q["payload"], q["source"]

loop.proposals(min_count=2)                    # recurring novel terms -> vocab promotion candidates
```

For each item, resolve it one of three ways:

```python
loop.store.resolve(qid, status="resolved_fixed",        note="fixed source row; re-ingested")
loop.store.resolve(qid, status="resolved_accept_widen", note="schema too strict; widening")
loop.store.resolve(qid, status="wontfix",               note="malformed source; not a schema bug")
```

Act on `proposals()` output — e.g. when `cytof` has tripped the `assay` enum
enough times, add it to `vocab.ASSAYS` (which forces the matching schema
`Literal` via the lockstep assert).

---

## Loop 3 — Schema change (the antifragile payoff)

When you change `schema.py` (widen a rule, add a vocab term, add a field):

1. Edit `schema.py` / `vocab.py`.
2. Re-run the regression suite assembled from history:
   ```python
   report = loop.replay()                      # re-validates every OPEN fixture vs current schema
   # report["now_pass"]  -> fixtures the change fixed; resolve them
   # report["still_fail"] -> guards that still hold
   ```
   Point the same idea at your **whole historical corpus** (not just fixtures)
   for full shadow-validation: confirm the change fixes what you intended and
   **nothing that used to pass now fails.**
3. Close the fixtures that now pass; bump `SCHEMA_VERSION` in `schema.py`.

This is the loop that makes the system get *stronger* with each stressor: a
failure mode, once fixed, can never silently regress.

**Worked example — v2.5.0 (the Rojas stressor).** The autogene-cevumeran PDAC
paper (Rojas et al. 2023) exposed two gaps: predicted MHC-I epitopes with no
printed IC50 (the class-I validator demanded a *parsed* affinity, so 232 real
epitope+allele claims fell out) and a recurrence-free-survival primary endpoint
(RANO `clinical_outcome` cannot express it). Both went through this loop:
(1) relax `_class_affinity_consistency` to require the `predicted_affinity`
*slot* but allow a value-less `unit="unknown"` Measurement when the source omits
the number, and add a paper-level `SurvivalOutcome` model + two vocab axes
(`SURVIVAL_ENDPOINTS`, `SURVIVAL_TIME_UNITS`); (2) guard regression proved the
relaxation was **scoped** — slot still mandatory, nM-on-class-II still rejected,
the new SurvivalOutcome guards fire — and that all five prior fixtures still
validate unchanged; (3) bump `SCHEMA_VERSION` → 2.5.0. The 232 MHC-I epitopes
and the RFS/OS/HR result that previously fell out of the contract are now
first-class. A textbook *widen-without-loosening*: nothing that used to pass
now fails, and the only new things that pass are exactly the two intended cases.

---

## SLURM deployment (spool → reducer)

The store is single-writer. **Never point N concurrent array tasks at it.**
Producers spool to per-task files (no contention); one reducer folds them in.

**Producer — inside each array task** (`extract_task.py`):

```python
from slurm_spool import Spooler                # stdlib only; no schema/store import
sp = Spooler(SPOOL_ROOT)                       # run_id/task_id auto-read from SLURM env vars

candidate = run_llm_extraction(paper)          # -> dict
sp.spool("ExtractedPaper", candidate, source={"pmid": pmid, "extractor": MODEL_ID})
```

**Reducer — one process after the array completes** (`reduce.py`):

```python
from quarantine import QuarantineStore, QuarantineLoop
from slurm_spool import reduce_spools

loop = QuarantineLoop(QuarantineStore(DB_PATH, FIXTURES_DIR))   # the SOLE writer
report = reduce_spools(SPOOL_ROOT, run_id=JOB_ID, loop=loop)
print(report)                                   # accepted / quarantined / warnings
```

**sbatch sketch:**

```bash
# 1) fan out extraction across papers
sbatch --array=0-199 --job-name=extract extract_task.sbatch
# 2) reduce ONCE, gated on the array finishing
sbatch --dependency=afterany:$EXTRACT_JOBID --job-name=reduce reduce.sbatch
```

Properties you get for free: producers can't corrupt the store; the reducer
takes an exclusive `flock` so a second reducer is refused; completed spool files
are atomically moved to `consumed/`, so re-running the reducer is a safe no-op;
quarantine capture is exactly-once (content-hash dedup). (Soft priors are
at-least-once — a reducer crash mid-file may double-count a few statistical
observations on re-run; harmless for an envelope.)

---

## Production checklist

- [ ] `DB_PATH` and `FIXTURES_DIR` are **durable** paths, not `/tmp`.
- [ ] `FIXTURES_DIR` is **under version control** — the regression corpus shows up in PRs.
- [ ] Every record carries `schema_version` + `extractor` id + `run_id` in `source` (drift forensics).
- [ ] Envelope priors are fed only from **reviewed / high-confidence** accepts — gate `envelope_observe` on `needs_review`/`confidence` so a systematic schema-valid bug can't poison "what normal looks like."
- [ ] Exactly **one reducer** per run; producers never import the store.
- [ ] Quarantine rate is monitored/alerted.
- [ ] Store + fixtures are backed up — they are institutional memory now.
- [ ] `vocab.py` stays the single source shared with the prompt builder (lockstep assert holds).
- [ ] On every `schema.py` change: `replay()` + corpus shadow-validation before deploy.

---

## Pitfalls

- **Partial writes.** Never DB-write a paper that only partially validated.
- **Silent schema edits.** Never change `schema.py` without `replay()`.
- **Prompt/schema drift.** Keep the Layer-2 prompt aligned to the schema; the `vocab` lockstep assert guards the enums, extend that discipline to new fields.
- **Concurrent writers.** One SQLite/DuckDB file + N writers = corruption. Use the spool/reducer.
- **Quarantine as a dumping ground.** The queue *is* the roadmap for where the schema and extractor need to improve — read it.
- **Poisoned priors.** Don't learn envelopes from unreviewed data.

---

## API quick reference

```python
# contract
schema.ExtractedPaper(**payload)                  # raises ValidationError if invalid
schema.Measurement.from_raw("0.80 (%rank)")       # lossless parse: value/unit + raw kept

# loop
loop.ingest(model_cls, payload, source) -> IngestResult(ok, instance, quarantine, warnings)
loop.replay() -> {"now_pass": [...], "still_fail": [...]}
loop.proposals(min_count=2) -> {"promote_vocab": [...]}

# store
store.list(status="open")                         # queue, sorted by seen_count desc
store.resolve(qid, status, note)
store.recurring_unknowns(min_count=2)

# slurm
Spooler(SPOOL_ROOT).spool(model_name, payload, source)   # producer (array task)
reduce_spools(SPOOL_ROOT, run_id, loop)                  # reducer (one process)
```

---

## Schema version history

The contract's version is **not hardcoded here** — it is read from
`schema.SCHEMA_VERSION` at the ingest boundary (above) and stamped onto every
record's `source`. This table mirrors the changelog in `schema.py`'s module
docstring; each bump went through Loop 3 (guard regression + no-loosening proof
+ all prior fixtures re-validated). Current: **v2.7.0**.

| Version | Change | Motivating stressor |
|---|---|---|
| 2.0 | Split single `peptides` into `immunizing_peptides` + `epitopes`; referential-integrity + per-patient count-reconciliation validators. | Keskin two-level pool/epitope structure |
| 2.1 | `Measurement` (unit-explicit, lossless `raw`, tier/method/source) replaces bare nM floats; affinities become provenanced claims. | Keskin Supp Table 5 mixing nM and %rank |
| 2.2 | `ConcomitantImmunosuppression` — an immunogenicity confounder recorded as fact, never conclusion. | Keskin dexamethasone during priming |
| 2.3 | `species` + general MHC nomenclature (murine H-2, `MhcAllele`); `vaccine_platform`. | Li et al. 2021 mouse DNA-vaccine cohorts |
| 2.4 | `PreclinicalEfficacy` — tumour-growth/survival by treatment arm; the preclinical analog of RANO `clinical_outcome`. | Li E0771/4T1.2 efficacy-by-arm |
| 2.4.1 | `readout`↔`result` consistency guard on `PreclinicalEfficacy`. | (hardening) survival-result-on-caliper-readout slips |
| **2.5.0** | **P14** — class-I epitope affinity rule relaxed to honor the lossless `Measurement` lane: a class-I epitope must carry a `predicted_affinity` *slot*, but its value may be absent (`unit="unknown"` + `raw`) when the source reports the epitope + allele without a number. **P15** — `SurvivalOutcome`, a paper-level human time-to-event outcome (RFS/OS/HR), the human analog of `PreclinicalEfficacy`; +2 vocab axes (`SURVIVAL_ENDPOINTS`, `SURVIVAL_TIME_UNITS`); lockstep now guards 21 axes. | Rojas et al. 2023 (autogene cevumeran PDAC): 232 predicted MHC-I epitopes without IC50; recurrence-free-survival primary endpoint |
| **2.6.0** | **P16** — `ResponseMagnitude`, a structured home for T-cell response size on every evidence row (parsed `value`+`unit` for SFC/1e6 · % of parent · stimulation index; an ordinal `grade` for `-/+/++/+++/++++`; lossless `raw`). A **sibling** of `Measurement`, deliberately not a reuse — magnitude units stay un-confusable from affinity units and an ordinal grade has a home a float `value` cannot give it. +2 vocab axes (`RESPONSE_MAGNITUDE_UNITS`, `RESPONSE_GRADES`); lockstep now guards 23 axes. | Rojas Fig 1g per-neoantigen SFC + Li mouse `+/++/+++` grades, both previously lost to `assay_detail` free text |
| **2.7.0** | **P17** — `CuratorNote`, an OPTIONAL paper-level container for curatorial commentary (`challenge`/`decision`/`caveat`/`highlight`). A **passive** layer ON TOP of the facts: the report prints it verbatim (no model call at render time), `needs_review` defaults **True** (inverted), and any `refs` MUST resolve to known `paper_local_id`s. +1 vocab axis (`CURATOR_NOTE_KINDS`); lockstep now guards 24 axes. Every existing record ships with it empty. | Recurring desire for the hand-built "challenges & decisions" ledger as a repeatable, hallucination-safe feature rather than free-typed prose |

**Open backlog** (tracked in `schema_known_limitations.md`, not yet implemented):
DNA/RNA encoded-count semantics (proposed `n_neoantigens_encoded` + platform-aware
reconciliation); raw tumour-growth curves are out of scope by design. *(ELISpot
response magnitude — resolved in v2.6.0.)*
