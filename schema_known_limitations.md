# antVacDB schema — known limitations ("the three remains")

**As of schema `v2.6.0`.** Files: `schema.py` (the contract) + `vocab.py` (controlled
vocabularies, kept in lockstep).

**Status:** Item 1 (response magnitude) is **RESOLVED in v2.6.0** — kept below as the
implemented record. Item 2 (DNA/RNA encoded-count) remains **READY** (not yet done).
Item 3 (raw growth curves) remains **WON'T DO** by design. So one item is genuinely
still open.

This is a **backlog**, not a plan. Item 2 is *ready to implement* if/when you want it.
Item 3 is a **deliberate scope decision — do NOT implement it** unless the revisit
condition is met. Each design below is a **recommended starting point**, not a finished
diff: validate it against the live file before committing (line numbers drift; symbol
names are stable — grep for those).

---

## Read this first — non-negotiable repo invariants

Any change here MUST follow the discipline that the rest of the schema was built on. A
change that skips these is wrong even if it "works":

1. **Vocab lockstep.** Every controlled-vocabulary `Literal` in `schema.py` MUST have a
   matching tuple in `vocab.py` AND a line in the `# VOCAB LOCKSTEP` block
   (`_assert_vocab(MyLiteral, "MY_TUPLE")`, currently ~line 779–797, 19 axes). If they
   disagree, the module raises on import. So: add the `Literal` → add the tuple → add the
   assert. No exceptions.
2. **Work in a copy, validate by running real data.** Edit a working copy, then
   `python3 -m py_compile schema.py` and `import schema` (the import runs all lockstep
   asserts). Then round-trip the real fixtures (below) through `ExtractedPaper(**json)`.
3. **Guard regression — prove no loosening.** After any change, re-run a guard suite that
   shows the *old* guards still fire (off-vocab rejected, name-like `paper_local_id`
   rejected, class-II+nM rejected, count reconciliation still catches a real mismatch,
   etc.). A schema change that silently widens what validates is a regression.
4. **Back-compat.** The existing fixtures MUST still validate unchanged. New fields are
   optional with sensible defaults so a minimal extractor keeps working.
5. **Version bump + changelog.** Bump `SCHEMA_VERSION` (~line 121), add a history line
   (the `# History:` comment) and a `Pxx` note in the module docstring. Per RUNBOOK
   "Loop 3": every schema change → regression + version bump before it's considered done.
   Update the RUNBOOK's `SCHEMA_VERSION` references too.
6. **Layer-2 prompt.** New fields that an extractor should populate also need the Layer-2
   extraction prompt updated (it is generated from this schema). Note it in the changelog.

### Test fixtures (same `outputs/` dir)
- `li_full_v2_4_extracted.json` — 3 cohorts (E0771 mouse, T412/4T1.2 mouse, GTB16 human),
  42 neoantigens, 10 evidence rows, 3 efficacy arms. Primary new-feature target.
- `li_gtb16_extracted.json` — human-only, 1 patient.
- `keskin_extracted_v2_2.json` — 8-patient human SLP trial; the back-compat canary.

### Standard validation snippet
```python
import importlib, vocab, schema
importlib.reload(vocab); importlib.reload(schema)
from schema import ExtractedPaper
import json; from pathlib import Path
for f in ["li_full_v2_4_extracted.json","li_gtb16_extracted.json","keskin_extracted_v2_2.json"]:
    ExtractedPaper(**json.loads((Path("outputs")/f).read_text()))   # must not raise
print(schema.SCHEMA_VERSION, "imports clean; fixtures validate")
```

---

## Item 1 — Immune-response magnitude is free-text  ·  STATUS: ✅ RESOLVED in v2.6.0

> **Shipped in v2.6.0 (P16).** `ResponseMagnitude` is implemented as a sibling of
> `Measurement` (NOT a reuse), attached as `ExtractedEvidence.magnitude`. Lockstep is now
> 23 axes (`RESPONSE_MAGNITUDE_UNITS`, `RESPONSE_GRADES`). Exercised on real data: Rojas
> ELISpot SFC (3 exact `value`+`sfc_per_1e6` for single-response patients; 21 lossless
> `raw` sets for polytopic patients, since Fig 1g labels magnitudes by patient not gene),
> and Li's 7 mouse `+/++/+++` grades migrated out of `assay_detail` into the ordinal lane.
> All prior fixtures still validate; guard regression confirmed the change is additive.
>
> **Deltas from the proposal below (decided during implementation):**
> - Units shipped as `sfc_per_1e6` / `percent_of_parent` / `stimulation_index` / `unknown`
>   — dropped the proposed separate `sfu_per_1e6` (SFC/SFU folded into one canonical unit;
>   `raw` preserves the paper's exact term) and renamed `percent_parent`→`percent_of_parent`,
>   `unspecified`→`unknown` (matches `Measurement`'s lossless-escape naming).
> - Added fields the proposal omitted: `background_subtracted: bool|None`,
>   `denominator: str|None`, `tier: MeasurementTier = "reported"`, `confidence`.
> - Added a minimal cross-field guard on `ExtractedEvidence`
>   (`_magnitude_outcome_consistency`): a positive-outcome row can't carry
>   `grade="negative"` — unambiguous-only, mirroring the `PreclinicalEfficacy` guard.
>
> The original proposal is retained below for provenance.

### Problem
ELISpot magnitude is the headline immunogenicity readout, but it has no structured home.
- Mouse (Tables S2/S3): ordinal grades `-` / `+` / `++` / `+++` / `++++`.
- Mouse/human figures: quantitative IFN-γ **SFU/SFC per 1e6 cells** (e.g. Fig 5, Fig 2).
Today this lives only in `ExtractedEvidence.assay_detail` (free text, `max_length=200`,
~line 842). It is not queryable or comparable across papers — you cannot ask "show
responses > 100 SFC/1e6" or "rank by magnitude."

### Where it lives
`class ExtractedEvidence(_Extracted)` (~line 800). The existing `Measurement` model
(~line 260) is the design template but should **not** be overloaded — its `unit` axis
(`MeasurementUnit` = `nM`/`pct_rank`/`unknown`, ~line 226) is affinity-specific. Adding
ELISpot units there would pollute the affinity axis. Make a sibling type instead.

### Recommended design
A small frozen value object mirroring `Measurement`'s philosophy (value + lossless `raw`
+ controlled unit + provenance), plus an ordinal grade so the `+`-scale survives.

```python
# new controlled axes (need vocab tuples + lockstep lines)
ResponseGrade = Literal["negative", "low", "moderate", "high", "very_high"]
ResponseMagnitudeUnit = Literal[
    "sfc_per_1e6", "sfu_per_1e6", "percent_parent", "stimulation_index", "unspecified",
]

class ResponseMagnitude(_Frozen):
    """Quantitative and/or ordinal magnitude of an immune readout. Lossless `raw`
    keeps the paper's exact token; `grade` carries semiquant scales."""
    value: float | None = Field(default=None, ge=0)     # e.g. 120 (SFC/1e6)
    unit: ResponseMagnitudeUnit = "unspecified"
    grade: ResponseGrade | None = None                  # mapped from -/+/++/+++/++++
    raw: Annotated[str, StringConstraints(max_length=64)] | None = None  # "++++" / "120 SFC/1e6"
    source: Provenance | None = None
```
Add to `ExtractedEvidence`: `magnitude: ResponseMagnitude | None = None`.

**Grade mapping (document for the extractor):** `-`→`negative`, `+`→`low`, `++`→`moderate`,
`+++`→`high`, `++++`→`very_high`. Keep `assay_detail` as-is (don't delete the prose).

**Optional consistency guard** (`@model_validator(mode="after")` on `ResponseMagnitude`):
if `value is not None` then `unit != "unspecified"` (a number needs units). Keep it
minimal — don't force `value`/`grade` to co-occur (a paper may report only one).

### Lockstep impact
2 new Literals → 2 vocab tuples (`RESPONSE_GRADES`, `RESPONSE_MAGNITUDE_UNITS`) → 2
`_assert_vocab` lines. The `ResponseMagnitude` model itself is not vocab.

### Acceptance criteria
- [ ] Imports clean; all (now 21) lockstep axes pass.
- [ ] Re-extract Li magnitudes into the new field: GTB16 → `value`+`sfc_per_1e6` (Fig 5
      SFC counts, background-subtracted); E0771/4T1.2 → `grade` from the Table S2/S3
      `+` symbols. `li_full_v2_4_extracted.json` re-validates with magnitudes attached.
- [ ] Guard regression: off-vocab `grade="huge"` and `unit="spots"` rejected; negative
      `value` rejected; `value` set with `unit="unspecified"` rejected (if guard added).
- [ ] Back-compat: all three fixtures validate with `magnitude` defaulting to `None`.
- [ ] `SCHEMA_VERSION` → `2.5.0`; changelog + RUNBOOK updated; Layer-2 prompt notes the
      new field + grade mapping.

### Notes / risks
- Low risk: purely additive, no existing validator touched.
- Decide whether `percent_parent` (% of CD8+/CD4+, from tetramer/dextramer) belongs here
  or stays separate — tetramer % is arguably a different readout than ELISpot counts.
  Keeping it in the same object with a distinct unit is fine for a first pass.

---

## Item 2 — DNA/RNA vaccine peptide-count semantics  ·  STATUS: READY (more invasive)

### Problem
`n_peptides_synthesized` (required, ~line 643) and `n_peptides_administered`
(~line 648) assume **discrete synthesized peptides**. For a DNA/RNA polyepitope vaccine,
neoantigens are **encoded in a construct**, not synthesized or administered as peptides.
We currently map encoded-count → these fields (e.g. Li GTB16/E0771/4T1.2 = 13/13/16),
which validates but is a documented semantic stretch. The per-paper reconciliation
validator `_peptide_counts_reconcile` (`ExtractedPaper`, ~line 1003) and the per-patient
`_check_cohort_arithmetic` (~line 679) both read these as peptide counts.

### Where it lives
`class ExtractedPatient` fields (~line 610–690) and `ExtractedPaper._peptide_counts_reconcile`
(~line 1003, reconciles per-patient IMP record count against `n_peptides_administered`,
falling back to `synthesized`; skips patients with 0 IMP records).

### Recommended design
Make the count basis explicit and platform-aware.

1. Add to `ExtractedPatient`:
   ```python
   # platform-neutral count of antigens carried by the construct (DNA/RNA vaccines),
   # distinct from peptides synthesized/administered (SLP/peptide vaccines).
   n_neoantigens_encoded: int | None = Field(default=None, ge=0)
   ```
2. **Relax** `n_peptides_synthesized` from required → `int | None = Field(default=None, ge=0)`.
   Rationale: for a DNA vaccine nothing is synthesized as peptide; forcing a number is the
   stretch we're removing. Back-compat: all existing records *provide* it, so they still
   validate — but this is a real loosening of a long-standing invariant, so call it out in
   the changelog and keep `_check_cohort_arithmetic` strict whenever the value IS present.
3. Make reconciliation **platform-aware**. In `_peptide_counts_reconcile`, choose the
   denominator per patient:
   - if `vaccine_platform in {"dna","rna"}` → reconcile IMP count against
     `n_neoantigens_encoded` (fallback to administered/synthesized if encoded is None);
   - else → current logic (`n_peptides_administered` else `synthesized`).
   `_check_cohort_arithmetic`: keep all existing inequalities but guard each on
   `is not None` so a `None` synthesized doesn't crash, and add
   `n_peptides_immunogenic <= n_neoantigens_encoded` when encoded is set.

### Lockstep impact
None **if** you keep this field-only (no new `Literal`). If you instead add a
`CountBasis` enum (`peptide` / `encoded`), that needs the full vocab+lockstep treatment —
not recommended for a first pass; platform already implies the basis.

### Acceptance criteria
- [ ] Imports clean; lockstep unchanged (or +1 axis if you add `CountBasis`).
- [ ] Re-extract Li: set `n_neoantigens_encoded` = 13 (E0771), 16 (4T1.2), 13 (GTB16);
      set `n_peptides_synthesized`/`administered` to `None`. `li_full_v2_4_extracted.json`
      re-validates and reconciliation passes against the encoded count.
- [ ] Guard regression: reconciliation STILL fails on a genuine mismatch (e.g. encoded=13
      but 14 IMP records for that patient); the existing peptide-platform fixtures
      (Keskin) reconcile exactly as before; the 6 historical guard fixtures still fail.
- [ ] Back-compat: Keskin (peptide platform, synthesized set) validates unchanged.
- [ ] `SCHEMA_VERSION` → `2.5.0` (or `2.6.0` if done after Item 1); changelog explicitly
      flags the `n_peptides_synthesized` required→optional loosening; RUNBOOK + Layer-2
      prompt updated.

### Notes / risks
- **Higher risk than Item 1**: it touches two validators (one cross-record) and relaxes a
  required field. Do the guard regression carefully — this is exactly the kind of change
  that can quietly let bad data through.
- Lower-risk alternative if you don't want to relax the required field: add
  `n_neoantigens_encoded` only, prefer it in reconciliation when present, and keep
  `n_peptides_synthesized` required (the DNA stretch persists but the reconciliation is
  honest). Pick based on how much you care about the stretch vs. the invariant.

---

## Item 3 — Raw tumor-growth curves  ·  STATUS: WON'T DO (by design)

### Decision
**Do not implement.** `preclinical_efficacy` (v2.4) deliberately stores the *interpretable
result + arm* (readout → result, by combination/setting), not the per-timepoint
tumor-volume time-series. This mirrors how the schema stores an affinity *value* rather
than a full binding titration: the database captures the comparable conclusion, not the
raw measurement stream.

### Why it's out of scope
- A growth curve is a different data model (per-arm, per-timepoint, per-animal series with
  error bars) that doesn't belong in an entity-level extraction ontology — it would
  dominate the record size and barely compare across papers (different units, schedules,
  endpoints).
- The interpretable signal it carries (did the tumor shrink / was it conditional on ICB)
  is *already captured* in `preclinical_efficacy.result` + `combination`.
- Cross-paper queries want "which vaccines suppressed growth with checkpoint blockade,"
  not raw caliper readings — the categorical field answers that; the curve does not.

### If this is ever revisited
Trigger: a concrete downstream need for quantitative growth dynamics (e.g. modeling TGI%
or survival hazard across the corpus). Even then, the right home is **not** this ontology
— it's an external linked dataset keyed by `pmid` + cohort (figure source-data / a
separate `measurements` table / the deposited repository, e.g. the paper's SRA/dbGaP
accessions), with `preclinical_efficacy` holding a pointer. Keep the ontology categorical.

---

## Suggested order & global checklist

**Order:** Item 1 ✅ shipped in v2.6.0 (additive, low-risk, high analytic payoff —
magnitude queries unlocked). Item 2 is the only one left to do (more invasive; touches
validators). Item 3: leave closed.

Item 1 shipped as its own `2.6.0` (after the v2.5.0 MHC-I/survival work). If Item 2 lands
next it would be `2.7.0` — fine as long as it passes the full loop independently.

**Per-item done-checklist:**
- [ ] `Literal`(s) + `vocab.py` tuple(s) + `_assert_vocab` line(s) all in lockstep.
- [ ] `py_compile` + `import` clean.
- [ ] New feature exercised on `li_full_v2_4_extracted.json`.
- [ ] Guard regression proves no loosening (old guards still fire).
- [ ] All three fixtures still validate (back-compat).
- [ ] `SCHEMA_VERSION` bumped; `# History:` + docstring `Pxx` note added; RUNBOOK
      `SCHEMA_VERSION` references updated.
- [ ] Layer-2 extraction prompt updated for any field an extractor must populate.
