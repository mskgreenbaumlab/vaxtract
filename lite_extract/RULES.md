# Cross-field rules (the part the schema digest leaves out)

These are the schema's `model_validator`/`field_validator` rules that a field-by-field digest
does not convey. In the 2026-06-11 A/B, an uncoached native agent produced great recall but
INVALID records, failing almost entirely on these. Stating them explicitly is most of what the
heavy MCP harness buys.

## Evidence (the #1 failure)
- Every `ExtractedEvidence` has a `target_kind` ∈ {peptide, epitope, pool, candidate}. Set
  **EXACTLY ONE** of `peptide_paper_id` / `epitope_paper_id` / `pool_paper_id` / `candidate_paper_id`
  — the one matching `target_kind` — and leave the **other three unset**. (`{target_kind}-target
  evidence needs exactly its own *_paper_id set and the other unset`.)
- `evidence outcome` and `magnitude.grade` must agree (don't mark `not_immunogenic` with a positive grade).
- Record reported **negatives** (`outcome=not_immunogenic`, `vaccine_induced=False`) at the grain the
  paper names them — they are completeness, not noise.
- **MANIFEST SCREENING READOUT goes in `screening_readouts`, NOT `evidence` (B1/#4).** A prediction/target
  MANIFEST with a per-row Response/"No response" column (e.g. Rojas 37165196's 232-target ELISpot table:
  200 "No response") is the **screening denominator**. Store EVERY such row as a `ScreeningReadout`
  (`manifest_outcome` = response|no_response|not_evaluable, bridged to its candidate/peptide) — the facts
  are KEPT, just off the immunogenicity axis. Do **NOT** emit `not_immunogenic` `evidence` rows for them.
  - **Discriminator:** a BULK per-target predicted-then-tested column (one row per screened target, no
    per-patient measured magnitude) = `screening_readouts`. A per-patient MEASURED response (an ELISpot
    SFC for a specific patient × peptide, a named responder/non-responder in Results) = `evidence`.
  - **Set `patient_paper_id`** when the manifest is per-patient (one row per patient × neoantigen, e.g.
    Rojas's `Patient number` column) — leave it `null` ONLY for a genuinely cohort-level manifest.
  - `evidence` still gets the POSITIVES and any negative the paper NAMES in Results/figures. A target that
    has a named `evidence` row must NOT also appear in `screening_readouts` (finalize drops the dup).
  - So Rojas = ~25 `evidence` rows (positives + named negatives) AND ~200 `screening_readouts` (the
    enumerated "No response" denominator) — both kept, neither conflated.
- A **pooled** readout = ONE `pool`-target evidence row + one `ExtractedPeptidePool`; emit per-peptide
  rows only when the source deconvolutes the pool. Figure-only/undeconvoluted ⇒ pool grain + `needs_review`.

## Referential integrity (no orphans, no dangling refs)
- Every `*_paper_id` reference must resolve to an entity that exists in the record:
  evidence→patient/peptide/epitope/pool/candidate; pool→patient + member peptides; epitope→parent peptide(s).
- No epitope with empty `parent_peptide_ids`.
- `paper_local_id` is a stable label, **unique within its entity type**, **not** patient initials, and must
  **not collide across entity types** (candidate ids can't reuse a peptide/epitope id).

## Epitope MHC class (MINT + TYPE — affirmative)
- **GRAIN — do not confuse an epitope with a long peptide.** A `MinimalEpitope` is the MINIMAL binder:
  a class-I epitope is an **8–11mer** (the corpus class-I epitopes are 8–12 aa), a class-II core is
  ~13–25mer. The long synthetic/vaccine peptide (15–30+mer) is an **`immunizing_peptide`**, NOT an
  `epitope`. If the paper predicts a minimal epitope per peptide (a NetMHC IC50/%rank table), MINT those
  minimal epitopes at minimal grain and link them to their parent peptide — do **not** re-list the long
  peptide in `epitopes`, and do **not** drop the minimal-epitope layer (a class-I "epitope" longer than
  14 aa or a peptide-table with zero minimal epitopes is a grain error that routes to needs_review).
- Every `MinimalEpitope.mhc_class` is MANDATORY and is `"I"` or `"II"`.
- **DO mint class-II epitopes** when the paper restricts a defined sequence to class II: a named
  **DR/DP/DQ** allele (human) or **I-A/I-E** (mouse) on `hla_allele`, OR a class-II tie in the prose
  (`class II–restricted`, tetramer/transfectant/blocking-antibody evidence) bound to that sequence.
  A class-II epitope with NEITHER a named allele NOR a class-II cue in `quoted_text` is REJECTED — do
  not emit a `"II"` label you cannot quote.
- **Class-I** epitope must carry a `predicted_affinity` slot. If the source gives no IC50, attach a
  lossless `Measurement(unit='unknown')` rather than omitting it. Class-I cannot carry a class-II-only field.
- **Class-II** epitope must **not** carry an nM affinity.
- Class-II HLA restriction may be a paired α/β heterodimer (`HLA-DPA1*01:03/DPB1*02:01`) — keep the pair.
- Mouse class-II alleles (`I-Ab`, `I-Ed`, `I-Ag7`, `H-2IAb`) are valid `hla_allele` values.

## Evidence T-cell subset + MHC class (affirmative)
- On EACH `ExtractedEvidence` you MAY set two ORTHOGONAL optional axes:
  - `t_cell_subset` ∈ {`cd4`, `cd8`, `bulk_or_unknown`}. Set `cd4`/`cd8` ONLY when the row's source
    states the compartment (`CD8+`, `cytotoxic`, `CD4+`, `helper`). Otherwise `bulk_or_unknown` (e.g. a
    bulk PBMC ELISpot). NEVER guess — a non-unknown value MUST have a verbatim cue in `quoted_text` /
    `assay_detail` / `hla_allele` or the record is REJECTED.
  - `mhc_class` ∈ {`class_i`, `class_ii`, `not_determined`}. Same rule: set `class_i`/`class_ii` only
    with a quoted restriction/class cue (named HLA-A/B/C or DR/DP/DQ, "class I/II-restricted", CD8/CD4).
- **One long peptide recognized by BOTH CD4 and CD8 = TWO evidence rows** (one `t_cell_subset='cd4'`,
  one `t_cell_subset='cd8'`), never one. Do NOT put subset on the peptide.
- Leave both fields unset (`null`) only if you did not process the subset/class question for that row.

## Measurements / magnitudes
- A **parsed** value requires a **known unit** (not `'unknown'`); `unit='unknown'` is allowed only for a
  raw/lossless token with no parsed number.
- Ranges: nM ≥ 0; pct_rank ∈ [0,100]; magnitude ≥ 0; percent_of_parent ∈ [0,100].

## NeoantigenMutation (genomic grain — NOT peptide-shaped)
- Real fields: `gene_symbol`, `genomic_change`, `hla_restrictions`, `status`
  (emerged|lost|retained|present), `clonality`, `cancer_cell_fraction`, per-timepoint VAF, `peptide_ref`.
  Needs **at least** `gene_symbol` OR `genomic_change`. Do **not** invent fields (extra inputs are rejected).
- Use this (not `immunizing_peptides`) for sequence-less named gene+protein-change targets, e.g. a
  companion/manifest-deferred paper's Supp-Table gene list.

## Patient / cohort / regimen
- **`cohort_kind` — TAG methodological arms, never DROP them (A1).** Each `ExtractedPatient` may set
  `cohort_kind` ∈ {`patient`, `tumor_model`, `model_antigen_validation`, `healthy_donor`, `other`}. A
  **methodological / pipeline-validation arm** — e.g. an HLA-A2-transgenic set using MODEL antigens
  (viral, gp100) to validate epitope prediction (Li 33879241 P1), NOT a disease-vaccination cohort — must
  be EXTRACTED and tagged `model_antigen_validation`; do not omit it. Real tumour-bearing animal cohorts
  (E0771, 4T1.2) = `tumor_model`; human disease cohorts = `patient`. Leave `null` only if genuinely
  un-determinable. (The "vaccinated patient" recall denominator excludes the methodological arms.)
- `adjuvant='none'` cannot also carry `adjuvant_detail`.
- `len(patients)` must equal the per-patient roster the record implies; `cohort_size` ≥ 1.
- **Peptide-count reconciliation:** `sum(patient.n_peptides_synthesized) == len(immunizing_peptides)`.
- `n_doses_received` ≤ planned priming+boost; `n_peptides_immunogenic` ≤ `n_peptides_administered`.
- **Concomitant therapy (axis #1) — combination is near-universal, DO capture it.** For each
  co-administered systemic anti-cancer agent, add a `ConcomitantTherapy` to
  `patient.concomitant_therapy`: `drug_class` ∈ {`checkpoint_inhibitor`,`chemotherapy`,`radiotherapy`,
  `targeted`,`other`} (the comparable axis), `agent` (verbatim), `timing` ∈ {`concurrent`,`sequential`,
  `unknown`}, optional `line`. Rojas = atezolizumab (checkpoint_inhibitor) + mFOLFIRINOX (chemotherapy).
  This is **distinct from `immunosuppression`** (an immunogenicity confounder). Never infer an agent.
- **Safety (axis #2).** Set the paper-level `safety_summary` headline facts: `max_related_grade`
  (CTCAE 1-5 — **distinct from the immunogenicity `grade`**), `any_grade3plus_related`,
  `n_patients_with_related_ae`, `irae_present`, `raw` (verbatim). Full per-event AE tables are deferred.
- **TMB / MSI (axis #3).** Set `patient.tmb_value` + `patient.tmb_raw` (verbatim, e.g. '12 mut/Mb') and
  `patient.msi_status` ∈ {`mss`,`msi_high`,`unknown`}. `null` = not reported — never infer.

## Survival / benefit
- A survival readout cannot carry a size/burden result; `pre_existing` outcome can't be `vaccine_induced=True`.

## Companion / manifest-absent papers
- If the per-neoantigen sequences are deferred to a prior publication, set `companion_paper_ref` and
  capture what IS reported (named gene/mutation targets → `neoantigen_mutations`; per-patient response
  counts and survival). Do **NOT** fabricate peptide sequences.
