"""Pydantic schema for LLM-extracted paper records.

This module is **the contract** between (LLM-based extraction) and Phase D (validation). Any LLM call must produce an
`ExtractedPaper` instance — anything that doesn't fit gets quarantined,
not coerced.

Design rules:
  - Strict types only. No `Any`. No untyped dicts.
  - Every field that points at a peptide / HLA / mutation gets a
    canonical-form regex so downstream `antvac.ids` calls are guaranteed
    to succeed.
  - Cohort arithmetic is enforced at the model level
    (`n_immunogenic <= n_tested`).
  - All extracted rows carry a `quoted_text` + `section_ref` pair so
    Phase D can round-trip them back to the source.
  - `model_config = ConfigDict(frozen=True)` — extracted records are
    immutable post-extraction; downstream code transforms them into
    DB rows but never mutates in place.

PubMed metadata (PMID, DOI, PMCID, NCT, journal, year) is fetched
deterministically by `pubmed.fetch`; the LLM never invents these fields.

===========================================================================
v2 PATCH — two-level peptides, MHC class, structured provenance/affinity.
Motivated by the Keskin (PMID 30568305) extraction, which exposed gaps the
single-`ExtractedPeptide` model could not represent honestly:

  P1. TWO peptide levels. The injected unit is a *long immunizing peptide*
      (what a pool contains); the predicted/known *minimal epitope* lives
      WITHIN it, carries the class-I affinity, and one long peptide yields
      several epitopes (tiling -> many-to-many). `ExtractedPeptide` is split
      into `ImmunizingPeptide` (pool member) and `MinimalEpitope`
      (`parent_peptide_ids` -> IMPs). `ExtractedPaper.peptides` becomes
      `immunizing_peptides` + `epitopes`. (Keskin: 103 IMPs, 105 epitopes.)
  P2. MHC CLASS is a first-class axis. NetMHCpan predicts class I only; the
      Keskin CD4 hits (GPC1, SHANK2, SVEP1) are class-II-restricted and
      carry NO class-I nM affinity. `MinimalEpitope.mhc_class` + a validator
      stop class-II epitopes from being given a class-I affinity.
  P3. AFFINITY units are heterogeneous. Supp Table 5 mixes nM and %rank in
      one column (ARHGAP35, a validated hit, is reported as %rank). v2.1
      models affinities as `Measurement` (value + explicit `unit` + lossless
      `raw` + method/tier/source) so %rank can never masquerade as nM, and the
      verbatim source token survives for retroactive re-parsing.
  P4. EVIDENCE can target either peptide level. `EvidenceTarget` gains
      `immunizing_peptide` and `epitope` (CD4 -> long peptide; CD8/class-I
      -> epitope), resolving the Keskin evidence-linkage gap.
  P5. MUTATION-SPECIFICITY. COX18 was recognized equally for mutant and WT
      (reactive but NOT a neoantigen response). `ExtractedEvidence.
      mutation_specific` captures it without overloading `outcome`.
  P6. ASSESSED vs NEGATIVE. Patients not tested for immunogenicity (Keskin
      Pt 1-3, no booster) are not the same as tested-negative (Pt 4-6).
      `ExtractedPatient.immunogenicity_assessed` separates them.
  P7. STRUCTURED, MULTI-SOURCE PROVENANCE + CONFIDENCE. A fact can have
      several sources of different kinds (dexamethasone: Supp Table 4a AND
      prose). `Provenance` + `confidence`/`needs_review` augment the legacy
      `quoted_text`/`section_ref` (kept, still required) so Phase D can rank
      what to review.
  P8. COUNT RECONCILIATION. Per-patient immunizing-peptide records must match
      `n_peptides_administered` (the manual Keskin 13/10/20/7/11/15/17/10
      check, now enforced).
  P9. STRUCTURED VARIANT. `variant_type` (missense/frameshift/...) + HGVS-ish
      `protein_change` sit alongside the free-form `mutation`; frameshift vs
      missense neoantigens differ mechanistically and the indel length cap
      already in `_AA_PATTERN` depends on knowing which is which.
  P10. PMCID added and pattern-validated.
  P11 (v2.2). CONCOMITANT IMMUNOSUPPRESSION. A known confounder of vaccine
      immunogenicity (Keskin: 6/8 on dexamethasone during priming failed to
      respond; the 2 who weren't, responded). `ExtractedPatient.
      immunosuppression` records the FACT (drug class, agent, dose, timing) —
      never the paper's causal conclusion. With many papers the corpus reveals
      the correlation (e.g. corticosteroid-during-priming vs
      n_peptides_immunogenic == 0) instead of trusting any one paper's claim.
  P12 (v2.3). SPECIES + MHC NOMENCLATURE + VACCINE PLATFORM. The schema was
      human-clinical by construction; Li et al. 2021 (Genome Medicine) is mostly
      preclinical MOUSE work delivered by a DNA vaccine — neither representable
      before. `ExtractedPatient.species` (+ `model_system`) tags the subject;
      the MHC allele type is generalized (`_MHC_PATTERN`, `MhcAllele`) to admit
      murine H-2 (H-2Kb/Db/Kd/Ld) so a mouse cohort keeps its real restriction
      instead of dropping it; `vaccine_platform` (+ detail) records the modality.
  P13 (v2.4). PRECLINICAL ANTITUMOR EFFICACY. `clinical_outcome` is RANO and
      cannot express animal tumour-growth/survival readouts; Li et al. report
      vaccine efficacy that is conditional on the arm (E0771: nil alone, growth
      suppression with anti-PD-L1). `ExtractedPatient.preclinical_efficacy` is a
      list of arm-level results (readout/result/setting/combination), recorded
      as fact. Human outcomes stay in clinical_outcome/best_response. v2.4.1
      adds a minimal guard: a survival result needs a survival readout, and a
      survival readout can't carry a size/burden result.

Back-compat note: this is a breaking change to the peptide contract
(`peptides` -> `immunizing_peptides` + `epitopes`; `EvidenceTarget` widened;
evidence `peptide_paper_id` -> `immunizing_peptide_paper_id`/
`epitope_paper_id`). The extraction prompt (Layer-2) and Phase D ingest must
be updated together. All new structured fields are optional and default-empty
so a minimal extractor that only fills the legacy fields still validates.
"""
from __future__ import annotations

import re
from collections import Counter
from typing import Annotated, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    StringConstraints,
    field_validator,
    model_validator,
)

# ---------------------------------------------------------------------------
# Schema version — single source of truth. Bump on every contract change and
# stamp it onto each record's `source` at ingest, so a stored row always
# carries the schema it was validated against (drift forensics; Loop 3 in
# RUNBOOK.md). History: 2.0 split-peptide + integrity; 2.1 Measurement/raw
# affinities; 2.2 ConcomitantImmunosuppression (immunogenicity confounder);
# 2.3 species + general MHC nomenclature (murine H-2) + vaccine_platform;
# 2.4 preclinical antitumor efficacy (tumour-growth/survival by arm);
# 2.4.1 readout<->result consistency guard on PreclinicalEfficacy;
# 2.5.0 (P14) class-I epitope affinity rule relaxed to honor the lossless
#   Measurement lane — a class-I epitope must carry a predicted_affinity SLOT
#   but its value may be absent (unit='unknown' + raw) when the source reports
#   the predicted epitope+allele without a number (Rojas Supp Table 5: 232
#   MHC-I epitopes had no printed IC50/%rank). (P15) SurvivalOutcome — the
#   human analog of PreclinicalEfficacy for time-to-event endpoints
#   (RFS/OS/HR), since RANO clinical_outcome can't express adjuvant survival.
# 2.6.0 (P16) ResponseMagnitude — structured ELISpot/ICS/tetramer response size
#   (SFC/1e6, % of parent, stimulation index, or an ordinal -/+/++/+++/++++
#   grade), replacing free-text magnitudes in evidence.assay_detail. A SIBLING
#   of Measurement, not a reuse: magnitude units stay un-confusable from affinity
#   units and an ordinal grade has a home a float value cannot give it.
# 2.7.0 (P17) CuratorNote — an OPTIONAL, paper-level container for curatorial
#   commentary (challenge/decision/caveat/highlight). It is NOT an extracted fact
#   and is never model-generated at render time; the report prints it verbatim.
#   needs_review defaults TRUE (inverted from extracted entities) and any `refs`
#   MUST resolve to known paper_local_ids — so commentary stays tied to real
#   entities and flagged for human sign-off. Every existing record has it empty.
# 2.7.1 (P18) MhcAllele now admits a class-II HLA-DP/DQ HETERODIMER PAIR
#   ('HLA-DPA1*01:03/DPB1*02:01'). DP/DQ restriction is an alpha+beta heterodimer;
#   the single-allele pattern could not hold the pair, so 29 Rojas class-II epitope
#   alleles were lost (17 unrepresentable pairs + 12 valid single-DRB neighbours
#   dropped by atomic add_table). Additive + back-compatible: every previously valid
#   single allele still validates; only the new paired form is now also accepted.
# 2.8.0 (P19) NeoantigenCandidate — the prediction -> selection -> outcome FUNNEL.
#   Records candidates at the same grain as an immunizing peptide but at an earlier
#   stage (candidate_status: predicted/selected/administered), each carrying its
#   prioritization `ranking_scores`. Those scores are a SIBLING value object
#   (PrioritizationScore + score_kind), NOT Measurement: MeasurementUnit is
#   affinity-only (nM/%rank/unknown), so a non-affinity funnel score (expression
#   TPM, pipeline quality score, pVACtools rank, agretopicity/DAI, VAF, clonality,
#   foreignness) forced onto Measurement would degrade to value=None+raw (a lossless
#   STRING, not a queryable number) and a quality-score vs a rank would be
#   type-indistinguishable — the same reason P16 made ResponseMagnitude a sibling.
#   A `selected_peptide_id` bridges a selected candidate to its ImmunizingPeptide
#   (whose evidence gives the outcome); a new `candidate` evidence target records
#   outcomes for candidates that were screened but not selected. Count
#   reconciliation is UNTOUCHED — candidates are not immunizing peptides — so the
#   funnel denominator can't dilute the IMP invariant. Turns the corpus into the
#   prediction benchmark the field is starved for. Additive/back-compatible:
#   candidates=[] default; every existing record still validates.
# 2.9.0 (P21) VaccineDelivery — per-patient delivery covariates (dose/adjuvant/
#   formulation/schedule + surgery->dose latency + doses-received) so a failed design
#   is distinguishable from a good design delivered under a response-limiting regimen.
#   A SIBLING of ConcomitantImmunosuppression (the vaccine's own administration, vs a
#   confounder against it), NOT an extension. One optional nested object on
#   ExtractedPatient (same per-patient grain as vaccine_platform) -> ZERO new section/
#   target/tool wiring. adjuvant = separately-added immunostimulant ONLY (named, not by
#   mechanism); an LNP/lipoplex is a formulation (RNA -> adjuvant='none'). Dose carries a
#   comparable number + a dose_basis tag (compare within-basis) + the lossless verbatim.
#   Additive/back-compatible: vaccine_delivery=None default; every existing record validates.
# 2.9.1 (P21) CohortLatency — paper-level cohort delivery/treatment latency distribution (e.g.
#   surgery->vaccine), the home for an UNLABELLED per-patient latency (Rojas Fig 1D) that the
#   per-patient VaccineDelivery field cannot hold without fabricating the patient<->value mapping.
#   median/benchmark are comparable numbers (within time_unit); range/list go in lossless raw.
#   Additive/back-compatible: cohort_latencies=[] default.
# 2.10.0 (P22) Response->outcome linkage + trial setting — roadmap item #3 (Nat Med 2026 asks for
#   "better proxies to estimate vaccine efficacy" + early-stage/(neo)adjuvant focus). TWO additive
#   pieces: (a) `trial_setting` — a PER-PATIENT clinical setting enum (adjuvant/neoadjuvant/
#   metastatic/...) so a survival/response endpoint's cross-trial comparability is no longer
#   confounded (an adjuvant RFS != a metastatic ORR); distinct from the PRECLINICAL EfficacySetting.
#   (b) `ClinicalBenefitSignal` — the response->benefit BRIDGE (epitope spreading, ctDNA dynamics,
#   intratumoral infiltration, antigen loss) that had no home: NOT ExtractedEvidence (its one-target
#   validator needs a peptide; these aren't peptide-specific), NOT SurvivalOutcome (not time-to-event),
#   NOT clinical_outcome (RANO only). Per-patient list on ExtractedPatient + optional paper-level list
#   on ExtractedPaper (mirrors VaccineDelivery/CohortLatency). NO peptide pointer (co-location +
#   associated_with_response is the link). Response->SURVIVAL linkage is UNCHANGED (already in
#   SurvivalOutcome.stratifier). `antigen_loss` readout doubles as the home for the deferred P20
#   clonality antigen-loss event. TimepointPhase relocated above ExtractedPatient for reuse.
#   Additive/back-compatible: trial_setting=None + clinical_benefit_signals=[] defaults.
# 2.11.0 (P20) Clonality + gene-level antigen dynamics — ranking item #4 (the dominant vaccine-failure
#   mechanism even when the response succeeds; Swanton Nat 2025). TWO additive pieces: (A) a new
#   `NeoantigenMutation` entity at the GENOMIC grain (gene/coords/VAF/HLA + clonality + emerged/lost
#   status), the funnel's upstream, so a mutation reported WITHOUT a peptide sequence (34903219 Table S6
#   'new neoantigen mutations in recurrent tumor') has a home — VAF is a list of per-timepoint VafPoint
#   (fraction [0,1] + lossless raw, NOT Measurement whose units are affinity-only). (B) per-mutation
#   `clonality`/`cancer_cell_fraction`/`wgd_timing` on `_PeptideCore` for peptides reporting clonality
#   inline. Vocab is literature-locked (clonal/subclonal headline, pre/post-WGD secondary; TRACERx
#   Nat 2023). Default ALWAYS unknown/None, NEVER clonal. HLA-LOH deferred (too sparse in vaccine papers
#   -> curator-note). NeoantigenMutation is NOT count-reconciled (mutations != IMP). `antigen_loss` on
#   ClinicalBenefitSignal stays the NARRATIVE claim; NeoantigenMutation is the tabulated per-mutation
#   fact. Additive/back-compatible: neoantigen_mutations=[] + the new peptide fields default None.
# 2.11.1 (P20.1) Evidence-completeness anchors — paper-level n_immunogenic_reported /
#   n_tested_negative_reported (the counts the PAPER states). A diagnostic over 11 Keskin + 15 Rojas
#   re-runs found the evidence-count variance is driven by enumeration granularity (negatives swing
#   6<->12), NOT duplicates — so the fix is an anchor + a canonical negative-granularity rule, not a
#   dedup merge. finalize nudges ONCE when recorded evidence materially disagrees with these anchors
#   (overridable allow_evidence_count_mismatch). Additive: both default None.
# 2.11.2 (P20.2) Clonality 0/1 coercion — Piece B's `clonality` field accepts a raw `Clonal` 0/1
#   indicator column (1->clonal, 0->subclonal) via a before-validator on both _PeptideCore and
#   NeoantigenMutation. Root cause (live, 39910301): the add_table mapping DSL is a direct
#   column->field copy with NO value transform, so a 0/1 `Clonal` column could never satisfy the
#   Literal and every clonality row was silently dropped (CCF landed, clonality stayed empty across
#   3 re-runs). Coerces ONLY an exact 0/1 (int/float/'0'/'1') + bool; a genuine clonal/subclonal/
#   unknown still validates and any other value still errors. NEVER infers clonality from CCF.
#   Additive/back-compatible (None default unchanged).
# 2.11.3 Override audit trail — paper-level `finalize_overrides_used: list[str]` records which SOFT
#   finalize guards were overridden to let a record through. Empty = clean. The live recall-guard test
#   (33064988 / 39972124) showed a turn-pressured agent OVERRIDES expensive completeness nudges rather
#   than recovering, so the override must leave a mark: the scale lane routes any non-empty record to
#   needs_review instead of clean extracted/. Additive/back-compatible ([] default).
# 2.11.4 Companion-paper reference — paper-level `companion_paper_ref: str | None`. Root cause (live,
#   39972124 "RNA neoantigen vaccines prime long-lived CD8+ T cells in PDAC"): a SECONDARY-ANALYSIS
#   paper that reuses a prior trial's neoantigens reports only COHORT COUNTS ("25 of 108 vaccine
#   neoantigens"; "we previously reported detailed characteristics of vaccine neoantigen selection¹")
#   and lists NO per-sequence manifest — the manifest + primary immunogenicity live in the companion
#   paper (here Rojas, PMID 37165196, already extracted at 232 IMP/464 epitopes). The agent set
#   n_selected_reported=108 / n_immunogenic_reported=25, could only name the few neoepitopes THIS paper
#   uses in functional/TCR experiments, and the peptide-recall + immunogenic-recall anchors fired → a
#   correct, complete secondary-paper extraction was forced to needs_review. Setting companion_paper_ref
#   (the prior paper's citation/PMID/DOI) asserts the manifest/characterization is deferred there, which
#   RELAXES those two recall anchors (the negative-grain anchor is untouched). Honest DB metadata too:
#   it links a secondary analysis to its primary trial report. Additive/back-compatible (None default).
# 2.11.5 Evidence-anchor REDESIGN + override TIERING (QC behavior; NO new fields, back-compatible). The
#   2026-06-09 iris validation batch routed ALL 7 papers to needs_review — the clean lane was empty and
#   the override-audit signal was noise. Root-caused to two false-firing nudges, both fixed in
#   agent_core._evidence_anchor_gap: (1) IMMUNOGENIC recall is now POOLING-AWARE — an immunogenic
#   target_kind='pool' row is expanded by its pool's member count before comparing to
#   n_immunogenic_reported (a pool of 7 counts as 7, not 1), since the paper counts members; (2) the
#   NEGATIVE side fires ONLY on material OVER-enumeration (>1.5x), never UNDER (the canonical grain rule
#   records only NAMED negatives — fewer than a cohort total is expected). Both sides also now exempt
#   under companion_paper_ref (completing 2.11.4). Plus override TIERING (agent_core.SOFT_OVERRIDES /
#   overrides_are_soft_only + scale/extract_one.py): a record whose overrides are ALL soft
#   (unknown_funnel_size / regimen_divergence — variance/metadata) stays in the clean lane with the
#   override still recorded; only a HARD override (or any unknown one) routes to needs_review. Validated
#   offline on the 7 batch records -> 5 clean / 2 needs_review (the 2 are genuine recall gaps).
# 2.12.0 (P23) CD4/CD8 SUBSET + evidence MHC CLASS — two ORTHOGONAL optional axes on ExtractedEvidence
#   (t_cell_subset: cd4|cd8|bulk_or_unknown; mhc_class: class_i|class_ii|not_determined). Subset is a
#   per-MEASUREMENT axis (CD4+CD8 to one long peptide = TWO rows), NOT on ImmunizingPeptide. Both NEVER
#   guessed: a non-unknown value needs a verbatim cue (ExtractedEvidence._subset_class_have_verbatim_token).
#   Class-II epitope minting is gated on a quoted restriction anchor (MinimalEpitope._class_ii_needs_
#   restriction_anchor). MinimalEpitope.mhc_class stays mandatory I/II. None (un-migrated) stays distinct
#   from the processed 'bulk_or_unknown'/'not_determined'. Additive/back-compatible: both default None.
#   Plus MinimalEpitope.mhc_class_inferred (bool, default False): an EXPLICIT "class was a heuristic call"
#   marker that exempts the anchor (must co-set needs_review) — used by deterministic table adapters so a
#   length-convention class is declared honestly, NOT laundered as a synthetic cue in quoted_text. Agent/
#   LLM extractions never set it, so the no-guess anchor stays strict for LLM output.
# 2.13.0 (A1) COHORT_KIND — ExtractedPatient.cohort_kind (patient|tumor_model|model_antigen_validation|
#   healthy_donor|other; default None=legacy). Tags a METHODOLOGICAL arm (model-antigen pipeline
#   validation, e.g. Li 33879241 P1) so it is queryably EXCLUDABLE from the "vaccinated patient"
#   denominator (agent_core._evidence_breadth_gap) instead of silently dropped. Distinct from species.
#   Additive/back-compatible (default None). Pairs with the B1 evidence-grain rule (manifest per-target
#   'No response' rows are a screening COUNT, not per-row evidence — see prompt_render/RULES.md).
# 2.13.1 FAITHFULNESS HARDENING (Fable review): (1) class-I MinimalEpitope >14 aa is a HARD reject
#   (_class_i_is_a_minimal_binder) — a deterministic long-peptide mislabel, not a soft nudge. (2) the
#   subset/class cue search drops sibling-provenance bleed (cue must be on the row's own quoted_text/
#   assay_detail/hla_allele) and the cd4/cd8 branches gain murine H-2 terms (a faithful mouse CD8/CD4 row
#   is no longer wrongly rejected). (3) cohort_kind 'model_antigen_validation'/'healthy_donor' now needs a
#   quoted cue (_methodological_cohort_kind_needs_cue) — the recall-denominator exclusion must be earned.
#   (4) _class_ii_minting_gap counts only NON-inferred class-II (a guess no longer silences the guard).
# 2.14.0 (#4) SCREENING BUCKET — new entity ScreeningReadout + ExtractedPaper.screening_readouts. A
#   prediction/target MANIFEST's per-target bulk readout (e.g. Rojas 232-target ELISpot: 200 "No
#   response") is STORED here (every stated fact kept) instead of collapsed to a count (the B1 data-loss)
#   OR dumped into `evidence` (~9x immunogenicity inflation). Screening is a SEPARATE list/table, never
#   UNIONed into evidence; manifest_outcome (response|no_response|not_evaluable) is distinct from
#   EvidenceOutcome. Precedence: a NAMED evidence row at the same (patient,target,assay) supersedes a
#   screening row (agent_core._screening_evidence_overlap drops the dup). Additive (default []).
# 2.15.0 TRIAL-CONTEXT axes — three additive, all-optional fields: (1) per-patient
#   ExtractedPatient.concomitant_therapy (list[ConcomitantTherapy]; co-administered ICB/chemo/RT/
#   targeted, sibling of immunosuppression; drug_class is the comparable axis); (2) paper-level
#   ExtractedPaper.safety_summary (SafetySummary; CTCAE-grade headline facts, DISTINCT from the
#   immunogenicity ResponseGrade; full per-event AE tables deferred); (3) per-patient tmb_value/
#   tmb_raw/msi_status. All default None/[] -> every pre-2.15 record still validates.
# 2.16.0 DATA_RESOLUTION — two additive paper-level fields that let a genuinely coarse-grained study be
#   admitted faithfully instead of held as under-extracted: (1) ExtractedPaper.data_resolution
#   (DataResolution enum, finest->coarsest; the grain the paper reports at); (2)
#   ExtractedPaper.peptide_manifest_present (the available-grain ANCHOR: was a per-sequence manifest in
#   the source?). finalize derives the ACHIEVED grain from content (agent_core.derive_data_resolution)
#   and, when it is coarser than per_sequence AND no manifest was available, exempts the per-sequence
#   recall/breadth anchors (n_selected_reported is then a CITED count, like companion_paper_ref). Both
#   default None -> every pre-2.16 record validates and gates exactly as before (conservative).
# ---------------------------------------------------------------------------
SCHEMA_VERSION = "2.16.0"

# ---------------------------------------------------------------------------
# Patterns — kept in lockstep with antvac.ids and the DB invariants.
# ---------------------------------------------------------------------------
# Length 8–50 admits Class I 8–14mers (the canonical MHC-I window),
# Class II long peptides (15–25mers; vaccine-style "long peptides"
# typically 25-mers per Cafri Methods: "Sequences composed of 25 aa with
# the mutation flanked by 12 normal aa on each side"), AND the indel/
# frameshift case the same Methods explicitly justifies: "For In/Del
# mutations ... all aa beyond the mutation until the first stop codon."
# Empirical ceiling in our 4 papers (cafri Sup Tables 2-4): KMT2C
# p.Y4586fs at 52 aa (Pt4303), SPATA31D1 p.G1456fs at 46 aa (Pt4289);
# all SNV peptides ≤26 aa. Cap of 55 covers the empirical max with
# small buffer. Tighten back if a future dataset uses shorter canonical
# sequences only; widen explicitly if a future indel exceeds 55.
_AA_PATTERN = r"^[ACDEFGHIKLMNPQRSTVWY]{8,55}$"
# antvac.ids.HLA_PATTERN allows trailing characters; we tighten here because
# extracted strings should be the canonical 4-digit allele exactly. Allow
# alphanumeric in the locus to admit Class II alleles (DRB1, DQA1, DPA1, ...)
# alongside the simpler A/B/C of Class I.
_HLA_PATTERN = r"^HLA-[A-Z][A-Z0-9]*\*\d{2,3}:\d{2,3}$"
# PATCH (v2.3): general MHC pattern — a SUPERSET of _HLA_PATTERN that also admits
# non-human MHC, so preclinical (esp. murine) cohorts can carry their real
# restriction instead of being silently dropped. Branches:
#   - human HLA class I & II   HLA-A*02:01, HLA-DRB1*04:01
#   - human class-II DP/DQ HETERODIMER pair  HLA-DPA1*01:03/DPB1*02:01   (v2.7.1)
#   - murine H-2 class I       H-2Kb, H-2Db, H-2Kd, H-2Ld   (K/D/L/Q + haplotype)
#   - murine H-2 class II      H-2IAb / I-Ab style
# Human HLA already validated under _HLA_PATTERN remains valid here (back-compat).
#
# PATCH (v2.7.1): class-II HLA-DP/DQ restriction is biologically a HETERODIMER of an
# alpha (DPA1/DQA1) and beta (DPB1/DQB1) chain; papers report the PAIR
# (e.g. 'HLA-DPA1*01:03/DPB1*02:01'). The single-allele branches above could not
# represent a pair, so such alleles were silently dropped (Rojas P28/P29: 17 epitopes,
# and atomic add_table then took 12 valid single-DRB neighbours down with them). The
# paired branch below stores the pair faithfully — same motivation as the v2.3 murine
# patch. DR is excluded on purpose: its alpha (DRA) is effectively monomorphic, so DR
# restriction is already reported (and valid) as the single beta allele 'HLA-DRB1*..'.
_MHC_PATTERN = (
    r"^("
    r"HLA-[A-Z][A-Z0-9]*\*\d{2,3}:\d{2,3}"                          # human HLA (class I & II)
    r"|HLA-D[PQ]A1\*\d{2,3}:\d{2,3}/D[PQ]B1\*\d{2,3}:\d{2,3}"       # human class-II DP/DQ pair
    r"|H-2[KDLQ][a-z]\d?"                                            # murine H-2 class I
    r"|H-2\s?I[AE][a-z]\d?"                                          # murine H-2 class II (incl. I-Ag7)
    r"|I-[AE][a-z]\d?"                                               # murine class II alt notation (incl. I-Ag7)
    r")$"
)
_PMID_PATTERN = r"^\d+$"
# PATCH (v2 P10): PMCID is a distinct identifier from PMID and is what PMC
# full-text / OA tooling keys on; fetched by pubmed.fetch, never invented.
_PMCID_PATTERN = r"^PMC\d+$"
_DOI_PATTERN = r"^10\.\d{4,9}/[-._;()/:A-Z0-9a-z]+$"
_NCT_PATTERN = r"^NCT\d{8}$"
# PATCH (v2 P9): permissive HGVS protein-change check. Deliberately loose —
# it only rejects obviously-non-HGVS strings (must start 'p.') while leaving
# the cross-paper-comparable axis to `variant_type`. The free-form `mutation`
# field is retained unchanged for the 'GENE:p.X' display form papers vary on.
_PROTEIN_CHANGE_PATTERN = r"^p\.\(?[A-Za-z0-9_*=>]+\)?$"
# UniProt accession: standard SwissProt-style (one of the two regular forms).
# FIX (audit Fix 1): the original lacked a wrapping group, so the top-level
# `|` left `^` anchoring only the first branch and `$` only the second —
# "P04637junktail" validated. Wrapping in (...) anchors the whole alternation.
_UNIPROT_PATTERN = (
    r"^("
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]"
    r"|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
    r")$"
)
# Q3-PHI allowlist: a recognizable paper-local patient label. Optional
# p/pt/patient/subject/case prefix, optional separator, 1-5 digits, optional
# single trailing letter. Rejects names/initials/surnames by construction.
# RECONCILED 2026-06 (Li et al. 2021): a second branch also accepts short
# de-identified study/model codes — a 1-6 letter token + 1-5 digits + optional
# trailing letter (e.g. 'GTB16', 'mel21', 'A1'). Still rejects names (no digit),
# initials (the explicit check below), and long MRNs (digits capped at 5).
_LOCAL_ID_OK = re.compile(
    r"^(?:p|pt|patient|subject|case)?[\s._-]*\d{1,5}[a-z]?$"
    r"|^[a-z]{1,6}\d{1,5}[a-z]?$",
    re.I,
)


# Reusable Annotated string types — keep field declarations short.
NonEmptyStr = Annotated[str, StringConstraints(min_length=1, strip_whitespace=True)]
ShortText = Annotated[str, StringConstraints(min_length=1, max_length=2000, strip_whitespace=True)]
# PATCH (v2.3): renamed from HlaAllele; validates any MHC system
# (human HLA or murine H-2) via _MHC_PATTERN. Field NAMES stay `hla_allele(s)`
# for back-compat, but the type now accepts non-human MHC.
MhcAllele = Annotated[str, StringConstraints(pattern=_MHC_PATTERN)]


def _canon_mhc(s: str) -> str:
    """PATCH (v2.3): canonicalize one MHC allele string. Bare murine class-I
    token ('Db','Kb','Kd','Ld') -> 'H-2Db'; uppercase human HLA; leave already-
    canonical H-2 / I-A forms untouched. Shared by the peptide- and patient-
    level normalizers so behaviour can't drift between them."""
    s = s.strip()
    m = re.match(r"^([KDLQ])([a-z]\d?)$", s)
    if m:
        return f"H-2{m.group(1)}{m.group(2)}"
    # PATCH (v2.7.1): class-II DP/DQ heterodimer pair 'HLA-DPA1*..:../DPB1*..:..' —
    # uppercase and collapse any spaces around the '/' so 'HLA-DPA1*01:03 / DPB1*02:01'
    # canonicalizes to the slash-joined form the paired pattern branch accepts.
    if "/" in s and s.upper().startswith("HLA"):
        return "/".join(part.strip() for part in s.upper().split("/"))
    if s.upper().startswith("HLA") or "*" in s:
        return s.upper()
    return s

# PATCH (v2): cross-cutting controlled axes.
MhcClass = Literal["I", "II"]
VariantType = Literal[
    "missense", "frameshift", "inframe_indel", "stop_loss",
    "splice", "fusion", "other",
]
ProvenanceKind = Literal["table", "figure", "prose", "lookup"]
# PATCH (v2.1; UPDATED v2.6): units a Measurement can carry — strictly
# binding-affinity units (nM, NetMHCpan %rank). The original v2.1 note floated
# migrating response magnitudes (SFC/1e6) ONTO Measurement; v2.6 decided against
# that — magnitudes live on the sibling `ResponseMagnitude` so affinity units
# stay un-confusable from spot counts, and an ordinal +/++/+++ grade gets a home
# a float `value` cannot give it. "unknown" is the lossless escape hatch — a raw
# token captured but not parsed to a known unit; preserved, never coerced.
MeasurementUnit = Literal["nM", "pct_rank", "unknown"]
# A prediction is NOT reported data is NOT a validated wet-lab readout. Keeping
# tier on the value preserves the predicted/T0 vs reported/validated distinction
# at the granularity of a single number (mirrors the Q9 evidence rule).
MeasurementTier = Literal["predicted", "reported", "validated"]
# ScoreKind (v2.8 P19) — the cross-paper-comparable metric TYPE of a funnel
# prioritization score (see PrioritizationScore below). Lives here with the other
# vocab Literals so the lockstep block can assert it. A funnel score is a DIFFERENT
# quantity from binding affinity (MeasurementUnit is affinity-only), which is why
# PrioritizationScore is a sibling of Measurement rather than a reuse — same
# rationale as ResponseMagnitude (P16).
ScoreKind = Literal[
    "affinity", "expression_tpm", "quality_score", "rank",
    "agretopicity_dai", "vaf", "clonality", "foreignness", "other",
]


class _Frozen(BaseModel):
    """Common config shared by every extracted model."""

    model_config = ConfigDict(frozen=True, str_strip_whitespace=True, extra="forbid")


# ---------------------------------------------------------------------------
# Provenance (v2 P7)
# ---------------------------------------------------------------------------
class Provenance(_Frozen):
    """One structured source pointer for an extracted fact.

    The legacy `quoted_text` + `section_ref` on each row remain the required,
    human-readable anchor. `Provenance` adds a *typed, repeatable* pointer so a
    fact backed by several places (e.g. Keskin dexamethasone in BOTH Supp
    Table 4a and the Results prose) records all of them, and so Phase D can
    weight a table-sourced value differently from a figure-read one.
    """

    kind: ProvenanceKind
    locator: NonEmptyStr            # 'Supp Table 4a', 'Fig 1b', 'Results ¶3', 'PubMed'
    quoted_text: ShortText | None = None


# ---------------------------------------------------------------------------
# Measurement (v2.1) — a unit-explicit, lossless, provenanced numeric claim
# ---------------------------------------------------------------------------
class Measurement(_Frozen):
    """A single quantitative measurement or prediction — unit-explicit and lossless.

    Replaces bare ``*_nm: float`` fields. Three reasons:

    1. UNIT IS DATA. nM and NetMHCpan %rank are different scales that live in
       the SAME source column (Keskin Supp Table 5 mixes them — ARHGAP35, a
       validated hit, is reported as %rank). Carrying the unit on the value
       makes "0.80 (%rank)" un-confusable with 0.80 nM, which collapses the
       former two-field P3 split into one well-typed object.
    2. LOSSLESS. ``raw`` keeps the verbatim source token regardless of whether
       parsing succeeded, so a future/better parser can re-derive ``value`` and
       retroactively upgrade already-stored records. ``value`` may be None
       (parse failed) as long as ``raw`` survives — a measurement is never
       silently dropped or coerced.
    3. PROVENANCED CLAIM. method + version + tier + source let contradictory
       claims about the same quantity (NetMHCpan v2.4 vs v4.1 vs the paper's
       own printed value) coexist side by side instead of overwriting.
    """

    value: float | None = None
    unit: MeasurementUnit
    # Lossless lane: the exact source token, e.g. "467.91", "0.80 (%rank)",
    # "<50 nM". REQUIRED whenever `value` is absent; recommended always.
    raw: ShortText | None = None
    method: NonEmptyStr | None = None                 # e.g. 'NetMHCpan'
    method_version: Annotated[str, StringConstraints(max_length=32)] | None = None
    tier: MeasurementTier = "predicted"
    source: Provenance | None = None
    confidence: int | None = Field(default=None, ge=1, le=5)

    @model_validator(mode="after")
    def _lossless_and_unit_consistent(self) -> Measurement:
        # Lossless rule: never store an empty measurement.
        if self.value is None and self.raw is None:
            raise ValueError(
                "Measurement carries neither a parsed value nor a raw token"
            )
        # A parsed number must have a known, range-valid unit.
        if self.value is not None:
            if self.unit == "unknown":
                raise ValueError("a parsed value requires a known unit, not 'unknown'")
            if self.unit == "nM" and self.value < 0:
                raise ValueError(f"nM value must be >= 0, got {self.value}")
            if self.unit == "pct_rank" and not (0.0 <= self.value <= 100.0):
                raise ValueError(f"pct_rank must be in [0, 100], got {self.value}")
        return self

    @classmethod
    def from_raw(
        cls,
        raw: str,
        *,
        method: str | None = None,
        method_version: str | None = None,
        tier: MeasurementTier = "predicted",
        source: Provenance | None = None,
        confidence: int | None = None,
    ) -> Measurement:
        """Reference parser: build a Measurement from a verbatim source token,
        preserving `raw` whether or not the numeric/unit parse succeeds.

        Deliberately conservative. The canonical parser may live in the
        extraction worker; keeping a default here makes the raw->parsed mapping
        testable and gives downstream a baseline. Unparseable tokens become
        unit='unknown', value=None — quarantined, never coerced. Improving this
        parser and re-running it over stored `raw` is the retroactive-upgrade
        loop the lossless lane exists for.
        """
        s = (raw or "").strip()
        m = re.match(r"-?\d+\.?\d*", s)
        if m is None:
            value: float | None = None
            unit: MeasurementUnit = "unknown"
        elif "rank" in s.lower() or "%" in s:
            value, unit = float(m.group()), "pct_rank"
        else:
            value, unit = float(m.group()), "nM"
        return cls(
            value=value, unit=unit, raw=raw, method=method,
            method_version=method_version, tier=tier, source=source,
            confidence=confidence,
        )


# ---------------------------------------------------------------------------
# ResponseMagnitude (v2.6 P16) — sibling of Measurement for T-cell response size
# ---------------------------------------------------------------------------
# A SEPARATE object, not a reuse of Measurement — three reasons the concrete
# ELISpot data forced (Keskin/Li/Rojas):
#   1. ORDINAL LANE. Murine vaccine tables grade responses -/+/++/+++/++++; that
#      is not a float and cannot live in Measurement.value. ResponseMagnitude
#      carries a parsed numeric lane (value+unit) AND an ordinal `grade` in ONE
#      object, so a response's size has one home whether the paper printed a
#      count or a symbol.
#   2. UN-CONFUSABLE UNITS. Magnitude units (SFC/1e6, % of parent, stimulation
#      index) are a different quantity from binding-affinity units (nM/%rank);
#      keeping them off MeasurementUnit means an affinity field can never
#      structurally accept a spot count — the same principle Measurement exists
#      for. (Supersedes the v2.1 "migrate SFC onto Measurement" note.)
#   3. LOSSLESS. Rojas Fig 1g gives per-neoantigen SFC per PATIENT without gene
#      labels, so polytopic responders' rows are value=None + a captured set.
ResponseMagnitudeUnit = Literal[
    "sfc_per_1e6", "percent_of_parent", "stimulation_index", "unknown",
]
ResponseGrade = Literal["negative", "low", "moderate", "high", "very_high"]
# CuratorNoteKind (v2.7 P17) — kept here with the other vocab Literals so the
# lockstep block below can assert it; the CuratorNote model lives just before
# ExtractedPaper (it needs the value-object/Provenance types defined first).
CuratorNoteKind = Literal["challenge", "decision", "caveat", "highlight"]


class ResponseMagnitude(_Frozen):
    """How big a measured T-cell response was — the structured home for the
    ELISpot / ICS / tetramer magnitude that previously fell into
    `ExtractedEvidence.assay_detail` free text. Attaches to ONE evidence row
    (one assay x timepoint).

    Any lane may be filled: a parsed `value` + `unit` (e.g. 2459 sfc_per_1e6);
    an ordinal `grade` (murine -/+/++/+++/++++ -> negative/low/moderate/high/
    very_high); or — losslessly — only `raw`, when the paper gives a magnitude
    that can't be parsed or mapped (a per-patient SFC set without gene labels).
    At least one of value / grade / raw must be present (never an empty size).
    """

    value: float | None = None
    unit: ResponseMagnitudeUnit = "unknown"
    grade: ResponseGrade | None = None
    # Lossless lane: verbatim source token, e.g. "2459", "++", ">2,000 SFC/1e6".
    # REQUIRED when neither value nor grade is given.
    raw: ShortText | None = None
    # ELISpot comparability factors — recorded as fact (papers vary), not assumed.
    background_subtracted: bool | None = None
    denominator: Annotated[str, StringConstraints(max_length=48)] | None = None  # '1e6 PBMCs', '1e6 splenocytes'
    tier: MeasurementTier = "reported"   # a magnitude is measured, never predicted
    source: Provenance | None = None
    confidence: int | None = Field(default=None, ge=1, le=5)

    @model_validator(mode="after")
    def _lossless_and_unit_consistent(self) -> ResponseMagnitude:
        # Never store an empty magnitude.
        if self.value is None and self.grade is None and self.raw is None:
            raise ValueError(
                "ResponseMagnitude carries no value, grade, or raw token"
            )
        # A parsed number needs a known, range-valid unit.
        if self.value is not None:
            if self.unit == "unknown":
                raise ValueError("a parsed magnitude value requires a known unit, not 'unknown'")
            if self.value < 0:
                raise ValueError(f"magnitude value must be >= 0, got {self.value}")
            if self.unit == "percent_of_parent" and self.value > 100:
                raise ValueError(f"percent_of_parent must be in [0, 100], got {self.value}")
        return self


# ---------------------------------------------------------------------------
# PrioritizationScore (v2.8 P19) — a SIBLING of Measurement, same rationale as
# ResponseMagnitude (P16): a funnel ranking score is a DIFFERENT quantity from a
# binding affinity. MeasurementUnit is affinity-only (nM/%rank/unknown), so a
# non-affinity score on Measurement would degrade to value=None+raw (a lossless
# STRING, not a queryable number) and a quality-score vs a rank from one tool
# would be type-indistinguishable. Giving scores their own typed home keeps
# Measurement clean and the funnel benchmark queryable.
# ---------------------------------------------------------------------------
class PrioritizationScore(_Frozen):
    """One prioritization/ranking score for a NeoantigenCandidate — a SIBLING of
    Measurement (not a reuse), same rationale as ResponseMagnitude (P16): funnel
    scores are a DIFFERENT quantity from binding affinity. MeasurementUnit is
    affinity-only, so a non-affinity score on Measurement would be value=None+raw
    (a lossless STRING, not a queryable number), destroying queryability; and a
    quality-score vs a rank from one tool would be type-indistinguishable.
    score_kind is the cross-paper-comparable metric TYPE; value is the queryable
    number; raw is the lossless verbatim token; method+method_version are
    load-bearing (scores aren't comparable across tool versions). At least one of
    value/raw must be present."""

    score_kind: ScoreKind
    name: NonEmptyStr | None = None
    value: float | None = None
    raw: ShortText | None = None
    method: NonEmptyStr | None = None
    method_version: Annotated[str, StringConstraints(max_length=32)] | None = None
    tier: MeasurementTier = "predicted"
    source: Provenance | None = None
    confidence: int | None = Field(default=None, ge=1, le=5)
    needs_review: bool = False

    @model_validator(mode="after")
    def _lossless(self) -> "PrioritizationScore":
        if self.value is None and self.raw is None:
            raise ValueError("PrioritizationScore carries neither a parsed value nor a raw token")
        return self


class _Extracted(_Frozen):
    """Base for every extracted *row* model.

    Carries the legacy provenance pair (kept REQUIRED for back-compat) plus the
    v2 optional structured provenance and an extraction-confidence signal.
    `ExtractedPaper` is a container and does NOT inherit this.
    """

    quoted_text: ShortText
    section_ref: NonEmptyStr  # 'Table 2', 'Supp Table S3', 'Results §1', etc.
    # PATCH (v2 P7): optional, multi-source, typed provenance. Empty list keeps
    # a minimal extractor (legacy quoted_text/section_ref only) valid.
    provenance: list[Provenance] = Field(default_factory=list, max_length=12)
    # PATCH (v2 P7): self-assessed extraction confidence, mirroring the decision
    # ledger's 1–5 scale (1 = settled/verified, 5 = judgment under ambiguity).
    # Phase D uses it to prioritize review; `needs_review` is the explicit flag
    # an extractor sets when it had to make a contestable call.
    confidence: int | None = Field(default=None, ge=1, le=5)
    needs_review: bool = False


# ---------------------------------------------------------------------------
# Clonality / antigen-dynamics vocab (schema v2.11 P20). Primary axis = clonal/subclonal (TRACERx);
# wgd_timing = secondary, mutation-relative; status = a mutation's antigen dynamics (emerged/lost).
# Default ALWAYS unknown/None, NEVER clonal. These power both the per-peptide clonality fields (below)
# and the gene-level NeoantigenMutation entity.
Clonality = Literal["clonal", "subclonal", "unknown"]
WgdTiming = Literal["pre_wgd", "post_wgd", "unknown"]
AntigenStatus = Literal["present", "emerged", "lost", "retained", "unknown"]


def _coerce_clonality(v: object) -> object:
    """Coerce a supplementary `Clonal` 0/1 indicator into the Clonality vocab.

    Source tables encode clonality as an integer flag (the universal `Clonal` column:
    1 = clonal, 0 = subclonal) -- but the `add_table` mapping DSL is a direct column->field
    copy with no value transform, so a raw 0/1 would never satisfy the Literal and the row
    would be dropped (the bug that left clonality empty on 39910301). Coerce ONLY an exact
    0/1 (int, float, or the strings '0'/'1'); everything else passes through untouched so a
    genuine 'clonal'/'subclonal'/'unknown' validates normally and an off-vocab value still
    errors. NEVER invents clonality from a CCF or any other field -- only an explicit flag.
    """
    if isinstance(v, bool):  # guard: bool is an int subclass; treat True/False as a flag
        return "clonal" if v else "subclonal"
    if isinstance(v, (int, float)) and v in (0, 1):
        return "clonal" if v == 1 else "subclonal"
    if isinstance(v, str) and v.strip() in ("0", "1"):
        return "clonal" if v.strip() == "1" else "subclonal"
    return v


# Peptide core (shared by both peptide levels) — v2 P1/P9
# ---------------------------------------------------------------------------
class _PeptideCore(_Extracted):
    """Fields common to immunizing peptides and minimal epitopes.

    `paper_local_id` is a stable identifier within the paper. PREFER the
    paper's own ID when one exists (e.g. Keskin Supp Table 5 carries
    '14362-007-IMP02' / '14362-007-EPT12A') — it is more traceable and
    survives re-extraction than a synthesized label. The canonical id
    (sha1-based) is computed at ingest via `antvac.ids`; the LLM never
    produces hashes.
    """

    paper_local_id: NonEmptyStr
    sequence: Annotated[str, StringConstraints(pattern=_AA_PATTERN, to_upper=True)]
    gene_symbol: NonEmptyStr | None = None
    uniprot_id: Annotated[str, StringConstraints(pattern=_UNIPROT_PATTERN)] | None = None
    is_neoantigen: bool
    # e.g. 'KRAS:G12D' or 'TP53:R175H'. Free-form because papers vary.
    mutation: Annotated[str, StringConstraints(max_length=64)] | None = None
    # PATCH (v2 P9): structured variant axis. `variant_type` is the
    # cross-paper-comparable key; `protein_change` is HGVS-ish (p.R175H,
    # p.L782fs). Frameshift vs missense matters: the _AA_PATTERN length ceiling
    # was justified specifically by frameshift indel lengths (see comment).
    variant_type: VariantType | None = None
    protein_change: (
        Annotated[str, StringConstraints(pattern=_PROTEIN_CHANGE_PATTERN, max_length=64)] | None
    ) = None
    # None = no HLA restriction published for this peptide. Synthetic long
    # vaccine peptides (e.g. Cafri's 23-25mers) are administered without a
    # single MHC restriction; forcing an allele here invites fabrication.
    # Class-I minimal epitopes still carry their (predicted/known) allele.
    hla_allele: MhcAllele | None = None
    # PATCH (v2.11 P20, Piece B): per-mutation clonality, for a peptide whose paper reports clonality
    # INLINE (no separate mutation-dynamics table). Default None = not reported, NEVER inferred clonal
    # (assuming clonality is the error that sinks a vaccine). When the peptide instead links to a
    # NeoantigenMutation, that row is the source of truth and these stay None (inherit at query time).
    clonality: Clonality | None = None
    cancer_cell_fraction: float | None = Field(default=None, ge=0, le=1)  # CCF; ~1 -> clonal
    wgd_timing: WgdTiming | None = None

    # FIX (audit Fix 2): StringConstraints applies `pattern` BEFORE `to_upper`,
    # so a lowercase/mixed-case sequence ("qtqkhldly") false-rejected against the
    # upper-only _AA_PATTERN. Normalize case in a before-validator instead of
    # relying on constraint ordering. Same treatment for HLA locus casing.
    @field_validator("sequence", mode="before")
    @classmethod
    def _normalize_sequence(cls, v: object) -> object:
        return v.upper().strip() if isinstance(v, str) else v

    @field_validator("hla_allele", mode="before")
    @classmethod
    def _normalize_mhc(cls, v: object) -> object:
        # PATCH (v2.3): MHC-aware normalization via the shared helper (expands
        # bare murine tokens, uppercases human HLA, leaves canonical H-2 as-is).
        return _canon_mhc(v) if isinstance(v, str) else v

    @field_validator("clonality", mode="before")
    @classmethod
    def _coerce_clonality(cls, v: object) -> object:
        # v2.11 P20 Piece B: accept a raw `Clonal` 0/1 column mapped via add_table.
        return _coerce_clonality(v)


# ---------------------------------------------------------------------------
# Immunizing peptide — the long, injected peptide that pools contain (v2 P1)
# ---------------------------------------------------------------------------
class ImmunizingPeptide(_PeptideCore):
    """One long peptide actually synthesized/administered — the *pool member*.

    This is the entity `ExtractedPeptidePool.member_peptide_ids` points at and
    `MinimalEpitope.parent_peptide_ids` points back to. For personalized
    neoantigen vaccines the peptide belongs to one patient; `patient_paper_id`
    records that (None for shared/off-the-shelf antigens). When set, it powers
    the per-patient count reconciliation in `ExtractedPaper` (v2 P8).
    """

    # PATCH (v2 P1/P8): personalized vaccines make peptides per-patient. Optional
    # so shared-antigen designs (paper-level library) still validate with None.
    patient_paper_id: NonEmptyStr | None = None


# ---------------------------------------------------------------------------
# NeoantigenCandidate (v2.8 P19) — the prediction->selection->outcome funnel
# ---------------------------------------------------------------------------
# Same peptide-shaped object as an ImmunizingPeptide, but recorded at an earlier
# funnel stage. The full ranked candidate set is the DENOMINATOR a prediction
# benchmark needs; the selected/administered subset is what actually became a
# vaccine. Candidates are NOT immunizing peptides and do NOT participate in
# _peptide_counts_reconcile — so recording a 300-candidate funnel cannot dilute
# the per-patient IMP invariant.
class NeoantigenCandidate(_PeptideCore):
    """A predicted neoantigen candidate with its prioritization scores and funnel
    status. `selected`/`administered` candidates SHOULD bridge to their
    ImmunizingPeptide via `selected_peptide_id` (whose evidence carries the
    immunogenicity outcome); a candidate that was screened but not selected gets
    its outcome through a `candidate`-target ExtractedEvidence row."""

    patient_paper_id: NonEmptyStr | None = None
    candidate_status: CandidateStatus
    # prioritization scores, using the PrioritizationScore SIBLING (score_kind +
    # value/raw + method + method_version + tier), NOT Measurement: a funnel score
    # (expression TPM, pipeline quality score, pVACtools rank, agretopicity/DAI,
    # VAF, clonality, foreignness) is a different quantity from binding affinity, so
    # score_kind keeps scores typed/comparable and a non-affinity score stays a
    # queryable number instead of degrading to a lossless string on MeasurementUnit.
    ranking_scores: list[PrioritizationScore] = Field(default_factory=list, max_length=20)
    # bridge to the chosen ImmunizingPeptide; resolution enforced at paper level.
    selected_peptide_id: NonEmptyStr | None = None

    @model_validator(mode="after")
    def _selected_link_requires_selected_status(self) -> NeoantigenCandidate:
        # a predicted-but-not-selected candidate cannot point at an immunizing
        # peptide; only selected/administered candidates may carry the bridge.
        if self.selected_peptide_id is not None and self.candidate_status == "predicted":
            raise ValueError(
                "selected_peptide_id is only valid when candidate_status is "
                "'selected' or 'administered', not 'predicted'"
            )
        return self


# ---------------------------------------------------------------------------
# Minimal epitope — the predicted/known short epitope within a long peptide
# ---------------------------------------------------------------------------
class MinimalEpitope(_PeptideCore):
    """A minimal epitope WITHIN one or more immunizing peptides (v2 P1/P2/P3).

    `parent_peptide_ids` is MANY-TO-MANY: overlapping/tiling long peptides can
    share an epitope (confirmed in Keskin Supp Table 5 — e.g. EPT2B sits in
    IMP02/03/04). Empty list = the paper reported the epitope without resolving
    which long peptide it came from (preserved honestly, not invented).

    `mhc_class` is mandatory because it governs which other fields are even
    meaningful: a class-I epitope (8–14mer) carries a NetMHCpan nM affinity or
    %rank; a class-II / CD4 epitope does NOT carry a class-I affinity and is
    typically restricted by HLA-DR/DQ/DP (often un-typed in the paper).

    Q9 (unchanged): in-silico prediction metadata are PREDICTIONS, never
    evidence — they must NEVER seed an ExtractedEvidence row. If the only thing
    known about an epitope is its predicted affinity, it has zero evidence rows
    (predicted/T0 tier), not a 'presented'/'positive' row.
    """

    # PATCH (v2 P1): which long peptide(s) this epitope is contained in.
    parent_peptide_ids: list[NonEmptyStr] = Field(default_factory=list, max_length=12)
    # PATCH (v2 P2): the missing axis.
    mhc_class: MhcClass
    # PATCH (v2.12): EXPLICIT marker that mhc_class was assigned by a deterministic HEURISTIC (e.g. a
    # length convention in a paper-specific table adapter) rather than read from a source-stated
    # restriction. Default False. An agent/LLM extraction NEVER sets this — so the class-II minting
    # anchor (_class_ii_needs_restriction_anchor) stays STRICT for LLM output (no length-guessing). A
    # record that sets it True is declaring "this is an audited heuristic call" and MUST also set
    # needs_review=True; the anchor is then exempt because the inference is named, not laundered through
    # quoted_text. This is the honest alternative to keyword-stuffing a synthetic class cue into a quote.
    mhc_class_inferred: bool = False
    wild_type_sequence: Annotated[str, StringConstraints(to_upper=True)] | None = None
    # PATCH (v2.1): affinities are Measurement claims (value + explicit unit +
    # lossless `raw` + method/version/tier/source), NOT bare nM floats. The
    # nM-vs-%rank split (former P3) is now carried by Measurement.unit, so a
    # "0.80 (%rank)" token can never be read as 0.80 nM. The lossless `raw`
    # lane lives inside the Measurement, enabling retroactive re-parsing.
    predicted_affinity: Measurement | None = None
    wild_type_affinity: Measurement | None = None

    @field_validator("wild_type_sequence", mode="before")
    @classmethod
    def _normalize_wt(cls, v: object) -> object:
        return v.upper().strip() if isinstance(v, str) else v

    @model_validator(mode="after")
    def _class_affinity_consistency(self) -> MinimalEpitope:
        # PATCH (v2 P2/P3, revised v2.1; RELAXED v2.5 P14): a class-I epitope is
        # a binding-defined claim, so it must carry a predicted_affinity SLOT —
        # but, honoring the lossless Measurement lane, that slot may be value-less
        # (unit='unknown' + a raw token / method) when the source reports the
        # predicted epitope + restricting allele WITHOUT a printed IC50/%rank
        # (Rojas Supp Table 5). Forcing a parsed number here previously discarded
        # 232 real MHC-I epitope+allele claims while the (also-number-less) MHC-II
        # epitopes passed — an asymmetry that penalized the more central CD8
        # restriction. The rule now requires the slot, and (if a value IS parsed)
        # that it be in a class-I unit. A class-II / CD4 epitope must still NOT
        # carry an nM-unit affinity (NetMHCpan nM is class-I only).
        if self.mhc_class == "I":
            pa = self.predicted_affinity
            if pa is None:
                raise ValueError(
                    f"class-I epitope {self.paper_local_id!r} needs a "
                    f"predicted_affinity (parsed nM/%rank, or a lossless "
                    f"unit='unknown' Measurement when the source omits the number)"
                )
            if pa.value is not None and pa.unit not in ("nM", "pct_rank"):
                raise ValueError(
                    f"class-I epitope {self.paper_local_id!r} has a parsed "
                    f"predicted_affinity that is not in nM or %rank"
                )
        else:  # class II
            if self.predicted_affinity is not None and self.predicted_affinity.unit == "nM":
                raise ValueError(
                    f"class-II epitope {self.paper_local_id!r} carries an nM "
                    f"predicted_affinity; NetMHCpan-style nM affinities are class-I only"
                )
        return self

    @model_validator(mode="after")
    def _class_ii_needs_restriction_anchor(self) -> MinimalEpitope:
        # PATCH (v2.12): minting a CLASS-II epitope is gated on a QUOTED RESTRICTION ANCHOR — either a
        # named DR/DP/DQ (or murine I-A/I-E) allele on hla_allele, OR a class-II-restriction phrase in
        # quoted_text (e.g. 'class II-restricted', 'HLA-DR', 'CD4'). Without one, a 'II' label is a
        # guess, not an extraction. Class-I is unaffected (its affinity slot is the anchor). Verified
        # against all reference_records: 0 legitimate class-II epitopes break (every one names an allele).
        # EXEMPTION: a record that explicitly declares mhc_class_inferred=True (a deterministic length/
        # heuristic call from a table adapter, always needs_review) is an AUDITED inference, not a silent
        # guess — it is exempt here. This keeps the rule fully strict for agent/LLM output (which never
        # sets the flag) without forcing the adapter to launder a synthetic cue into quoted_text.
        if self.mhc_class == "II" and not self.mhc_class_inferred:
            hay = (self.quoted_text or "").lower() + " " + (self.hla_allele or "").lower()
            if not self.hla_allele and not re.search(
                r"class[ \-]?ii|hla-?d|\bdr\b|\bdp\b|\bdq\b|drb|dpa|dpb|dqa|dqb|cd4|h-2i|i-[ae]", hay
            ):
                raise ValueError(
                    f"class-II epitope {self.paper_local_id!r} has neither a named DR/DP/DQ (or "
                    f"I-A/I-E) allele nor a class-II restriction cue in quoted_text — a class-II "
                    f"label requires a quoted restriction anchor (do not guess class II)"
                )
        # the inferred-class escape hatch is only honest if it routes to a human: enforce review.
        if self.mhc_class_inferred and not self.needs_review:
            raise ValueError(
                f"epitope {self.paper_local_id!r} sets mhc_class_inferred=True (a heuristic class call) "
                f"but needs_review is not set — an inferred MHC class MUST be flagged for review"
            )
        return self

    @model_validator(mode="after")
    def _class_i_is_a_minimal_binder(self) -> MinimalEpitope:
        # PATCH (v2.13, Fable review #1): a class-I MinimalEpitope is a MINIMAL binder — the canonical
        # MHC-I window is 8-11mers, and the whole gold corpus class-I epitopes are 8-12 aa. A class-I
        # "epitope" longer than 14 aa is not an epitope at all, it is a LONG IMMUNIZING PEPTIDE mislabeled
        # (the lite-lane Li failure: 19-29mer peptides emitted as class-I epitopes). This is a DETERMINISTIC
        # mislabel, not a judgement call, so it is a HARD reject (move it to immunizing_peptides and mint
        # the predicted minimal epitope). Class-II cores are legitimately long and are NOT constrained here.
        if self.mhc_class == "I" and len(self.sequence or "") > 14:
            raise ValueError(
                f"class-I epitope {self.paper_local_id!r} has a {len(self.sequence)}-aa sequence — a "
                f"class-I MinimalEpitope is an 8-11mer minimal binder (>14 aa is a long immunizing "
                f"peptide mislabeled as an epitope; record it in immunizing_peptides instead)"
            )
        return self


ClinicalOutcome = Literal["CR", "PR", "SD", "PD", "not_reported"]

# PATCH (v2.3): subject species + vaccine platform. The schema was human-clinical
# by construction; Li et al. 2021 is mostly preclinical mouse work delivered by a
# DNA vaccine, neither of which the schema could represent. `species` tags the
# subject/cohort; `vaccine_platform` records the modality (DNA vs peptide vs ...).
Species = Literal["human", "mouse", "rat", "non_human_primate", "other"]
# PATCH (v2.13 A1): kind of cohort an ExtractedPatient represents. DISTINCT from species — a mouse
# subject may be a real tumor_model cohort OR a model_antigen_validation arm. Lets the "vaccinated
# patient" denominator exclude methodological arms (model_antigen_validation/healthy_donor) without
# dropping the data. None = un-tagged (legacy/back-compat).
CohortKind = Literal["patient", "tumor_model", "model_antigen_validation", "healthy_donor", "other"]
# PATCH (v2.16 DATA_RESOLUTION): the finest immunogenicity grain the PAPER reports at, finest->coarsest.
# Admits genuinely coarse papers faithfully (tagged + queryable) instead of holding them as under-
# extracted by the per-sequence recall anchor. None = un-tagged (legacy/back-compat -> treated as
# per_sequence-expected; conservative). Paired with peptide_manifest_present (the available-grain anchor).
DataResolution = Literal[
    "per_sequence", "per_mutation", "per_target_gene", "cohort_summary", "clinical_only",
]
VaccinePlatform = Literal[
    "dna", "rna", "synthetic_long_peptide", "short_peptide",
    "dendritic_cell", "viral_vector", "other", "unspecified",
]
# PATCH (v2.10 P22): clinical trial setting — the disease/treatment context, a PER-PATIENT covariate
# (like vaccine_platform) so a combined-cohort paper isn't collapsed. None = not reported (never
# inferred); 'other' = a stated-but-unlisted setting. EfficacySetting (prophylactic/therapeutic) is
# the PRECLINICAL axis and is distinct — this is the human-trial enrollment context.
ClinicalTrialSetting = Literal[
    "adjuvant", "neoadjuvant", "perioperative", "metastatic",
    "locally_advanced", "recurrent", "other",
]
# PATCH (v2.10 P22): TimepointPhase relocated here from the Evidence vocab block so the
# ClinicalBenefitSignal value-object (below, before ExtractedPatient) reuses it without a forward ref.
TimepointPhase = Literal[
    "pre_vaccine", "on_treatment", "post_vaccine", "memory", "unspecified",
]


# PATCH (v2.11 P20, Piece A): gene-level antigen dynamics. The candidate funnel starts at a predicted
# PEPTIDE; a tumour mutation is one grain higher (mutation -> predicted peptides -> candidate -> IMP).
# NeoantigenMutation models that mutation directly so a paper reporting mutations WITHOUT peptide
# sequences (34903219 Table S6) has a home, and so clonality/VAF live on their natural carrier.
class VafPoint(_Frozen):
    """One variant-allele-frequency reading at a timepoint (e.g. primary vs recurrent tumour). VAF is
    a fraction in [0,1]; `raw` keeps the verbatim token (and any range/qualifier). At least one of
    value / raw must be present (never an empty VAF point)."""

    timepoint_phase: TimepointPhase | None = None
    timepoint_label: ShortText | None = None          # 'primary tumor', 'recurrent tumor', 'week 24'
    value: float | None = Field(default=None, ge=0, le=1)
    raw: ShortText | None = None

    @model_validator(mode="after")
    def _vaf_nonempty(self) -> VafPoint:
        if self.value is None and self.raw is None:
            raise ValueError("VafPoint carries neither a value nor a raw token")
        return self


class NeoantigenMutation(_Extracted):
    """A tumour neoantigen mutation at the GENOMIC grain (gene / coords / VAF / HLA) — the upstream of
    the candidate funnel, for mutations a paper reports WITHOUT a peptide sequence (34903219 Table S6
    'New neoantigen mutations in recurrent tumor'). Captures antigen DYNAMICS: `status='emerged'`
    (primary VAF 0 -> recurrent >0) and `'lost'` (the Swanton relapse-loss failure mode — a correctly
    clonal neoantigen the tumour deleted). NOT peptide-shaped, so a sequence-less mutation has a home;
    `peptide_ref` bridges to an IMP/candidate when the paper gives the peptide too (clonality then
    inherits from here). Clonality default is None, NEVER clonal."""

    paper_local_id: NonEmptyStr
    patient_paper_id: NonEmptyStr
    gene_symbol: NonEmptyStr | None = None
    genomic_change: ShortText | None = None           # verbatim, e.g. 'chr1:40705080 A>T' (no fixed coord schema)
    variant_type: VariantType | None = None
    hla_restrictions: list[MhcAllele] = Field(default_factory=list, max_length=20)
    clonality: Clonality | None = None
    cancer_cell_fraction: float | None = Field(default=None, ge=0, le=1)
    wgd_timing: WgdTiming | None = None
    vaf: list[VafPoint] = Field(default_factory=list, max_length=12)
    status: AntigenStatus | None = None
    peptide_ref: NonEmptyStr | None = None            # IMP/candidate paper_local_id, when the paper links one
    source: Provenance | None = None

    @model_validator(mode="after")
    def _mutation_identifiable(self) -> NeoantigenMutation:
        # never an anonymous mutation — at least name the gene or the genomic change
        if not (self.gene_symbol or self.genomic_change):
            raise ValueError("NeoantigenMutation needs at least gene_symbol or genomic_change")
        return self

    @field_validator("hla_restrictions", mode="before")
    @classmethod
    def _canon_restrictions(cls, v: object) -> object:
        if isinstance(v, list):
            return [_canon_mhc(a) if isinstance(a, str) else a for a in v]
        return v

    @field_validator("clonality", mode="before")
    @classmethod
    def _coerce_clonality(cls, v: object) -> object:
        # v2.11 P20: accept a raw `Clonal` 0/1 column mapped via add_table (genomic-grain).
        return _coerce_clonality(v)


# PATCH (v2.2): concomitant immunosuppression vocab. Captures a known
# immunogenicity confounder as a FACT, never the paper's causal conclusion.
ImmunosuppressantClass = Literal[
    "corticosteroid", "chemotherapy", "calcineurin_inhibitor", "mtor_inhibitor",
    "antimetabolite", "anti_tnf", "other_immunosuppressant", "unspecified",
]
ImmunosuppressionTiming = Literal[
    "pre_vaccine", "during_priming", "during_boost", "throughout",
    "post_vaccine", "unspecified",
]


class ConcomitantImmunosuppression(_Frozen):
    """One immunosuppressive exposure during the vaccine window, recorded as a
    FACT — never the paper's interpretation of its effect.

    Keskin Supp Table 4a: 6/8 patients on dexamethasone 2-4 mg/day during
    vaccine priming  ->  agent_class='corticosteroid', agent='dexamethasone',
    dose='2-4 mg/day', timing='during_priming'. The schema does NOT store
    "this caused the response to fail" — that correlation is derived across the
    corpus (steroid-during-priming vs n_peptides_immunogenic == 0), not asserted.
    """

    agent_class: ImmunosuppressantClass               # cross-paper-comparable axis
    agent: NonEmptyStr | None = None                  # specific drug, e.g. 'dexamethasone'
    dose: Annotated[str, StringConstraints(max_length=64)] | None = None
    timing: ImmunosuppressionTiming = "unspecified"
    source: Provenance | None = None


# PATCH (v2.15 — trial-context axis #1): concomitant systemic anti-cancer therapy. A SIBLING of
# ConcomitantImmunosuppression (same fact-not-conclusion pattern), but it records the CO-ADMINISTERED
# anti-cancer treatment (ICB / chemo / RT / targeted), NOT an immunogenicity confounder. `drug_class`
# is the cross-paper-comparable axis (its single highest-value field). NOT CombinationClass: that
# Literal is a preclinical ARM descriptor (monotherapy / plus_*), the wrong grain for a per-drug fact.
ConcomitantDrugClass = Literal[
    "checkpoint_inhibitor", "chemotherapy", "radiotherapy", "targeted", "other",
]
ConcomitantTherapyTiming = Literal["concurrent", "sequential", "unknown"]


class ConcomitantTherapy(_Frozen):
    """One co-administered systemic anti-cancer therapy during the vaccine course, recorded as a
    FACT — never the paper's interpretation of its contribution to outcome.

    Rojas et al. 2023 (autogene cevumeran PDAC): vaccine given with atezolizumab and mFOLFIRINOX ->
    two records: drug_class='checkpoint_inhibitor', agent='atezolizumab', timing='concurrent'; and
    drug_class='chemotherapy', agent='mFOLFIRINOX', timing='sequential'. The schema does NOT store
    "the checkpoint inhibitor drove the response" — that is derived across the corpus, not asserted.
    """

    drug_class: ConcomitantDrugClass                  # cross-paper-comparable axis (highest value)
    agent: NonEmptyStr | None = None                  # specific drug/regimen, e.g. 'atezolizumab'
    timing: ConcomitantTherapyTiming = "unknown"      # concurrent vs sequential with the vaccine
    line: Annotated[str, StringConstraints(max_length=40)] | None = None  # e.g. '1L', 'second-line'
    source: Provenance | None = None


# PATCH (v2.9 P21): vaccine DELIVERY covariates — a SIBLING of ConcomitantImmunosuppression
# (same fact-not-conclusion pattern), NOT a subtype. ConcomitantImmunosuppression records a
# confounder working AGAINST the vaccine; this records the vaccine's own administration so a
# 'bad design' can be told apart from a 'good design delivered under a response-limiting regimen'
# (Keskin dexamethasone + Pt1-3 no boost; Rojas surgery->dose latency; EVX-01 dose<->magnitude).
# adjuvant = a SEPARATELY-ADDED immunostimulant ONLY, named as papers report it (mechanism is
# derivable later via a name->TLR lookup, not inferred here); an LNP/lipoplex is a FORMULATION,
# not an adjuvant -> formulation_detail, and an RNA vaccine takes adjuvant='none'.
AdjuvantClass = Literal[
    "poly_iclc", "montanide", "gm_csf", "cpg", "none", "other", "unspecified",
]
DoseScheduleBasis = Literal["per_peptide", "per_pool", "total", "unspecified"]


class VaccineDelivery(_Frozen):
    """How the vaccine was administered for ONE patient — sibling of
    ConcomitantImmunosuppression, recorded as FACT not conclusion. Every field optional;
    `unknown`/None means 'not reported', never inferred (absent adjuvant != 'none'). Trial-constant
    regimen fields repeat across a paper's patients (same grain as `vaccine_platform`); the two
    per-patient fields genuinely vary. Dose is captured as a comparable number + a `dose_basis`
    tag (compare only within-basis) + the lossless verbatim, so a per-peptide µg is never silently
    averaged against a total."""

    # --- regimen (usually trial/arm-constant) ---
    adjuvant: AdjuvantClass | None = None                 # separately-added immunostimulant; RNA -> 'none'
    adjuvant_detail: Annotated[str, StringConstraints(max_length=200)] | None = None
    formulation_detail: Annotated[str, StringConstraints(max_length=200)] | None = None  # 'RNA-lipoplex, IV'
    dose_amount_raw: Annotated[str, StringConstraints(max_length=120)] | None = None     # LOSSLESS, e.g. '300 µg/peptide'
    dose_per_peptide_ug: float | None = Field(default=None, ge=0)  # the comparable number, when cleanly stated
    dose_basis: DoseScheduleBasis = "unspecified"
    n_priming_doses: int | None = Field(default=None, ge=0)
    n_boost_doses: int | None = Field(default=None, ge=0)
    schedule_detail: Annotated[str, StringConstraints(max_length=300)] | None = None     # verbatim regimen
    # --- per-patient delivery (genuinely varies) ---
    weeks_surgery_to_first_dose: float | None = Field(default=None, ge=0)   # Rojas latency
    n_doses_received: int | None = Field(default=None, ge=0)                # Keskin Pt1-3 -> 0 boosts
    source: Provenance | None = None

    @model_validator(mode="after")
    def _delivery_consistency(self) -> VaccineDelivery:
        # received cannot exceed planned, when both are known (mirrors cohort arithmetic).
        if self.n_priming_doses is not None or self.n_boost_doses is not None:
            planned = (self.n_priming_doses or 0) + (self.n_boost_doses or 0)
            if self.n_doses_received is not None and self.n_doses_received > planned:
                raise ValueError(
                    f"n_doses_received ({self.n_doses_received}) exceeds planned priming+boost "
                    f"doses ({planned})"
                )
        if self.adjuvant == "none" and self.adjuvant_detail:
            raise ValueError("adjuvant='none' cannot carry adjuvant_detail")
        return self


# PATCH (v2.10 P22): the response->clinical-benefit BRIDGE — the "better proxies to estimate vaccine
# efficacy" the Nat Med 2026 roadmap asks for. These benefit events have NO existing home: not
# ExtractedEvidence (its exactly-one-target validator requires a peptide; epitope spreading is BY
# DEFINITION reactivity to antigens NOT in the vaccine, ctDNA/infiltration have no peptide key), not
# SurvivalOutcome (not a time-to-event statistic), not clinical_outcome (RANO tumour-burden only).
# Response->SURVIVAL linkage is ALREADY expressible via SurvivalOutcome.stratifier='vaccine response'
# — this captures the NON-survival bridge only. `antigen_loss` doubles as the home for the deferred
# P20 clonality doc's antigen-loss-at-relapse event.
BenefitReadout = Literal[
    "epitope_spreading", "ctdna_dynamics", "tumor_infiltration", "antigen_loss", "other",
]
BenefitDirection = Literal[
    "increased", "decreased", "cleared", "persisted",
    "lost", "detected", "not_detected", "unchanged",
]


class ClinicalBenefitSignal(_Frozen):
    """One response->clinical-benefit bridge observation, recorded as FACT — never the paper's causal
    story. Lives per-patient on ExtractedPatient (epitope spreading / ctDNA are per-subject) AND,
    optionally, paper-level on ExtractedPaper for a cohort-aggregate signal (e.g. a trial-wide
    infiltration summary) — mirroring the VaccineDelivery (per-patient) / CohortLatency (paper-level)
    split. NO peptide pointer by design: the driving immune response already lives in
    ExtractedEvidence; co-location + `associated_with_response` is the link, and duplicating it would
    re-introduce the join fan-out the unique-id guards prevent.

    Keskin GBM: intratumoral neoantigen-specific T cells at relapse ->
      readout='tumor_infiltration', direction='detected', timepoint_phase='post_vaccine',
      associated_with_response=True.
    Rojas PDAC: ctDNA cleared in responders -> readout='ctdna_dynamics', direction='cleared'.
    """

    readout: BenefitReadout                            # cross-paper-comparable axis (WHAT)
    direction: BenefitDirection                        # the result/direction
    timepoint_phase: TimepointPhase | None = None      # REUSED Evidence vocab
    timepoint_label: Annotated[str, StringConstraints(max_length=64)] | None = None  # verbatim ('at relapse')
    magnitude: ResponseMagnitude | None = None         # REUSED measurement shape, when quantified
    # Did the PAPER tie this signal to vaccine response? A reported fact, NOT our inference.
    associated_with_response: bool | None = None
    note: Annotated[str, StringConstraints(max_length=300)] | None = None
    source: Provenance | None = None

    @model_validator(mode="after")
    def _readout_direction_consistency(self) -> ClinicalBenefitSignal:
        # Minimal, UNAMBIGUOUS-only guard (cf. PreclinicalEfficacy / ResponseMagnitude). 'cleared' is
        # a ctDNA/disease concept; it cannot describe antigen_loss or epitope_spreading, where loss/
        # gain is the axis. Everything else is a legitimate fact and is NOT constrained.
        if self.direction == "cleared" and self.readout in ("antigen_loss", "epitope_spreading"):
            raise ValueError(
                f"benefit direction 'cleared' is incoherent with readout={self.readout!r} "
                f"(use 'lost'/'increased'/'detected' etc.)"
            )
        return self


# PATCH (v2.4): preclinical antitumor-efficacy vocab. The schema's
# `clinical_outcome` is RANO (CR/PR/SD/PD) — meaningful for human solid-tumour
# response but unable to express the tumour-growth / survival readouts of animal
# studies. These axes describe one EFFICACY EXPERIMENT ARM as a FACT (what was
# measured, the result, the setting, what it was combined with) — never the
# paper's mechanistic interpretation.
EfficacyReadout = Literal[
    "tumor_growth", "survival", "tumor_free_fraction", "metastasis", "other",
]
EfficacyResult = Literal[
    "no_effect", "partial_inhibition", "growth_inhibition",
    "tumor_regression", "complete_response", "prolonged_survival", "not_reported",
]
EfficacySetting = Literal["prophylactic", "therapeutic", "unspecified"]
CombinationClass = Literal[
    "monotherapy", "plus_checkpoint_inhibitor", "plus_chemotherapy",
    "plus_radiotherapy", "plus_other", "unspecified",
]


class PreclinicalEfficacy(_Frozen):
    """One antitumor-efficacy experiment arm, as reported — the preclinical
    analog of `clinical_outcome`. Recorded as a FACT (arm -> readout -> result),
    not the paper's causal story.

    A list because efficacy depends on the treatment ARM. Li et al. 2021, E0771:
      - vaccine alone           -> readout='tumor_growth', result='no_effect',
                                   combination='monotherapy'
      - vaccine + anti-PD-L1    -> result='growth_inhibition',
                                   combination='plus_checkpoint_inhibitor'
    4T1.2: vaccine alone -> 'partial_inhibition' (note: model is anti-PD-L1
    resistant; combination not yet tested). For HUMAN cohorts the primary
    outcome belongs in `clinical_outcome`/`best_response`; leave this empty.
    """

    readout: EfficacyReadout                          # WHAT was measured
    result: EfficacyResult                            # the direction/magnitude
    setting: EfficacySetting = "unspecified"          # prophylactic vs therapeutic
    combination: CombinationClass = "monotherapy"     # arm: alone vs +ICB/+chemo
    combination_detail: Annotated[str, StringConstraints(max_length=120)] | None = None
    comparator: Annotated[str, StringConstraints(max_length=120)] | None = None  # e.g. 'control vector DNA'
    n_animals: int | None = Field(default=None, ge=0)
    timepoint: Annotated[str, StringConstraints(max_length=64)] | None = None
    statistic: Annotated[str, StringConstraints(max_length=120)] | None = None    # e.g. 'P=0.0381, ANOVA'
    note: Annotated[str, StringConstraints(max_length=300)] | None = None
    source: Provenance | None = None

    @model_validator(mode="after")
    def _readout_result_consistency(self) -> PreclinicalEfficacy:
        # PATCH (v2.4.1): flag only the UNAMBIGUOUS readout/result mismatches.
        # 'no_effect' / 'not_reported' / 'complete_response' are readout-agnostic
        # and always allowed; this is deliberately minimal so it never fires on a
        # sensible record — it just stops a survival result on a caliper readout
        # (and vice-versa), which would otherwise be a silent data-entry error.
        size_results = {"partial_inhibition", "growth_inhibition", "tumor_regression"}
        if self.result == "prolonged_survival" and self.readout != "survival":
            raise ValueError(
                f"result 'prolonged_survival' requires readout='survival', "
                f"got readout={self.readout!r}"
            )
        if self.readout == "survival" and self.result in size_results:
            raise ValueError(
                f"survival readout cannot carry a size/burden result "
                f"{self.result!r} (expected a survival-type result, e.g. "
                f"'prolonged_survival', 'no_effect', 'complete_response')"
            )
        return self


# PATCH (v2.5 P15): clinical time-to-event outcome vocab + model. The HUMAN
# analog of PreclinicalEfficacy. `clinical_outcome` is RANO (CR/PR/SD/PD),
# which is undefined for an ADJUVANT trial with no measurable disease — Rojas
# et al. 2023 (autogene cevumeran PDAC) reports recurrence-free survival as its
# primary result, with no RANO response to record. These axes let a survival
# endpoint be stored as a FACT (endpoint -> arm -> median/HR), not a conclusion.
SurvivalEndpoint = Literal[
    "rfs", "os", "landmark_rfs", "dfs", "pfs", "efs", "ttr", "other",
]
SurvivalTimeUnit = Literal["months", "weeks", "days", "years"]


class SurvivalOutcome(_Frozen):
    """One reported time-to-event outcome — the human analog of
    PreclinicalEfficacy, for endpoints RANO cannot express (RFS/OS/...).

    Modeled as ONE record per (endpoint x arm). A two-arm comparison is two
    records: the experimental arm carries the comparison stats (hazard_ratio,
    p_value, comparator_label); the reference arm is a plain median record.
    Lives at the PAPER level — a between-group comparison spans patient groups,
    not a single patient. Per-patient IPD (if a paper ever reports it) would be
    a separate per-patient addition; this captures what trials actually print.

    Rojas RFS-by-response:
      - arm='vaccine responders', n=8, not_reached=True, comparator='non-responders',
        hazard_ratio=0.08, hr_ci='0.01-0.4', p_value=0.003, test='log-rank',
        stratifier='vaccine response', median_followup_value=18.0
      - arm='non-responders', n=8, median_value=13.4   (reference arm)
    """

    endpoint: SurvivalEndpoint                         # cross-paper-comparable axis
    arm_label: NonEmptyStr | None = None               # 'vaccine responders'; None = whole cohort
    n_patients: int | None = Field(default=None, ge=0)
    median_value: float | None = Field(default=None, ge=0)   # in `time_unit`
    time_unit: SurvivalTimeUnit = "months"
    not_reached: bool = False                          # median NOT reached (distinct from None=not reported)
    events: int | None = Field(default=None, ge=0)     # recurrences / deaths
    median_followup_value: float | None = Field(default=None, ge=0)  # in `time_unit`
    # Comparison block — populated on the EXPERIMENTAL arm of a stratified analysis.
    comparator_label: NonEmptyStr | None = None        # reference arm, e.g. 'non-responders'
    hazard_ratio: float | None = Field(default=None, ge=0)
    hr_ci: Annotated[str, StringConstraints(max_length=40)] | None = None     # '0.01-0.4'
    p_value: float | None = Field(default=None, ge=0, le=1)
    statistical_test: Annotated[str, StringConstraints(max_length=48)] | None = None  # 'log-rank'
    stratifier: Annotated[str, StringConstraints(max_length=80)] | None = None        # 'vaccine response'
    note: Annotated[str, StringConstraints(max_length=300)] | None = None
    source: Provenance | None = None

    @model_validator(mode="after")
    def _survival_consistency(self) -> SurvivalOutcome:
        # A median was either reached (a number) or not (not_reached) — never both.
        if self.not_reached and self.median_value is not None:
            raise ValueError(
                "SurvivalOutcome: not_reached=True cannot coexist with a parsed "
                "median_value"
            )
        # A hazard ratio is a between-arm statistic; it needs a named comparator.
        if self.hazard_ratio is not None and self.comparator_label is None:
            raise ValueError(
                "SurvivalOutcome: hazard_ratio requires a comparator_label "
                "(a HR compares this arm to a reference arm)"
            )
        return self


# PATCH (v2.9.1 P21): COHORT-level delivery/treatment latency — the paper-level home for a
# manufacturing/logistics timing reported across patients (Rojas Fig 1D: per-patient "time to mRNA
# vaccine" values, UNLABELLED, + a benchmark target). Distinct from SurvivalOutcome (a time-to-
# clinical-EVENT, not a process latency) and from VaccineDelivery.weeks_surgery_to_first_dose (use
# THAT only when the paper labels a latency per patient; record the unlabelled cohort distribution
# HERE). The benchmark-vs-achieved framing is the manufacturing-feasibility signal the field cares
# about. Additive/back-compatible: cohort_latencies=[] default.
LatencyMetric = Literal[
    "surgery_to_first_vaccine", "surgery_to_first_treatment",
    "diagnosis_to_first_vaccine", "other",
]


class CohortLatency(_Frozen):
    """One cohort-level delivery/treatment latency distribution (paper-level), recorded as FACT.
    `median_value`/`benchmark_value` are the comparable numbers (compare within `time_unit`);
    range/IQR and the raw value list go losslessly in `raw`. At least one of median_value / raw
    must be present (never an empty latency)."""

    metric: LatencyMetric                              # cross-paper-comparable axis
    n_patients: int | None = Field(default=None, ge=0)
    median_value: float | None = Field(default=None, ge=0)   # in `time_unit`
    time_unit: SurvivalTimeUnit = "weeks"
    benchmark_value: float | None = Field(default=None, ge=0) # the target (Rojas Fig 1D = 9)
    raw: ShortText | None = None                       # 'median 9.4, range 7.4-11' / the value list
    note: Annotated[str, StringConstraints(max_length=300)] | None = None
    source: Provenance | None = None

    @model_validator(mode="after")
    def _latency_nonempty(self) -> CohortLatency:
        if self.median_value is None and self.raw is None:
            raise ValueError("CohortLatency carries neither a median_value nor a raw token")
        return self


# PATCH (v2.15 — trial-context axis #3): per-patient tumour mutational burden + microsatellite status.
MsiStatus = Literal["mss", "msi_high", "unknown"]


class ExtractedPatient(_Extracted):
    """One trial subject as reported in the paper.

    `paper_local_id` is whatever the paper itself uses ('P1', 'Patient
    3', 'pt07'). The fully-qualified `patient_id` written to the DB is
    `pmid:<PMID>:<paper_local_id>` — assembled at ingest time, not
    by the LLM.

    Patient identifiers must be paper-local labels — never names,
    initials, or medical record numbers. The Phase D validators reject
    fields that look name-like.
    """

    paper_local_id: NonEmptyStr
    indication: NonEmptyStr  # e.g. 'metastatic melanoma', 'PDAC'
    # PATCH (v2.10 P22): clinical trial setting (adjuvant/neoadjuvant/metastatic/...). The cross-trial
    # covariate that gates endpoint comparability (an adjuvant RFS is not a metastatic ORR). Per-patient
    # grain so a combined-cohort paper isn't collapsed. None = not reported, NEVER inferred.
    trial_setting: ClinicalTrialSetting | None = None
    # PATCH (v2.3): species of the subject/cohort. Defaults to 'human' (the
    # schema's origin and all pre-v2.3 records); mouse/other MUST be set
    # explicitly. A preclinical cohort is modeled as a patient with species set
    # and `model_system` naming the strain/cell-line/transgenic.
    species: Species = "human"
    model_system: NonEmptyStr | None = None   # e.g. 'E0771 / C57BL/6', 'HHD II transgenic'
    # PATCH (v2.13 A1): kind of cohort. None = un-tagged (legacy). Set 'model_antigen_validation' for a
    # methodological pipeline-validation arm using MODEL antigens (e.g. Li 33879241 P1: HLA-A2-transgenic
    # viral/gp100 optimization) — KEEP and TAG it, never drop it; the "vaccinated patient" denominator
    # (agent_core._evidence_breadth_gap) then excludes it. 'tumor_model' = a tumour-bearing animal cohort
    # (E0771/4T1.2); 'patient' = a real disease cohort. Distinct from species (a mouse can be either).
    cohort_kind: CohortKind | None = None
    # PATCH (v2.3): vaccine modality. For single-platform trials this is the
    # same across patients; per-patient placement also handles multi-platform
    # papers (Li et al. compares DNA vs SLP) and per-subject clarity.
    vaccine_platform: VaccinePlatform = "unspecified"
    vaccine_platform_detail: Annotated[str, StringConstraints(max_length=200)] | None = None
    # PATCH (v2.9 P21): vaccine delivery covariates (dose/adjuvant/formulation/schedule/latency),
    # so a failed design can be separated from one delivered under a response-limiting regimen.
    # One optional nested object (same per-patient grain as vaccine_platform); None = not reported.
    vaccine_delivery: VaccineDelivery | None = None
    # Cap admits full Class I + Class II typing: up to 6 Class I (2x A/B/C)
    # plus DRB1/3/4/5, DQA1/DQB1, DPA1/DPB1 — Rojas Pt29 alone has 13.
    # v2.3: items may be human HLA or murine H-2 (MhcAllele).
    hla_alleles: list[MhcAllele] = Field(default_factory=list, max_length=20)
    # Synthesized = peptides made for this patient (paper's library size).
    # Required because every paper reports at least this number, even if it's
    # just the per-row count in the supplementary library table.
    n_peptides_synthesized: int = Field(ge=0)
    # Administered = peptides actually given as the vaccine. Often == synthesized,
    # but Keskin Pt 4-8 (and likely other future papers) show a real distinction
    # when only a subset of synthesized peptides was administered. None when
    # the paper doesn't unambiguously report (preserves the gap honestly).
    n_peptides_administered: int | None = Field(default=None, ge=0)
    n_peptides_immunogenic: int = Field(ge=0)
    # PATCH (v2 P6): was the patient evaluated for vaccine-peptide immunogenicity
    # at all? Keskin Pt 1-3 received no booster and were never assessed — their
    # n_peptides_immunogenic=0 means "not tested", NOT "tested, none worked"
    # (Pt 4-6). None = the paper doesn't say. False REQUIRES n_immunogenic == 0.
    immunogenicity_assessed: bool | None = None
    # `clinical_outcome` is the OVERALL/final RANO status (Keskin: PD for all 8).
    clinical_outcome: ClinicalOutcome | None = None
    # PATCH (v2): best RANO response at any timepoint, when the paper reports it
    # per patient (distinct from final status; many trials report only one).
    best_response: ClinicalOutcome | None = None
    # PATCH (v2.2): concomitant immunosuppression during the vaccine window — a
    # known immunogenicity confounder, stored as fact not conclusion. Empty list
    # = none reported. (Keskin Pt 1-6: corticosteroid during priming; Pt 7-8: none.)
    immunosuppression: list[ConcomitantImmunosuppression] = Field(
        default_factory=list, max_length=12
    )
    # PATCH (v2.4): preclinical antitumor efficacy — the preclinical analog of
    # clinical_outcome, which (being RANO) can't express tumour-growth / survival
    # readouts. A LIST because efficacy depends on the treatment arm (Li E0771:
    # vaccine alone -> no_effect; vaccine + anti-PD-L1 -> growth_inhibition).
    # Human cohorts use clinical_outcome/best_response and leave this empty.
    preclinical_efficacy: list[PreclinicalEfficacy] = Field(
        default_factory=list, max_length=12
    )
    # PATCH (v2.10 P22): the per-patient response->clinical-benefit bridge (epitope spreading, ctDNA,
    # intratumoral infiltration, antigen loss). The non-survival 'efficacy proxy' axis; response->
    # survival linkage stays in SurvivalOutcome.stratifier. Empty = none reported. Nested on the patient
    # -> rides the existing patients add-path (no new SECTION_MODEL/target/tool wiring).
    clinical_benefit_signals: list[ClinicalBenefitSignal] = Field(
        default_factory=list, max_length=20
    )
    # PATCH (v2.15 — trial-context axis #1): co-administered systemic anti-cancer therapy (ICB/chemo/RT/
    # targeted) during the vaccine course — sibling of `immunosuppression`, recorded as FACT not
    # conclusion. Near-universal in modern trials (Rojas = vaccine + atezolizumab + mFOLFIRINOX). Empty
    # list = none reported.
    concomitant_therapy: list[ConcomitantTherapy] = Field(
        default_factory=list, max_length=12
    )
    # PATCH (v2.15 — trial-context axis #3): tumour mutational burden + microsatellite status. The
    # genomic context that gates neoantigen-load interpretation. None = not reported, NEVER inferred.
    tmb_value: float | None = Field(default=None, ge=0)               # comparable number, when cleanly stated
    tmb_raw: Annotated[str, StringConstraints(max_length=64)] | None = None  # LOSSLESS, e.g. '12 mut/Mb'
    msi_status: MsiStatus | None = None

    @model_validator(mode="after")
    def _check_cohort_arithmetic(self) -> ExtractedPatient:
        # Synthesized is the upper bound: a peptide must have been made before
        # it can be administered or scored immunogenic.
        if self.n_peptides_immunogenic > self.n_peptides_synthesized:
            raise ValueError(
                f"n_peptides_immunogenic ({self.n_peptides_immunogenic}) > "
                f"n_peptides_synthesized ({self.n_peptides_synthesized})"
            )
        if (
            self.n_peptides_administered is not None
            and self.n_peptides_administered > self.n_peptides_synthesized
        ):
            raise ValueError(
                f"n_peptides_administered ({self.n_peptides_administered}) > "
                f"n_peptides_synthesized ({self.n_peptides_synthesized})"
            )
        # PATCH (v2 P6): "not assessed" cannot coexist with a positive count.
        if self.immunogenicity_assessed is False and self.n_peptides_immunogenic != 0:
            raise ValueError(
                f"immunogenicity_assessed=False but n_peptides_immunogenic="
                f"{self.n_peptides_immunogenic} (not-tested must be 0)"
            )
        return self

    @model_validator(mode="after")
    def _methodological_cohort_kind_needs_cue(self) -> ExtractedPatient:
        # PATCH (v2.13, Fable review #6): tagging a cohort 'model_antigen_validation' or 'healthy_donor'
        # REMOVES it from the "vaccinated patient" recall denominator (agent_core._evidence_breadth_gap),
        # so it can SILENCE an under-extraction. That exclusion must be earned: require a verbatim cue in
        # the cohort's own text (quoted_text / indication / model_system). 'patient' / 'tumor_model' /
        # 'other' / None assert nothing extra and are exempt.
        if self.cohort_kind in ("model_antigen_validation", "healthy_donor"):
            hay = " ".join(t for t in (self.quoted_text, self.indication, self.model_system) if t).lower()
            cue = (r"model[ \-]?antigen|model system|pipeline|optimi[sz]|validation|well[ \-]?characteri[sz]|"
                   r"transgenic|polyepitope|surrogate"
                   if self.cohort_kind == "model_antigen_validation"
                   else r"healthy (donor|volunteer|control|subject)|\bhd\d|normal donor")
            if not re.search(cue, hay):
                raise ValueError(
                    f"cohort {self.paper_local_id!r} is tagged cohort_kind={self.cohort_kind!r} but no "
                    f"cue for it appears in quoted_text/indication/model_system — a methodological tag "
                    f"that excludes the cohort from recall checks must be quotable (do not guess it)"
                )
        return self

    @field_validator("paper_local_id")
    @classmethod
    def _no_real_names(cls, v: str) -> str:
        # Q3-PHI (resolved: allowlist). A denylist of name shapes always leaks
        # (it missed 'JM Smith', 'Smith', lowercase names). Instead REQUIRE the
        # id to look like a paper-local label and reject anything else — a leaked
        # patient name is unrecoverable, a false quarantine is not.
        # Accepted: '4251', 'P1', 'Pt 8', 'Patient 3', 'pt07', 'subject12', etc.
        # NOTE: rejects alpha-prefixed labels like 'mel-21' / 'A1'. If your
        # papers use those, widen _LOCAL_ID_OK (documented tradeoff).
        if re.search(r"\b[A-Z]\.\s*[A-Z]\.", v):
            raise ValueError(f"paper_local_id looks like initials: {v!r}")
        if not _LOCAL_ID_OK.match(v.strip()):
            raise ValueError(
                f"paper_local_id {v!r} is not a recognizable paper-local label "
                f"(expected e.g. 'P1', 'Patient 3', 'pt07', '4251')"
            )
        return v

    @field_validator("hla_alleles", mode="before")
    @classmethod
    def _normalize_allele_list(cls, v: object) -> object:
        # PATCH (v2.3): canonicalize each allele (expand bare murine tokens etc.)
        # via the shared helper, so a mouse cohort's H-2 type validates.
        if isinstance(v, list):
            return [_canon_mhc(a) if isinstance(a, str) else a for a in v]
        return v


# ---------------------------------------------------------------------------
# Evidence
# ---------------------------------------------------------------------------
# Controlled vocabularies are imported from vocab.py — the SINGLE source of
# truth shared with the Layer-2 prompt builder, so prompt and schema cannot
# drift (resolves QV; see SCHEMA.md). `Literal[*tuple]` unpacks the vocab.
from typing import get_args as _get_args  # noqa: E402  (local import keeps vocab optional)
from . import vocab as _vocab  # noqa: E402

EvidenceAssay = Literal[
    "elispot", "tetramer", "ics", "flow_cytometry", "ms", "binding",
    "tcr_reporter", "vaccine_response", "clinical_response",
]
EvidenceOutcome = Literal[
    "immunogenic", "not_immunogenic", "pre_existing", "not_assessed",
    "positive", "negative", "presented",
]
# PATCH (v2 P4): evidence can target a long immunizing peptide (CD4 response to
# the processed long peptide), a minimal epitope (CD8/class-I), or a pool. This
# directly resolves the Keskin linkage gap where CD4 hits had no clean epitope
# key and had to attach at the long-peptide level.
EvidenceTarget = Literal["immunizing_peptide", "epitope", "pool", "candidate"]
# CandidateStatus (v2.8 P19) — kept here with the other vocab Literals so the
# lockstep block can assert it; the NeoantigenCandidate model lives above (it is
# referenced as a forward ref there and resolved at model build).
CandidateStatus = Literal["predicted", "selected", "administered"]
# PATCH (v2.12): two orthogonal OPTIONAL axes on ExtractedEvidence. t_cell_subset records the
# responding compartment; mhc_class records the restricting MHC class. DISTINCT from MinimalEpitope's
# mandatory mhc_class (Literal["I","II"]) — an evidence row's class is the assay-attributed restriction
# (often not_determined for a bulk ELISpot), a different claim from an epitope's intrinsic class. Never
# guessed: set only with a verbatim cue (enforced by _subset_class_have_verbatim_token below).
TCellSubset = Literal["cd4", "cd8", "bulk_or_unknown"]
EvidenceMhcClass = Literal["class_i", "class_ii", "not_determined"]
AssayStimulation = Literal["ex_vivo", "in_vitro_expanded", "unspecified"]
AssayAntigenFormat = Literal[
    "peptide_pulsed", "minigene", "autologous_tumor", "tetramer", "unspecified",
]
# TimepointPhase moved UP to the patient-context Literals (near VaccinePlatform) in v2.10 P22 so the
# ClinicalBenefitSignal value-object (defined before ExtractedPatient) can reuse it without a forward
# ref; its lockstep assert still lives in the block below. (definition: see ~VaccinePlatform.)

# ---------------------------------------------------------------------------
# VOCAB LOCKSTEP — every controlled-vocabulary Literal MUST equal its vocab.py
# tuple, or import fails loudly (by design; resolves QV). vocab.py is the single
# source shared with the Layer-2 prompt builder; editing one side without the
# other is caught HERE at import, not silently at extraction time. Adding a new
# Literal field? add its tuple to vocab.py and one line below.
# NOTE: this consolidated block replaces the original two scattered asserts and
# now also guards EvidenceTarget — whose drift, undetected precisely because it
# was unguarded, is what prompted centralizing every axis.
# ---------------------------------------------------------------------------
def _assert_vocab(literal: object, tup_name: str) -> None:
    have = set(_get_args(literal))
    want = set(getattr(_vocab, tup_name))
    if have != want:
        raise AssertionError(
            f"schema Literal drifted from vocab.{tup_name}: "
            f"schema-only={sorted(have - want)}  vocab-only={sorted(want - have)}"
        )


_assert_vocab(EvidenceOutcome,         "OUTCOMES")
_assert_vocab(EvidenceAssay,           "ASSAYS")
_assert_vocab(EvidenceTarget,          "EVIDENCE_TARGETS")
_assert_vocab(AssayStimulation,        "ASSAY_STIMULATION")
_assert_vocab(AssayAntigenFormat,      "ASSAY_ANTIGEN_FORMAT")
_assert_vocab(TimepointPhase,          "TIMEPOINT_PHASE")
_assert_vocab(ClinicalTrialSetting,    "CLINICAL_TRIAL_SETTINGS")
_assert_vocab(BenefitReadout,          "BENEFIT_READOUTS")
_assert_vocab(BenefitDirection,        "BENEFIT_DIRECTIONS")
_assert_vocab(Clonality,               "CLONALITY_CLASSES")
_assert_vocab(WgdTiming,               "WGD_TIMINGS")
_assert_vocab(AntigenStatus,           "ANTIGEN_STATUSES")
_assert_vocab(MhcClass,                "MHC_CLASSES")
_assert_vocab(VariantType,             "VARIANT_TYPES")
_assert_vocab(ProvenanceKind,          "PROVENANCE_KINDS")
_assert_vocab(MeasurementUnit,         "MEASUREMENT_UNITS")
_assert_vocab(MeasurementTier,         "MEASUREMENT_TIERS")
_assert_vocab(ImmunosuppressantClass,  "IMMUNOSUPPRESSANT_CLASSES")
_assert_vocab(ImmunosuppressionTiming, "IMMUNOSUPPRESSION_TIMING")
_assert_vocab(Species,                 "SPECIES")
_assert_vocab(CohortKind,              "COHORT_KINDS")
_assert_vocab(DataResolution,          "DATA_RESOLUTIONS")
_assert_vocab(VaccinePlatform,         "VACCINE_PLATFORMS")
_assert_vocab(EfficacyReadout,         "EFFICACY_READOUTS")
_assert_vocab(EfficacyResult,          "EFFICACY_RESULTS")
_assert_vocab(EfficacySetting,         "EFFICACY_SETTINGS")
_assert_vocab(CombinationClass,        "COMBINATION_CLASSES")
_assert_vocab(SurvivalEndpoint,        "SURVIVAL_ENDPOINTS")
_assert_vocab(SurvivalTimeUnit,        "SURVIVAL_TIME_UNITS")
_assert_vocab(ResponseMagnitudeUnit,   "RESPONSE_MAGNITUDE_UNITS")
_assert_vocab(ResponseGrade,           "RESPONSE_GRADES")
_assert_vocab(CuratorNoteKind,         "CURATOR_NOTE_KINDS")
_assert_vocab(CandidateStatus,         "CANDIDATE_STATUSES")
_assert_vocab(TCellSubset,             "T_CELL_SUBSETS")
_assert_vocab(EvidenceMhcClass,        "EVIDENCE_MHC_CLASSES")
_assert_vocab(ScoreKind,               "SCORE_KINDS")
_assert_vocab(AdjuvantClass,           "ADJUVANT_CLASSES")
_assert_vocab(DoseScheduleBasis,       "DOSE_SCHEDULE_BASES")
_assert_vocab(LatencyMetric,           "LATENCY_METRICS")
_assert_vocab(ConcomitantDrugClass,    "CONCOMITANT_DRUG_CLASSES")
_assert_vocab(ConcomitantTherapyTiming, "CONCOMITANT_THERAPY_TIMING")
_assert_vocab(MsiStatus,               "MSI_STATUSES")


class ExtractedEvidence(_Extracted):
    """One evidence row linking a patient × target to an assay outcome.

    The (assay, stimulation, antigen_format, timepoint) tuple identifies a
    distinct *measurement* — two timepoints of the same assay on the same
    peptide×patient are two evidence rows, not one. See the v2 design doc.
    """

    patient_paper_id: NonEmptyStr   # ExtractedPatient.paper_local_id
    # Q8 (pooled evidence) + v2 P4: an evidence row targets EXACTLY ONE of
    # {immunizing_peptide, epitope, pool}. target_kind discriminates; the
    # validator enforces exactly-one. Default keeps the common epitope case terse.
    target_kind: EvidenceTarget = "epitope"
    immunizing_peptide_paper_id: NonEmptyStr | None = None  # iff target_kind == "immunizing_peptide"
    epitope_paper_id: NonEmptyStr | None = None             # iff target_kind == "epitope"
    pool_paper_id: NonEmptyStr | None = None                # iff target_kind == "pool"
    candidate_paper_id: NonEmptyStr | None = None           # iff target_kind == "candidate" (v2.8 P19)
    # None = the assay didn't establish (or the paper didn't report) an HLA
    # restriction for this measurement. A patient×peptide reactivity (e.g. an
    # ELISpot on a long vaccine peptide) is valid evidence even without a
    # named allele; HLA is stored for display/audit only (not in the DB key).
    hla_allele: MhcAllele | None = None
    # QV: provenance for the allele ('figure-4C', 'results-text', 'supp-table-2').
    hla_source: Annotated[str, StringConstraints(max_length=64)] | None = None
    assay: EvidenceAssay
    outcome: EvidenceOutcome
    # QV: separates vaccine-INDUCED from pre-existing without overloading outcome.
    # pre_existing rows must have vaccine_induced=False (enforced below).
    vaccine_induced: bool | None = None
    # PATCH (v2 P5): was the response MUTATION-SPECIFIC? Keskin COX18 was
    # recognized equally for mutant and wild-type — reactive but NOT a
    # neoantigen-specific response. False flags exactly that; True = preferential
    # mutant reactivity; None = the paper didn't test/report the WT comparison.
    mutation_specific: bool | None = None
    # v2 dimensions — optional, default None (back-compatible).
    assay_stimulation: AssayStimulation | None = None
    assay_antigen_format: AssayAntigenFormat | None = None
    timepoint_phase: TimepointPhase | None = None
    # Verbatim-ish timepoint wording for display/audit (e.g. "16 weeks",
    # "year 3-4.5"); timepoint_phase is the cross-paper-comparable axis.
    timepoint_label: Annotated[str, StringConstraints(max_length=64)] | None = None
    # Free-text tail for exotic conditions that don't fit the enums above
    # (e.g. "minigene response blocked by anti-DR antibody").
    assay_detail: Annotated[str, StringConstraints(max_length=200)] | None = None
    # PATCH (v2.6 P16): structured response magnitude (SFC/1e6, % of parent,
    # stimulation index, or an ordinal grade) — replaces stuffing the size into
    # assay_detail free text. None = the paper reported the response without a
    # magnitude (or it isn't applicable, e.g. a binding/MS row).
    magnitude: ResponseMagnitude | None = None
    # PATCH (v2.12): responding T-cell compartment. CD4 vs CD8 is a per-MEASUREMENT axis, NOT a peptide
    # axis — "both CD4 and CD8 respond to one long peptide" is TWO evidence rows, one per subset, never a
    # field on ImmunizingPeptide. None = un-populated; 'bulk_or_unknown' = processed but the paper stated
    # no subset cue (e.g. a bulk PBMC ELISpot). Set cd4/cd8 ONLY with a verbatim cue (validator below).
    t_cell_subset: TCellSubset | None = None
    # PATCH (v2.12): MHC class that restricted THIS measurement. None = un-populated; 'not_determined' =
    # processed, no class cue. DISTINCT from MinimalEpitope.mhc_class (an epitope's intrinsic I/II). Set
    # class_i/class_ii ONLY with a verbatim cue (validator below). Provenance tier (reported_restriction
    # vs inferred_from_subset) is carried by the backfill/ETL, not stored on the row.
    mhc_class: EvidenceMhcClass | None = None

    @model_validator(mode="after")
    def _exactly_one_target(self) -> ExtractedEvidence:
        # PATCH (v2 P4): generalized to three target kinds.
        # PATCH (v2.8 P19): + candidate target.
        chosen = {
            "immunizing_peptide": self.immunizing_peptide_paper_id,
            "epitope": self.epitope_paper_id,
            "pool": self.pool_paper_id,
            "candidate": self.candidate_paper_id,
        }
        want = chosen[self.target_kind]
        others = [v for k, v in chosen.items() if k != self.target_kind]
        if not want or any(others):
            raise ValueError(
                f"{self.target_kind}-target evidence needs exactly its own "
                f"*_paper_id set and the other two unset"
            )
        return self

    @model_validator(mode="after")
    def _pre_existing_not_vaccine_induced(self) -> ExtractedEvidence:
        # protocol hard rule 2: pre_existing means NOT vaccine-induced.
        if self.outcome == "pre_existing" and self.vaccine_induced is True:
            raise ValueError("pre_existing outcome cannot have vaccine_induced=True")
        return self

    @model_validator(mode="after")
    def _magnitude_outcome_consistency(self) -> ExtractedEvidence:
        # PATCH (v2.6 P16): minimal, UNAMBIGUOUS-only guard (cf. the
        # PreclinicalEfficacy readout<->result guard). A response the row calls
        # positive cannot simultaneously be graded 'negative'. Everything else
        # (e.g. not_immunogenic with a small nonzero count below threshold) is a
        # legitimate fact and is NOT constrained.
        if (
            self.magnitude is not None
            and self.magnitude.grade == "negative"
            and self.outcome in ("immunogenic", "positive")
        ):
            raise ValueError(
                f"evidence outcome={self.outcome!r} but magnitude.grade='negative' "
                f"(a positive response cannot have a negative-graded magnitude)"
            )
        return self

    @model_validator(mode="after")
    def _subset_class_have_verbatim_token(self) -> ExtractedEvidence:
        # PATCH (v2.12; tightened v2.13): a NON-unknown subset/class is a CLAIM about THIS measurement and
        # must be backed by a verbatim token on the row's OWN fields (quoted_text / assay_detail / the named
        # hla_allele). 'bulk_or_unknown' / 'not_determined' / None assert nothing and are exempt.
        # v2.13 faithfulness fixes (Fable review): (1) DROP the sibling-provenance bleed — a cue in a
        # provenance quote may be about a DIFFERENT target on the same row, so it no longer counts; the cue
        # must be on quoted_text/assay_detail/hla_allele. (2) Add the murine H-2 terms to the cd8/cd4
        # branches so a faithful "H-2Kb-restricted, CD8" (or "I-Ab, CD4") row is NOT wrongly rejected.
        hay = " ".join(
            t for t in (self.quoted_text, self.assay_detail, self.hla_allele) if t
        ).lower()
        if self.t_cell_subset == "cd4" and not re.search(r"cd4|helper|\bth\b|class[ \-]?ii|hla-?d|h-2i|i-[ae]", hay):
            raise ValueError(
                f"t_cell_subset='cd4' but no CD4/helper/class-II cue in the row's verbatim text"
            )
        if self.t_cell_subset == "cd8" and not re.search(r"cd8|cytotoxic|\bctl\b|class[ \-]?i\b|hla-?[abc]\*|h-2[kdlq]", hay):
            raise ValueError(
                f"t_cell_subset='cd8' but no CD8/cytotoxic/class-I cue in the row's verbatim text"
            )
        if self.mhc_class == "class_ii" and not re.search(r"class[ \-]?ii|hla-?d|\bdr\b|\bdp\b|\bdq\b|cd4|h-2i|i-[ae]", hay):
            raise ValueError(
                f"mhc_class='class_ii' but no class-II restriction cue in the row's verbatim text"
            )
        if self.mhc_class == "class_i" and not re.search(r"class[ \-]?i\b|hla-?[abc]\*|cd8|h-2[kdlq]", hay):
            raise ValueError(
                f"mhc_class='class_i' but no class-I restriction cue in the row's verbatim text"
            )
        return self


# ---------------------------------------------------------------------------
# Peptide pool (Q8 — pooled-peptide papers, e.g. Keskin pools A-D)
# ---------------------------------------------------------------------------
class ExtractedPeptidePool(_Extracted):
    """A named immunizing pool within a paper (e.g. Keskin's pool 'A'..'D').

    Pools are per-patient here. A pool-target ExtractedEvidence row references
    `paper_local_id`. peptide <-> pool is MANY-TO-MANY (a peptide can sit in
    multiple pools — confirmed in Keskin Supp Table 5): `member_peptide_ids`
    holds every peptide in the pool; the DB join table is built from it.

    PATCH (v2 P1): `member_peptide_ids` now references `ImmunizingPeptide`
    paper_local_ids (the long, injected peptides) — NOT minimal epitopes. A pool
    physically contains long peptides; epitopes live inside those.
    """

    paper_local_id: NonEmptyStr            # pool label, e.g. 'A', 'PoolC'
    patient_paper_id: NonEmptyStr          # pools are per-patient
    member_peptide_ids: list[NonEmptyStr] = Field(min_length=1)


# ---------------------------------------------------------------------------
# ScreeningReadout (v2.14, #4) — the SCREENING BUCKET
# ---------------------------------------------------------------------------
ManifestOutcome = Literal["response", "no_response", "not_evaluable"]
_assert_vocab(ManifestOutcome, "MANIFEST_OUTCOMES")


class ScreeningReadout(_Extracted):
    """One per-target row of a SCREENING/prediction MANIFEST (e.g. Rojas 37165196's 232-target ELISpot
    table). DISTINCT from ExtractedEvidence on purpose: a manifest's bulk "No response" denominator is a
    SCREENING fact, not a per-measurement immunogenicity claim, and dumping it into `evidence` inflates the
    immunogenicity axis ~9x. Storing it HERE keeps every stated fact (root-cause fix for the B1 data-loss)
    while the immunogenicity axis stays at the named-result grain.

    Separation, not a flag: screening lives in its OWN list/table and is NEVER UNIONed into evidence.
    PRECEDENCE: a target with a NAMED ExtractedEvidence row at the same (patient, target, assay) is the
    real claim and must NOT also appear here (a finalize de-dup drops the screening row) — see
    agent_core._screening_evidence_overlap.
    """

    patient_paper_id: NonEmptyStr | None = None   # cohort-level manifest may have no per-patient grain
    target_kind: EvidenceTarget = "candidate"     # manifests usually enumerate candidates/peptides
    immunizing_peptide_paper_id: NonEmptyStr | None = None
    epitope_paper_id: NonEmptyStr | None = None
    pool_paper_id: NonEmptyStr | None = None
    candidate_paper_id: NonEmptyStr | None = None
    assay: EvidenceAssay
    manifest_outcome: ManifestOutcome             # response | no_response | not_evaluable (the column value)
    magnitude: ResponseMagnitude | None = None
    hla_allele: MhcAllele | None = None
    timepoint_phase: TimepointPhase | None = None
    timepoint_label: Annotated[str, StringConstraints(max_length=64)] | None = None

    @model_validator(mode="after")
    def _exactly_one_target(self) -> ScreeningReadout:
        chosen = {
            "immunizing_peptide": self.immunizing_peptide_paper_id,
            "epitope": self.epitope_paper_id,
            "pool": self.pool_paper_id,
            "candidate": self.candidate_paper_id,
        }
        want = chosen[self.target_kind]
        others = [v for k, v in chosen.items() if k != self.target_kind]
        if not want or any(others):
            raise ValueError(
                f"{self.target_kind}-target screening row needs exactly its own *_paper_id set "
                f"and the others unset"
            )
        return self


# ---------------------------------------------------------------------------
# Paper
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# CuratorNote (v2.7 P17) — curatorial commentary, NOT an extracted fact
# ---------------------------------------------------------------------------
# Design intent (anti-hallucination): this is a PASSIVE container. The renderer
# prints it verbatim and makes NO model call; nothing is auto-generated at render
# time. A note is authored deliberately (by a human, or a separate scoped pass
# that reads only one paper's JSON + its source), never free-typed from a long
# chat context. Three structural safeguards keep it honest:
#   - `needs_review` defaults TRUE (inverted from every extracted entity), so a
#     note renders with an "unverified" marker until a human clears it.
#   - `refs`, when given, MUST resolve to known paper_local_ids (enforced by
#     ExtractedPaper._curator_refs_resolve) — commentary is tied to real entities.
#   - it carries no `quoted_text`/`section_ref` extraction contract; it is clearly
#     a curatorial layer ON TOP of the facts, not a fact.
# (CuratorNoteKind Literal is defined above with the other vocab Literals.)


class CuratorNote(_Frozen):
    """One piece of curatorial commentary about the paper. Optional, paper-level,
    and rendered verbatim — see the anti-hallucination note above. Empty for every
    record unless a human/scoped pass deliberately adds one."""

    kind: CuratorNoteKind
    text: ShortText                                   # authored prose, printed verbatim
    # paper_local_ids this note concerns (patients / peptides / epitopes / pools).
    # MUST resolve (ExtractedPaper._curator_refs_resolve); empty = a paper-level note.
    refs: list[NonEmptyStr] = Field(default_factory=list, max_length=40)
    confidence: int | None = Field(default=None, ge=1, le=5)
    # INVERTED default: a note is unverified until a human signs off.
    needs_review: bool = True
    source: Provenance | None = None


# PATCH (v2.15 — trial-context axis #2): minimal paper-level safety / adverse-event summary. CTCAE
# severity GRADE here is DISTINCT from the response-magnitude `grade` (ResponseGrade) — that one is an
# immunogenicity tier, this one is clinical-toxicity severity (CTCAE 1-5). Full per-event AdverseEvent
# tables are DEFERRED (a future per-event entity); this captures the headline safety facts every trial
# reports so they aren't lost. Every field optional; None = not reported, NEVER inferred.
class SafetySummary(_Frozen):
    """Headline safety facts for ONE paper, recorded as FACT. DEFERRED: full per-event AE tables
    (grade-by-event, attribution-by-event) get a dedicated AdverseEvent entity later."""

    max_related_grade: int | None = Field(default=None, ge=1, le=5)   # highest treatment-related CTCAE grade
    any_grade3plus_related: bool | None = None                        # any related grade >= 3 AE reported
    n_patients_with_related_ae: int | None = Field(default=None, ge=0)
    irae_present: bool | None = None                                  # immune-related AE reported (combo trials)
    raw: Annotated[str, StringConstraints(max_length=600)] | None = None  # lossless verbatim safety sentence
    source: Provenance | None = None

    @model_validator(mode="after")
    def _safety_consistency(self) -> SafetySummary:
        # Minimal, UNAMBIGUOUS-only guard (cf. ClinicalBenefitSignal). If a max related grade is stated
        # AND any_grade3plus_related is explicitly False, grade must be < 3 — otherwise the two facts
        # contradict. Everything else (None on either side) is a legitimate not-reported state.
        if (
            self.max_related_grade is not None
            and self.any_grade3plus_related is False
            and self.max_related_grade >= 3
        ):
            raise ValueError(
                f"max_related_grade={self.max_related_grade} contradicts "
                f"any_grade3plus_related=False"
            )
        return self


class ExtractedPaper(_Frozen):
    """Top-level container: one paper, all extracted entities."""

    pmid: Annotated[str, StringConstraints(pattern=_PMID_PATTERN)]
    doi: Annotated[str, StringConstraints(pattern=_DOI_PATTERN)] | None = None
    # PATCH (v2 P10): PMC identifier (full-text/OA key), fetched not invented.
    pmcid: Annotated[str, StringConstraints(pattern=_PMCID_PATTERN)] | None = None
    nct_id: Annotated[str, StringConstraints(pattern=_NCT_PATTERN)] | None = None
    journal: NonEmptyStr
    year: int = Field(ge=2000, le=2030)
    title: NonEmptyStr
    # Q1: cohort_size = TREATED/vaccinated patients (the operational definition).
    cohort_size: int = Field(ge=1)
    # Q1: optional enrolled count when it differs from treated (e.g. Keskin
    # 10 enrolled / 8 vaccinated); lets the tolerance check compare against the
    # right denominator without losing the enrolled number.
    n_enrolled: int | None = Field(default=None, ge=1)
    indication_summary: NonEmptyStr   # e.g. 'metastatic melanoma + PDAC'
    patients: list[ExtractedPatient] = Field(default_factory=list)
    # PATCH (v2 P1): the single `peptides` list is split into the two real
    # levels. `immunizing_peptides` are pool members; `epitopes` are the
    # predicted/known minimal epitopes within them.
    immunizing_peptides: list[ImmunizingPeptide] = Field(default_factory=list)
    epitopes: list[MinimalEpitope] = Field(default_factory=list)
    pools: list[ExtractedPeptidePool] = Field(default_factory=list)
    evidence: list[ExtractedEvidence] = Field(default_factory=list)
    # PATCH (v2.14 #4): the SCREENING BUCKET — per-target rows of a prediction/target MANIFEST
    # (every screened target's bulk readout, e.g. Rojas's 200 "No response"). KEPT SEPARATE from
    # `evidence` so the stated facts aren't lost (the B1 fix) without inflating the immunogenicity axis.
    # NEVER UNION this into evidence counts. Empty for papers with no enumerated screening manifest.
    screening_readouts: list[ScreeningReadout] = Field(default_factory=list, max_length=5000)
    # PATCH (v2.5 P15): paper-level time-to-event outcomes (RFS/OS/...). Empty
    # for measurable-disease trials that report RANO in clinical_outcome; used
    # by adjuvant / time-to-event trials whose primary result RANO can't hold.
    survival_outcomes: list[SurvivalOutcome] = Field(default_factory=list, max_length=40)
    # PATCH (v2.9.1 P21): cohort-level delivery/treatment latency distributions (e.g. surgery->
    # vaccine), the home for an UNLABELLED per-patient latency the per-patient field can't hold.
    cohort_latencies: list[CohortLatency] = Field(default_factory=list, max_length=12)
    # PATCH (v2.10 P22): paper-level home for a COHORT-AGGREGATE benefit signal (e.g. a trial-wide
    # infiltration/ctDNA summary not attributable to one patient) — sibling of the per-patient
    # ExtractedPatient.clinical_benefit_signals, mirroring CohortLatency vs per-patient VaccineDelivery.
    clinical_benefit_signals: list[ClinicalBenefitSignal] = Field(default_factory=list, max_length=20)
    # PATCH (v2.15 — trial-context axis #2): minimal paper-level safety/AE summary (CTCAE-grade facts,
    # distinct from the immunogenicity ResponseGrade). Full per-event AE tables deferred. None = not
    # reported.
    safety_summary: SafetySummary | None = None
    # PATCH (v2.7 P17): optional curatorial commentary, rendered verbatim. Empty
    # for every record unless deliberately authored. NOT an extracted fact.
    curator_notes: list[CuratorNote] = Field(default_factory=list, max_length=60)
    # PATCH (v2.8 P19): the candidate funnel (predicted -> selected -> administered).
    # Large by nature (hundreds of predicted candidates); not count-reconciled.
    candidates: list[NeoantigenCandidate] = Field(default_factory=list, max_length=5000)
    # PATCH (v2.11 P20): gene-level tumour neoantigen mutations (the funnel's genomic upstream), for
    # mutations reported WITHOUT a peptide (34903219 Table S6). Carries clonality/VAF/HLA + emerged/lost
    # status. Additive/back-compatible: neoantigen_mutations=[] default. Not count-reconciled (mutations
    # are not immunizing peptides), so the IMP invariant is untouched.
    neoantigen_mutations: list[NeoantigenMutation] = Field(default_factory=list, max_length=5000)
    # PATCH (v2.8 P19): funnel-completeness honesty signals — the candidate counts
    # the PAPER STATES (denominators), so a truncated funnel (50 of 322 predicted)
    # is distinguishable from a complete one and any "selection rate" /
    # "immunogenic-per-predicted" computed downstream isn't a silently wrong
    # denominator. No hard validator — a paper may publish only a sample.
    n_predicted_reported: int | None = Field(default=None, ge=0)
    n_selected_reported: int | None = Field(default=None, ge=0)
    # PATCH (v2.11.1 P20.1): evidence-completeness anchors — the COUNTS THE PAPER STATES for immunogenic
    # responses and for tested-negative (non-)responses. Diagnostic (Keskin 6<->12 swing) found the
    # run-to-run evidence variance is driven by how finely a run enumerates findings (esp. negatives),
    # not by duplicates. These honesty anchors let finalize compare the recorded evidence against the
    # paper's own count and nudge ONCE on a material gap (the funnel n_predicted_reported analogue). Not
    # hard-validated — pooled responses legitimately collapse the immunogenic count, hence a soft nudge.
    n_immunogenic_reported: int | None = Field(default=None, ge=0)
    n_tested_negative_reported: int | None = Field(default=None, ge=0)
    # PATCH (v2.11.4): companion/primary-paper reference for a SECONDARY-ANALYSIS paper. Set this (to the
    # prior paper's citation / PMID / DOI) ONLY when the paper EXPLICITLY defers its neoantigen selection
    # and primary immunogenicity characterization to a previously-published companion paper and reports
    # only COHORT COUNTS here (no per-sequence manifest) — e.g. 39972124 ("25 of 108 vaccine neoantigens";
    # "we previously reported detailed characteristics of vaccine neoantigen selection¹") whose manifest
    # lives in Rojas (PMID 37165196). When set, the peptide-recall and immunogenic-recall anchors are
    # RELAXED (n_selected_reported / n_immunogenic_reported are CITED cohort counts, not targets this
    # paper enumerates), so a complete secondary-paper extraction finalizes clean instead of false
    # needs_review. NOT an escape hatch for "I couldn't find the table" — only for a genuine deferral the
    # paper itself states. None = standalone paper. The negative-grain evidence anchor is unaffected.
    companion_paper_ref: str | None = Field(default=None, max_length=300)
    # PATCH (v2.16): paper-level DATA RESOLUTION — the finest immunogenicity grain THIS paper reports
    # at (see DataResolution; finest->coarsest). A genuinely coarse paper (per_mutation/cohort_summary
    # /clinical_only, e.g. 39762422 autogene cevumeran's 40 gene+mutation neoantigens + survival, no
    # peptide sequences) is admitted faithfully + tagged instead of held by the per-sequence recall
    # anchor. The agent sets it to the finest grain ACTUALLY reported; finalize cross-checks against the
    # achieved grain derived from content (agent_core.derive_data_resolution). None = un-tagged (legacy
    # -> treated as per_sequence-expected; conservative).
    data_resolution: DataResolution | None = None
    # The available-grain ANCHOR: True iff a peptide/epitope SEQUENCE manifest/table was present in the
    # source (paper or supplement). Distinguishes faithfully-coarse (no finer grain offered -> admit)
    # from under-extracted (a manifest existed but wasn't fully parsed -> still flag). None = unknown
    # (legacy). NOT an escape hatch: set it truthfully from source inspection.
    peptide_manifest_present: bool | None = None
    # PATCH (v2.11.3): audit trail of which SOFT finalize guards were OVERRIDDEN to let this record
    # through (e.g. 'allow_sparse_evidence', 'allow_peptide_count_mismatch'). Empty = no override used =
    # all completeness/consistency nudges passed clean. Non-empty means the record finalized DESPITE a
    # fired guard, so it is NOT clean-silver: the scale lane routes it to needs_review and a human must
    # confirm the override was legitimate (the live test showed a turn-pressured agent overrides rather
    # than do expensive recall, so this is the QC signal that keeps such records from landing silently).
    finalize_overrides_used: list[str] = Field(default_factory=list, max_length=20)

    @model_validator(mode="after")
    def _unique_local_ids(self) -> ExtractedPaper:
        # PATCH (v2): paper_local_id must be unique within each entity type, or
        # cross-references become ambiguous and the DB join silently fans out.
        for label, ids in (
            ("patients", [p.paper_local_id for p in self.patients]),
            ("immunizing_peptides", [i.paper_local_id for i in self.immunizing_peptides]),
            ("epitopes", [e.paper_local_id for e in self.epitopes]),
            ("pools", [p.paper_local_id for p in self.pools]),
            ("candidates", [c.paper_local_id for c in self.candidates]),
            ("neoantigen_mutations", [m.paper_local_id for m in self.neoantigen_mutations]),
        ):
            dupes = [k for k, n in Counter(ids).items() if n > 1]
            if dupes:
                raise ValueError(f"duplicate paper_local_id in {label}: {dupes}")
        return self

    @model_validator(mode="after")
    def _curator_refs_resolve(self) -> ExtractedPaper:
        # PATCH (v2.7 P17): a curator note's refs (when given) must point at real
        # entities in THIS paper — commentary can't reference things that aren't
        # extracted. Empty refs = a legitimate paper-level note.
        known = (
            {p.paper_local_id for p in self.patients}
            | {i.paper_local_id for i in self.immunizing_peptides}
            | {e.paper_local_id for e in self.epitopes}
            | {p.paper_local_id for p in self.pools}
            | {c.paper_local_id for c in self.candidates}
        )
        for note in self.curator_notes:
            bad = [r for r in note.refs if r not in known]
            if bad:
                raise ValueError(
                    f"curator_note ({note.kind}) references unknown paper_local_id(s): {bad}"
                )
        return self

    @model_validator(mode="after")
    def _cross_reference_check(self) -> ExtractedPaper:
        """Every evidence row resolves to a known patient + target; every pool's
        members and every epitope's parents resolve to known immunizing peptides."""
        patient_ids = {p.paper_local_id for p in self.patients}
        imp_ids = {i.paper_local_id for i in self.immunizing_peptides}
        epitope_ids = {e.paper_local_id for e in self.epitopes}
        pool_ids = {p.paper_local_id for p in self.pools}
        candidate_ids = {c.paper_local_id for c in self.candidates}

        # PATCH (v2.8 P19): candidate ids must be disjoint from every other entity
        # id — a candidate bridges to an IMP by id and curator refs resolve against
        # the union, so a collision would make the join ambiguous.
        overlap = candidate_ids & (patient_ids | imp_ids | epitope_ids | pool_ids)
        if overlap:
            raise ValueError(
                f"candidate paper_local_id(s) collide with other entities: {sorted(overlap)}"
            )

        # PATCH (v2.8 P19): a selected/administered candidate's bridge must resolve
        # to a real immunizing peptide (the funnel join can't dangle).
        for c in self.candidates:
            if c.selected_peptide_id is not None and c.selected_peptide_id not in imp_ids:
                raise ValueError(
                    f"candidate {c.paper_local_id!r} selected_peptide_id "
                    f"{c.selected_peptide_id!r} does not resolve to a known immunizing peptide"
                )

        for ev in self.evidence:
            if ev.patient_paper_id not in patient_ids:
                raise ValueError(
                    f"evidence references unknown patient_paper_id "
                    f"{ev.patient_paper_id!r}"
                )
            if ev.target_kind == "immunizing_peptide":
                if ev.immunizing_peptide_paper_id not in imp_ids:
                    raise ValueError(
                        f"evidence references unknown immunizing_peptide_paper_id "
                        f"{ev.immunizing_peptide_paper_id!r}"
                    )
            elif ev.target_kind == "epitope":
                if ev.epitope_paper_id not in epitope_ids:
                    raise ValueError(
                        f"evidence references unknown epitope_paper_id "
                        f"{ev.epitope_paper_id!r}"
                    )
            elif ev.target_kind == "pool":  # (Q8)
                if ev.pool_paper_id not in pool_ids:
                    raise ValueError(
                        f"evidence references unknown pool_paper_id "
                        f"{ev.pool_paper_id!r}"
                    )
            else:  # candidate target (v2.8 P19)
                if ev.candidate_paper_id not in candidate_ids:
                    raise ValueError(
                        f"evidence references unknown candidate_paper_id "
                        f"{ev.candidate_paper_id!r}"
                    )

        # PATCH (v2 P1): pool membership now resolves to IMPs, not epitopes.
        for pool in self.pools:
            if pool.patient_paper_id not in patient_ids:
                raise ValueError(
                    f"pool {pool.paper_local_id!r} references unknown patient "
                    f"{pool.patient_paper_id!r}"
                )
            for mid in pool.member_peptide_ids:
                if mid not in imp_ids:
                    raise ValueError(
                        f"pool {pool.paper_local_id!r} references unknown member "
                        f"immunizing peptide {mid!r}"
                    )

        # PATCH (v2 P1): each epitope's parent long-peptides must resolve, and
        # an epitope's patient (via its IMP) must be consistent when set.
        for epi in self.epitopes:
            for pid in epi.parent_peptide_ids:
                if pid not in imp_ids:
                    raise ValueError(
                        f"epitope {epi.paper_local_id!r} references unknown "
                        f"parent immunizing peptide {pid!r}"
                    )

        # PATCH (v2 P1): immunizing peptides with a patient link must resolve.
        for imp in self.immunizing_peptides:
            if imp.patient_paper_id is not None and imp.patient_paper_id not in patient_ids:
                raise ValueError(
                    f"immunizing peptide {imp.paper_local_id!r} references unknown "
                    f"patient {imp.patient_paper_id!r}"
                )
        return self

    @model_validator(mode="after")
    def _peptide_counts_reconcile(self) -> ExtractedPaper:
        """PATCH (v2 P8): per-patient immunizing-peptide records must match the
        reported count. This is the manual Keskin check (distinct IMPs per
        patient == Supp Table 3b 'peptides in vaccine': 13/10/20/7/11/15/17/10)
        promoted to an invariant. Skipped for a patient with zero IMP records
        (peptides not extracted) so partial extractions still validate.
        """
        by_patient = Counter(
            i.patient_paper_id for i in self.immunizing_peptides if i.patient_paper_id
        )
        patients = {p.paper_local_id: p for p in self.patients}
        for pid, n_records in by_patient.items():
            pt = patients.get(pid)
            if pt is None:
                continue  # unresolved link already caught by _cross_reference_check
            expected = (
                pt.n_peptides_administered
                if pt.n_peptides_administered is not None
                else pt.n_peptides_synthesized
            )
            if n_records != expected:
                raise ValueError(
                    f"patient {pid!r}: {n_records} immunizing-peptide records but "
                    f"reported administered/synthesized count is {expected}"
                )
        return self

    @model_validator(mode="after")
    def _cohort_size_within_tolerance(self) -> ExtractedPaper:
        """Sum of patient cohort entries should be ~ paper-level cohort_size.

        Tolerance: ±5 patients OR 20% of cohort_size (whichever is larger).
        Real trials report enrollment + treated counts that differ by
        multiple patients (Rojas BNT122 PDAC: 19 enrolled, 16 treated).
        The validator is meant to catch the LLM truncating a cohort
        ("100 patients" but only 5 enumerated), not the natural
        enrollment-vs-treated gap.
        """
        if not self.patients:
            return self
        tolerance = max(5, self.cohort_size // 5)
        if abs(len(self.patients) - self.cohort_size) > tolerance:
            raise ValueError(
                f"len(patients)={len(self.patients)} differs from "
                f"cohort_size={self.cohort_size} by >{tolerance}"
            )
        return self
