"""Single source of truth for antVac controlled vocabularies.

Both `schema.py` (Layer-1 validation) and the Layer-2 prompt builder
(`build_gold_outcomes.py`) MUST draw their allowed outcome/assay values from
here, so the prompt and the schema can never silently drift. See
SCHEMA.md "Vocabulary reconciliation (QV)".

How to use:
  - schema.py:  `EvidenceOutcome = Literal[*OUTCOMES]`  (etc.)
  - build_gold_outcomes.py: inject `", ".join(OUTCOMES)` into the protocol prompt
    instead of a hand-written `immunogenic|not_immunogenic|...` string, and
    replace the open-ended assay `...` with `", ".join(ASSAYS)`.
  - tests assert the protocol's allowed values == these tuples.

DECISIONS BAKED IN (resolved 2026-05; see HANDOFF.md decisions-log):
  - outcome spelling is `not_immunogenic` (underscore) — matches the protocol
    and all existing Layer-2 data; the old schema `non-immunogenic` (hyphen) was
    the outlier and is removed.
  - `flow_cytometry` is a first-class assay (4-1BB/OX40 activation-marker flow is
    distinct from `ics` intracellular cytokine staining).
  - `tcr_reporter` added (Jurkat/reporter-line IL-2 specificity assay).
"""
from __future__ import annotations

# Per-peptide / per-measurement outcome vocabulary.
OUTCOMES: tuple[str, ...] = (
    "immunogenic",       # vaccine-induced T-cell response detected (named in a quote)
    "not_immunogenic",   # tested, no vaccine-induced response (incl. mut==WT)
    "pre_existing",      # reactivity present PRE-vaccine (vaccine_induced=false)
    "not_assessed",      # included but immunogenicity not reported
    "positive",          # generic positive (binding / MS hit)
    "negative",          # generic negative
    "presented",         # MS-eluted / surface-presented
)

# Assay vocabulary. CLOSED set — the Layer-2 prompt must list exactly these,
# not an open-ended "...". Tier hints retained from the original schema.
ASSAYS: tuple[str, ...] = (
    "elispot",           # IFN-gamma ELISpot
    "tetramer",          # tetramer / dextramer
    "ics",               # intracellular cytokine staining
    "flow_cytometry",    # 4-1BB / OX40 activation-marker flow (NOT ics)
    "ms",                # MS-immunopeptidomics
    "binding",           # MHC binding affinity / stability
    "tcr_reporter",      # Jurkat / reporter-line IL-2 specificity assay
    "vaccine_response",  # post-vaccination T-cell response (paper-level)
    "clinical_response", # vaccine-induced clinical benefit
)

# Granularity of an evidence row's target. RECONCILED 2026-06 with schema v2.2:
# the single "peptide" target split into the two real levels — a long
# "immunizing_peptide" (what a pool contains; CD4 / long-peptide responses) and
# a minimal "epitope" (class-I; CD8 responses). Was ("peptide", "pool"); the
# schema changed and this tuple did not — the drift that motivated guarding ALL
# vocab axes (see schema.py "VOCAB LOCKSTEP").
# v2.8 (P19): "candidate" added so a screened-but-not-selected NeoantigenCandidate
# can carry its own outcome evidence (the funnel's rejected arm).
EVIDENCE_TARGETS: tuple[str, ...] = ("immunizing_peptide", "epitope", "pool", "candidate")

# Per-target outcome in a SCREENING manifest (schema v2.14, #4). A prediction/target manifest enumerates
# every screened target with a bulk readout column (e.g. Rojas 232-target ELISpot: 200 "No response").
# Those are STORED as ScreeningReadout rows (NOT immunogenicity ExtractedEvidence) so the stated facts are
# kept without inflating the immunogenicity axis. MUST equal the schema Literal (VOCAB LOCKSTEP).
MANIFEST_OUTCOMES: tuple[str, ...] = ("response", "no_response", "not_evaluable")

# T-cell subset + MHC class on an evidence row (schema v2.12). Two ORTHOGONAL axes added to
# ExtractedEvidence: which T-cell compartment responded (CD4 vs CD8) and which MHC class restricted
# the response. Both are OPTIONAL (None = the field was never populated / un-migrated); the processed-
# but-no-cue value is an explicit member ('bulk_or_unknown' / 'not_determined'), kept DISTINCT from
# None so a back-filled NULL is never confused with a deliberate "no cue found" call. A subset/class is
# NEVER guessed: it is set only when the paper states a cue (verbatim token enforced in quoted_text).
# MUST equal their schema Literals (VOCAB LOCKSTEP).
T_CELL_SUBSETS: tuple[str, ...] = ("cd4", "cd8", "bulk_or_unknown")
EVIDENCE_MHC_CLASSES: tuple[str, ...] = ("class_i", "class_ii", "not_determined")

ASSAY_STIMULATION: tuple[str, ...] = ("ex_vivo", "in_vitro_expanded", "unspecified")
ASSAY_ANTIGEN_FORMAT: tuple[str, ...] = (
    "peptide_pulsed", "minigene", "autologous_tumor", "tetramer", "unspecified",
)
TIMEPOINT_PHASE: tuple[str, ...] = (
    "pre_vaccine", "on_treatment", "post_vaccine", "memory", "unspecified",
)

# Layer-2 immunogenicity-scoring subset.
# build_gold_outcomes.py renders these into layer2_protocol.md placeholders.
# Per (patient, peptide) outcome scoring uses ONLY these 4 — the other 3
# OUTCOMES values (positive/negative/presented) belong to Layer-1 evidence
# (binding / MS / generic), not vaccine-induced T-cell response scoring.
LAYER2_OUTCOMES: tuple[str, ...] = (
    "immunogenic",
    "not_immunogenic",
    "pre_existing",
    "not_assessed",
)

# Layer-2 assay subset: T-cell functional readouts only.
# Excluded from ASSAYS: ms (presentation, not response), binding (affinity),
# vaccine_response (paper-level summary), clinical_response (clinical benefit).
LAYER2_ASSAYS: tuple[str, ...] = (
    "elispot",
    "tetramer",
    "ics",
    "flow_cytometry",
    "tcr_reporter",
)


# ===========================================================================
# Structural / value vocabularies added with schema v2.x. RECONCILED 2026-06:
# these were introduced directly in schema.py and had no single source here —
# now centralized so the Layer-2 prompt and the schema cannot drift on them.
# Each MUST equal its schema Literal (enforced by schema.py "VOCAB LOCKSTEP").
# ===========================================================================

# Peptide / epitope structural axes (schema v2.0 two-level model).
# MHC class governs which fields are meaningful: class-I minimal epitopes carry
# a NetMHCpan nM/%rank affinity; class-II / CD4 epitopes do not.
MHC_CLASSES: tuple[str, ...] = ("I", "II")

# Variant class of a neoantigen. Frameshift vs missense matters mechanistically
# (and sets the peptide-length expectation); kept as the cross-paper axis
# alongside the free-form `mutation` and HGVS-ish `protein_change`.
VARIANT_TYPES: tuple[str, ...] = (
    "missense", "frameshift", "inframe_indel", "stop_loss", "splice", "fusion", "other",
)

# Measurement (unit-explicit, lossless value-claim) axes (schema v2.1).
# A parsed number must carry a known unit; "unknown" is the lossless escape
# hatch for a captured-but-unparsed raw token (never silently coerced).
MEASUREMENT_UNITS: tuple[str, ...] = ("nM", "pct_rank", "unknown")
# A prediction is not reported data is not a validated wet-lab readout.
MEASUREMENT_TIERS: tuple[str, ...] = ("predicted", "reported", "validated")

# Structured provenance source kind (schema v2 P7). 'lookup' = deterministic
# fetch (e.g. PubMed), distinct from text/table/figure read from the paper.
PROVENANCE_KINDS: tuple[str, ...] = ("table", "figure", "prose", "lookup")

# Concomitant immunosuppression (schema v2.2) — a known immunogenicity
# confounder, stored as FACT not conclusion (Keskin: corticosteroid during
# priming vs n_peptides_immunogenic==0). agent_class is the comparable axis.
IMMUNOSUPPRESSANT_CLASSES: tuple[str, ...] = (
    "corticosteroid", "chemotherapy", "calcineurin_inhibitor", "mtor_inhibitor",
    "antimetabolite", "anti_tnf", "other_immunosuppressant", "unspecified",
)
IMMUNOSUPPRESSION_TIMING: tuple[str, ...] = (
    "pre_vaccine", "during_priming", "during_boost", "throughout",
    "post_vaccine", "unspecified",
)

# Concomitant systemic anti-cancer therapy (schema v2.15 — trial-context axis #1). The
# CO-ADMINISTERED anti-cancer treatment a vaccinated patient received (checkpoint inhibitor /
# chemo / RT / targeted), recorded as FACT — never the paper's efficacy attribution. Distinct
# from ConcomitantImmunosuppression (a confounder working AGAINST the vaccine) and from
# PreclinicalEfficacy.combination (a preclinical ARM label). Near-universal in modern trials
# (Rojas = autogene cevumeran + atezolizumab + mFOLFIRINOX). `drug_class` is the comparable axis.
# MUST equal their schema Literals (VOCAB LOCKSTEP).
CONCOMITANT_DRUG_CLASSES: tuple[str, ...] = (
    "checkpoint_inhibitor", "chemotherapy", "radiotherapy", "targeted", "other",
)
CONCOMITANT_THERAPY_TIMING: tuple[str, ...] = ("concurrent", "sequential", "unknown")

# Tumour mutational burden / microsatellite status (schema v2.15 — trial-context axis #3).
# Per-patient genomic context that gates neoantigen-load interpretation. None on the field =
# not reported; 'unknown' is a STATED-but-uncharacterized status (not a not-reported escape).
# MUST equal its schema Literal (VOCAB LOCKSTEP).
MSI_STATUSES: tuple[str, ...] = ("mss", "msi_high", "unknown")


# Subject species + vaccine platform (schema v2.3). The schema was human-clinical
# by construction; these let preclinical (esp. murine) DNA-vaccine cohorts be
# represented. MUST equal their schema Literals (VOCAB LOCKSTEP).
SPECIES: tuple[str, ...] = ("human", "mouse", "rat", "non_human_primate", "other")
# Kind of cohort an ExtractedPatient represents (schema v2.13 A1). Separates real vaccinated DISEASE
# cohorts (which count toward "vaccinated patients" / recall denominators) from METHODOLOGICAL arms
# (e.g. an HLA-A2-transgenic model-antigen pipeline-validation set using viral/gp100 epitopes, Li
# 33879241 P1) so the latter is queryably EXCLUDABLE instead of silently dropped. None = un-tagged
# (legacy/back-compat). MUST equal the schema Literal (VOCAB LOCKSTEP).
COHORT_KINDS: tuple[str, ...] = (
    "patient", "tumor_model", "model_antigen_validation", "healthy_donor", "other",
)
# Paper-level DATA RESOLUTION (schema v2.16): the finest immunogenicity grain a study REPORTS at,
# ordered FINEST -> COARSEST. Lets a genuinely coarse paper (gene+mutation or cohort-summary only,
# e.g. 39762422 autogene cevumeran) be admitted faithfully + tagged instead of held as under-
# extracted by the per-sequence recall anchor. None = un-tagged (legacy/back-compat). MUST equal the
# schema Literal (VOCAB LOCKSTEP).
DATA_RESOLUTIONS: tuple[str, ...] = (
    "per_sequence", "per_mutation", "per_target_gene", "cohort_summary", "clinical_only",
)
VACCINE_PLATFORMS: tuple[str, ...] = (
    "dna", "rna", "synthetic_long_peptide", "short_peptide",
    "dendritic_cell", "viral_vector", "other", "unspecified",
)

# Vaccine DELIVERY covariates (schema v2.9 P21). adjuvant = a SEPARATELY-ADDED
# immunostimulant ONLY, named as papers report it (NOT by mechanism — derive TLR class
# later via a name->TLR lookup). An LNP/lipoplex is a FORMULATION, not an adjuvant; an RNA
# vaccine takes adjuvant='none'. dose_basis tags how a dose number is expressed so it is
# compared only within-basis. MUST equal their schema Literals (VOCAB LOCKSTEP).
ADJUVANT_CLASSES: tuple[str, ...] = (
    "poly_iclc", "montanide", "gm_csf", "cpg", "none", "other", "unspecified",
)
DOSE_SCHEDULE_BASES: tuple[str, ...] = ("per_peptide", "per_pool", "total", "unspecified")

# Clinical trial setting (schema v2.10 P22). The DISEASE/treatment context a patient was enrolled
# in — the cross-trial covariate that gates whether a survival/response endpoint is even comparable
# (an adjuvant RFS is not a metastatic ORR; the Nat Med 2026 roadmap names early-stage/(neo)adjuvant
# focus an explicit priority). PER-PATIENT grain (like vaccine_platform) so a combined-cohort paper
# (e.g. metastatic melanoma + adjuvant PDAC) is captured per subject, not collapsed to one paper
# label. None = not reported, NEVER inferred; 'other' is a STATED-but-unlisted setting (not a
# not-reported escape). MUST equal its schema Literal (VOCAB LOCKSTEP).
CLINICAL_TRIAL_SETTINGS: tuple[str, ...] = (
    "adjuvant",          # post-resection, no measurable disease (Rojas PDAC)
    "neoadjuvant",       # pre-resection
    "perioperative",     # spanning surgery (neoadjuvant + adjuvant)
    "metastatic",        # measurable metastatic disease
    "locally_advanced",  # unresectable locoregional disease
    "recurrent",         # treated at recurrence/relapse (some Keskin GBM)
    "other",
)

# Cohort-level delivery/treatment latency metric (schema v2.9.1 P21). The paper-level home for an
# unlabelled per-patient latency distribution (Rojas Fig 1D surgery->vaccine). MUST equal its
# schema Literal (VOCAB LOCKSTEP).
LATENCY_METRICS: tuple[str, ...] = (
    "surgery_to_first_vaccine", "surgery_to_first_treatment",
    "diagnosis_to_first_vaccine", "other",
)


# Preclinical antitumor efficacy (schema v2.4). The preclinical analog of the
# RANO clinical_outcome axis; describes one efficacy experiment ARM as a fact.
# MUST equal their schema Literals (VOCAB LOCKSTEP).
EFFICACY_READOUTS: tuple[str, ...] = (
    "tumor_growth", "survival", "tumor_free_fraction", "metastasis", "other",
)
EFFICACY_RESULTS: tuple[str, ...] = (
    "no_effect", "partial_inhibition", "growth_inhibition",
    "tumor_regression", "complete_response", "prolonged_survival", "not_reported",
)
EFFICACY_SETTINGS: tuple[str, ...] = ("prophylactic", "therapeutic", "unspecified")
COMBINATION_CLASSES: tuple[str, ...] = (
    "monotherapy", "plus_checkpoint_inhibitor", "plus_chemotherapy",
    "plus_radiotherapy", "plus_other", "unspecified",
)


# Clinical time-to-event outcome (schema v2.5). The HUMAN analog of
# PreclinicalEfficacy: clinical_outcome is RANO (CR/PR/SD/PD), which is
# undefined for an ADJUVANT trial with no measurable disease (Rojas et al. 2023
# autogene cevumeran PDAC — primary result is recurrence-free survival, not a
# RANO response). `SURVIVAL_ENDPOINTS` is the cross-paper-comparable axis; the
# median/HR/P live as parsed values on the SurvivalOutcome object.
# MUST equal their schema Literals (VOCAB LOCKSTEP).
SURVIVAL_ENDPOINTS: tuple[str, ...] = (
    "rfs",          # recurrence-free survival
    "os",           # overall survival
    "landmark_rfs", # RFS measured from a landmark (e.g. last priming dose)
    "dfs",          # disease-free survival
    "pfs",          # progression-free survival
    "efs",          # event-free survival
    "ttr",          # time to recurrence
    "other",
)
SURVIVAL_TIME_UNITS: tuple[str, ...] = ("months", "weeks", "days", "years")


# Response-magnitude axes (schema v2.6). The structured home for ELISpot/ICS/
# tetramer magnitude that previously fell into evidence.assay_detail free text.
# A SIBLING of the Measurement (affinity) axes, deliberately NOT merged: magnitude
# units are a different quantity from binding-affinity units, and an ordinal
# +/++/+++ grade cannot live in a float Measurement.value. MUST equal their
# schema Literals (VOCAB LOCKSTEP).
RESPONSE_MAGNITUDE_UNITS: tuple[str, ...] = (
    "sfc_per_1e6",        # spot-forming cells/units per 1e6 cells (ELISpot); raw preserves SFC vs SFU
    "percent_of_parent",  # % of a gated parent population (tetramer+/cytokine+/marker+ of CD8+ etc.)
    "stimulation_index",  # ratio stimulated:unstimulated
    "unknown",            # lossless escape hatch (raw kept, value unparsed/unmapped)
)
# Ordinal semiquantitative grade. Maps the murine vaccine-table scale
# - / + / ++ / +++ / ++++  ->  negative / low / moderate / high / very_high.
RESPONSE_GRADES: tuple[str, ...] = ("negative", "low", "moderate", "high", "very_high")


# Curator-note kinds (schema v2.7). A PASSIVE container for human/curatorial
# commentary on a paper (challenges, decisions-under-uncertainty, caveats,
# highlights) — NOT an extracted fact and NOT model-generated at render time. The
# kind axis lets the report group/colour notes and lets analysis filter them.
# MUST equal its schema Literal (VOCAB LOCKSTEP).
CURATOR_NOTE_KINDS: tuple[str, ...] = ("challenge", "decision", "caveat", "highlight")


# NeoantigenCandidate funnel stage (schema v2.8 P19). A candidate is the same
# peptide-shaped object as an ImmunizingPeptide but recorded at an earlier funnel
# stage: predicted (in the ranked list) -> selected (chosen for the vaccine) ->
# administered (actually given). MUST equal its schema Literal (VOCAB LOCKSTEP).
CANDIDATE_STATUSES: tuple[str, ...] = ("predicted", "selected", "administered")

# PrioritizationScore kinds (schema v2.8 P19). The cross-paper-comparable metric
# TYPE of a funnel ranking score — a DIFFERENT quantity from binding affinity, so
# it lives on its own sibling (PrioritizationScore), not on Measurement. Covers
# the metrics real pipelines rank on (NetMHCpan affinity, expression TPM, pipeline
# quality score, pVACtools rank, agretopicity/DAI, VAF, clonality, foreignness).
# MUST equal its schema Literal (VOCAB LOCKSTEP).
SCORE_KINDS: tuple[str, ...] = (
    "affinity", "expression_tpm", "quality_score", "rank",
    "agretopicity_dai", "vaf", "clonality", "foreignness", "other",
)


# Clinical-benefit-signal axes (schema v2.10 P22). The response->clinical-benefit BRIDGE — the
# "better proxies to estimate vaccine efficacy" the Nat Med 2026 roadmap asks for. These are the
# benefit events that have NO existing home: epitope spreading, ctDNA dynamics, intratumoral T-cell
# infiltration, antigen loss at relapse. NOT ExtractedEvidence (its exactly-one-target validator
# requires a peptide; these are not peptide-specific — epitope spreading is BY DEFINITION reactivity
# to antigens NOT in the vaccine), NOT SurvivalOutcome (not a time-to-event statistic), NOT
# clinical_outcome (RANO tumour-burden category only). Recorded as FACT, never the paper's causal
# story. `antigen_loss` deliberately doubles as the home for the deferred P20 clonality doc's
# antigen-loss-at-relapse event. MUST equal their schema Literals (VOCAB LOCKSTEP).
BENEFIT_READOUTS: tuple[str, ...] = (
    "epitope_spreading",   # reactivity broadened to non-vaccine antigens
    "ctdna_dynamics",      # circulating tumour DNA trajectory (clearance/recurrence)
    "tumor_infiltration",  # T-cell infiltration into the tumour (e.g. intratumoral CD8)
    "antigen_loss",        # targeted antigen/HLA lost at relapse (P20 bridge)
    "other",
)
# Direction/result of a benefit signal. Cross-readout because the SAME verbs describe most signals
# (ctDNA cleared, infiltration increased, antigen lost); a readout<->direction sanity guard lives on
# the model. MUST equal its schema Literal (VOCAB LOCKSTEP).
BENEFIT_DIRECTIONS: tuple[str, ...] = (
    "increased", "decreased", "cleared", "persisted",
    "lost", "detected", "not_detected", "unchanged",
)


# Clonality / antigen-dynamics axes (schema v2.11 P20). The literature's primary axis is
# clonal-vs-subclonal (TRACERx, Nat 2023); whole-genome-doubling timing is the secondary,
# mutation-relative refinement (Swanton "prioritize pre-WGD neoantigens"). `status` records a
# mutation's antigen DYNAMICS across timepoints — `emerged` (primary VAF 0 -> recurrent >0, e.g.
# 34903219 Table S6 'new neoantigen mutations in recurrent tumor') and `lost` (the Swanton relapse-
# loss failure mode) are the two that matter. Default is ALWAYS unknown/None, NEVER clonal — assuming
# clonality is exactly the error that sinks a vaccine. MUST equal their schema Literals (VOCAB LOCKSTEP).
CLONALITY_CLASSES: tuple[str, ...] = ("clonal", "subclonal", "unknown")
WGD_TIMINGS: tuple[str, ...] = ("pre_wgd", "post_wgd", "unknown")
ANTIGEN_STATUSES: tuple[str, ...] = ("present", "emerged", "lost", "retained", "unknown")
