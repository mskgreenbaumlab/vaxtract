#!/usr/bin/env python3
"""antVacDB pre-deploy gate (RUNBOOK Loop 3 / production checklist).

Runs the two halves the RUNBOOK requires before any schema deploy, plus an
integrity check:

  PART 1  BACK-COMPAT       every record in the corpus still validates under the
                            current schema  ("nothing that used to pass now fails")
  PART 2  GUARD REGRESSION  a battery of known-bad inputs that MUST be rejected
                            ("the schema still catches what it should")
  PART 3  INTEGRITY         module imports (runs all vocab-lockstep asserts) and
                            stamps SCHEMA_VERSION

Exit code 0 iff every corpus record validates AND every guard fires. Usable in CI.

Usage:  python3 predeploy_gate.py [CORPUS_DIR]   (default: reference_records)

NOTE: the corpus here is the fixture set in CORPUS_DIR. In production, point
CORPUS_DIR at (a snapshot/export of) the full historical store so the
shadow-validation is against everything ever ingested, not just the fixtures.
"""
from __future__ import annotations
import sys, json, glob, os, importlib

CORPUS_DIR = sys.argv[1] if len(sys.argv) > 1 else "reference_records"
sys.path.insert(0, CORPUS_DIR)   # import the DEPLOYED schema/vocab (source of truth)
import vocab, schema  # noqa
from schema import (
    ExtractedPaper, ExtractedPatient, ImmunizingPeptide, MinimalEpitope,
    ExtractedEvidence, Measurement, ResponseMagnitude, SurvivalOutcome,
    PreclinicalEfficacy, CuratorNote, ClinicalBenefitSignal,
    NeoantigenMutation, VafPoint,
)

fail = 0
def ok(msg):  print(f"  \033[32mPASS\033[0m  {msg}")
def bad(msg):
    global fail; fail += 1; print(f"  \033[31mFAIL\033[0m  {msg}")

# --- minimal valid building blocks (mutated below to trigger guards) ----------
def imp(pid="i1", patient="P1"):
    return ImmunizingPeptide(paper_local_id=pid, sequence="SLLQHLIGL", is_neoantigen=True,
        patient_paper_id=patient, quoted_text="q", section_ref="s")
def patient(**kw):
    base = dict(paper_local_id="P1", indication="X", n_peptides_synthesized=2,
                n_peptides_immunogenic=1, quoted_text="q", section_ref="s")
    base.update(kw); return ExtractedPatient(**base)
def evidence(**kw):
    base = dict(patient_paper_id="P1", target_kind="immunizing_peptide",
                immunizing_peptide_paper_id="i1", assay="elispot", outcome="immunogenic",
                quoted_text="q", section_ref="s")
    base.update(kw); return ExtractedEvidence(**base)
def paper(**kw):
    base = dict(pmid="1", journal="J", year=2023, title="t", cohort_size=1,
                indication_summary="x", patients=[patient()])
    base.update(kw); return ExtractedPaper(**base)

def must_reject(label, fn):
    try:
        fn(); bad(f"{label}  (constructed without error — guard did NOT fire)")
    except Exception:
        ok(f"{label}  -> rejected")

# =============================================================================
print(f"\n=== antVacDB PRE-DEPLOY GATE ===  corpus={CORPUS_DIR}")
print(f"=== PART 3  INTEGRITY ===")
naxes = open(os.path.join(CORPUS_DIR, "schema.py")).read().count("_assert_vocab(") - 1
print(f"  schema {schema.SCHEMA_VERSION} imported clean; {naxes} vocab-lockstep axes pass")

print(f"\n=== PART 1  BACK-COMPAT  (every corpus record must validate) ===")
files = sorted(glob.glob(os.path.join(CORPUS_DIR, "*extracted*.json")))
if not files: bad("no corpus records found")
for f in files:
    try:
        p = ExtractedPaper(**json.loads(open(f).read()))
        ok(f"{os.path.basename(f):34} patients={len(p.patients)} epitopes={len(p.epitopes)} "
           f"evidence={len(p.evidence)} survival={len(p.survival_outcomes)}")
    except Exception as e:
        bad(f"{os.path.basename(f)} :: {str(e)[:90]}")

print(f"\n=== PART 2  GUARD REGRESSION  (known-bad inputs must be rejected) ===")
print("  -- controlled-vocab (Literal) axes --")
must_reject("off-vocab evidence outcome",        lambda: evidence(outcome="bogus"))
must_reject("off-vocab evidence assay",          lambda: evidence(assay="qpcr"))
must_reject("off-vocab survival endpoint",       lambda: SurvivalOutcome(endpoint="ttx"))
must_reject("off-vocab response grade",          lambda: ResponseMagnitude(grade="huge", raw="x"))
must_reject("off-vocab response unit",           lambda: ResponseMagnitude(value=1, unit="spots", raw="1"))
print("  -- identity / PHI guards --")
must_reject("name-like patient paper_local_id",  lambda: patient(paper_local_id="John Smith"))
print("  -- peptide / epitope guards --")
must_reject("class-I epitope missing affinity slot", lambda: MinimalEpitope(
    paper_local_id="e", sequence="CRLQKVAAL", is_neoantigen=True, quoted_text="q",
    section_ref="s", mhc_class="I"))
must_reject("class-II epitope carrying nM affinity", lambda: MinimalEpitope(
    paper_local_id="e", sequence="ELAGIGILTV", is_neoantigen=True, quoted_text="q",
    section_ref="s", mhc_class="II", predicted_affinity=Measurement(unit="nM", value=5.0, raw="5")))
print("  -- evidence guards --")
must_reject("evidence with two targets set",     lambda: evidence(epitope_paper_id="e1"))
must_reject("pre_existing + vaccine_induced=True", lambda: evidence(outcome="pre_existing", vaccine_induced=True))
must_reject("positive outcome + grade negative", lambda: evidence(
    outcome="immunogenic", magnitude=ResponseMagnitude(grade="negative", raw="-")))
print("  -- Measurement / ResponseMagnitude value guards --")
must_reject("Measurement value + unit unknown",  lambda: Measurement(unit="unknown", value=1.0, raw="1"))
must_reject("magnitude empty (no value/grade/raw)", lambda: ResponseMagnitude())
must_reject("magnitude value + unit unknown",    lambda: ResponseMagnitude(value=1.0, unit="unknown", raw="1"))
must_reject("magnitude percent_of_parent > 100", lambda: ResponseMagnitude(value=150, unit="percent_of_parent", raw="150%"))
print("  -- SurvivalOutcome guards --")
must_reject("survival not_reached + median_value", lambda: SurvivalOutcome(endpoint="rfs", not_reached=True, median_value=13.4))
must_reject("survival hazard_ratio w/o comparator", lambda: SurvivalOutcome(endpoint="rfs", hazard_ratio=0.5))
must_reject("survival p_value > 1",              lambda: SurvivalOutcome(endpoint="rfs", p_value=1.5))
print("  -- PreclinicalEfficacy guard --")
must_reject("efficacy prolonged_survival on growth readout", lambda: PreclinicalEfficacy(
    readout="tumor_growth", result="prolonged_survival"))
must_reject("efficacy size-result on survival readout", lambda: PreclinicalEfficacy(
    readout="survival", result="growth_inhibition"))
print("  -- ExtractedPaper integrity guards --")
must_reject("duplicate immunizing_peptide paper_local_id", lambda: paper(
    immunizing_peptides=[imp("i1"), imp("i1")], patients=[patient(n_peptides_synthesized=2, n_peptides_administered=2)]))
must_reject("count reconciliation mismatch (3 IMP vs administered=2)", lambda: paper(
    patients=[patient(n_peptides_synthesized=3, n_peptides_administered=2)],
    immunizing_peptides=[imp("i1"), imp("i2"), imp("i3")]))
must_reject("cohort size far out of tolerance",  lambda: paper(cohort_size=100, patients=[patient()]))
must_reject("evidence referencing unknown patient", lambda: paper(
    immunizing_peptides=[imp("i1")], evidence=[evidence(patient_paper_id="P9")]))
print("  -- trial_setting + ClinicalBenefitSignal guards (v2.10 P22) --")
must_reject("patient off-vocab trial_setting", lambda: patient(trial_setting="phase_2"))
must_reject("benefit signal off-vocab readout", lambda: ClinicalBenefitSignal(readout="vibes", direction="increased"))
must_reject("benefit 'cleared' on antigen_loss readout", lambda: ClinicalBenefitSignal(
    readout="antigen_loss", direction="cleared"))
print("  -- evidence-count anchor guards (v2.11.1 P20.1) --")
must_reject("negative n_immunogenic_reported", lambda: paper(n_immunogenic_reported=-1))
must_reject("negative n_tested_negative_reported", lambda: paper(n_tested_negative_reported=-3))
print("  -- NeoantigenMutation guards (v2.11 P20) --")
must_reject("mutation with neither gene nor genomic_change", lambda: NeoantigenMutation(
    paper_local_id="M1", patient_paper_id="P1", quoted_text="q", section_ref="s"))
must_reject("mutation off-vocab clonality", lambda: NeoantigenMutation(
    paper_local_id="M1", patient_paper_id="P1", gene_symbol="RLF", clonality="branchy",
    quoted_text="q", section_ref="s"))
must_reject("mutation off-vocab status", lambda: NeoantigenMutation(
    paper_local_id="M1", patient_paper_id="P1", gene_symbol="RLF", status="vanished",
    quoted_text="q", section_ref="s"))
must_reject("VafPoint empty (no value/raw)", lambda: VafPoint())
must_reject("VafPoint value > 1",            lambda: VafPoint(value=1.5))
must_reject("peptide cancer_cell_fraction > 1", lambda: ImmunizingPeptide(
    paper_local_id="i1", sequence="SLLQHLIGL", is_neoantigen=True, cancer_cell_fraction=1.5,
    quoted_text="q", section_ref="s"))
print("  -- CuratorNote guards (v2.7) --")
must_reject("curator note off-vocab kind",       lambda: CuratorNote(kind="win", text="x"))
must_reject("curator note empty text",           lambda: CuratorNote(kind="decision", text=""))
must_reject("curator note dangling ref",         lambda: paper(
    curator_notes=[CuratorNote(kind="highlight", text="re NOPE", refs=["NOPE"])]))

# also confirm a fully-valid loaded record with the NEW fields round-trips (sanity)
print("  -- positive control (valid record with new fields) --")
try:
    paper(patients=[patient(n_peptides_synthesized=1, trial_setting="adjuvant",
                clinical_benefit_signals=[ClinicalBenefitSignal(
                    readout="ctdna_dynamics", direction="cleared", associated_with_response=True)])],
          immunizing_peptides=[imp("i1")],
          evidence=[evidence(magnitude=ResponseMagnitude(value=2459, unit="sfc_per_1e6", raw="2459"))],
          survival_outcomes=[SurvivalOutcome(endpoint="rfs", arm_label="responders", not_reached=True,
                comparator_label="non-responders", hazard_ratio=0.08, p_value=0.003)],
          clinical_benefit_signals=[ClinicalBenefitSignal(readout="tumor_infiltration", direction="detected")],
          neoantigen_mutations=[NeoantigenMutation(paper_local_id="MUT1", patient_paper_id="P1",
                gene_symbol="RLF", genomic_change="chr1:40705080 A>T", status="emerged", clonality="subclonal",
                hla_restrictions=["HLA-A*11:01"], vaf=[VafPoint(timepoint_label="primary", value=0.0),
                VafPoint(timepoint_label="recurrent", value=0.166)], quoted_text="RLF", section_ref="Table S6")],
          curator_notes=[CuratorNote(kind="decision", text="example note", refs=["P1"])])
    ok("valid record w/ trial_setting + benefit_signal + neoantigen_mutation + magnitude + survival + curator_note constructs")
except Exception as e:
    bad(f"valid record rejected :: {str(e)[:90]}")

# =============================================================================
print(f"\n=== GATE RESULT ===")
if fail == 0:
    print("  \033[32mALL CHECKS PASSED — safe to deploy schema "
          f"{schema.SCHEMA_VERSION}\033[0m")
    sys.exit(0)
else:
    print(f"  \033[31m{fail} CHECK(S) FAILED — DO NOT DEPLOY\033[0m")
    sys.exit(1)
