#!/usr/bin/env python3
"""
gold_spec.py — the GOLD ANSWER-KEY contract for antVacDB extraction evals.

DELIBERATELY SEPARATE from schema.py. The gold key is ground truth a HUMAN
adjudicates from the source paper; it must not be coupled to (or derived from,
for grading purposes) any model's extraction, or you measure "does the agent
agree with the model that built the key" instead of "does the agent match
reality." A gold key is only valid for grading when `verified=True`, set by a
human who checked it against the source.

It captures the *checkable, countable* facts — not the full ExtractedPaper. Each
list item has a stable identity (used by score_extraction.py for set matching).
"""
from __future__ import annotations
from typing import Literal, Optional
from pydantic import BaseModel, Field


class GoldNeoantigen(BaseModel):
    patient: str                       # patient id as printed in the paper
    gene: str
    mutation: Optional[str] = None     # e.g. "G12V"; None for blank/frameshift-labelled
    immunogenic: bool                  # individually immunogenic (a single-target response)
    # identity = (patient, gene, norm(mutation))


class GoldEpitope(BaseModel):
    sequence: str
    mhc_class: Literal["I", "II"]
    allele: Optional[str] = None
    # identity = (norm(sequence), mhc_class)


class GoldMagnitude(BaseModel):
    patient: str
    gene: Optional[str] = None         # target neoantigen; None for pool/patient-level
    mutation: Optional[str] = None
    value: Optional[float] = None      # numeric magnitude, if the source gives a per-target number
    unit: Optional[str] = None         # sfc_per_1e6 / percent_of_parent / stimulation_index
    grade: Optional[str] = None        # negative/low/moderate/high/very_high (ordinal source)
    # identity = (patient, gene, norm(mutation))


class GoldSurvival(BaseModel):
    endpoint: str                      # rfs / os / landmark_rfs / ...
    arm: str                           # arm label
    median: Optional[float] = None     # months (None if not reached or not reported)
    not_reached: bool = False
    hazard_ratio: Optional[float] = None
    # identity = (endpoint, norm(arm))


class GoldKey(BaseModel):
    pmid: str
    paper: str = ""                    # human label, e.g. "Rojas 2023"
    source_refs: str = ""              # where the human read it (tables/figures), for audit
    # GRADING GATE — score_extraction refuses to grade against an unverified key.
    verified: bool = False
    verified_by: Optional[str] = None

    cohort_size: Optional[int] = None
    n_enrolled: Optional[int] = None
    # patient id -> count of individually immunogenic neoantigens
    per_patient_immunogenic: dict[str, int] = Field(default_factory=dict)

    neoantigens: list[GoldNeoantigen] = Field(default_factory=list)
    epitopes: list[GoldEpitope] = Field(default_factory=list)
    magnitudes: list[GoldMagnitude] = Field(default_factory=list)
    survival: list[GoldSurvival] = Field(default_factory=list)

    notes: str = ""

    model_config = {"extra": "forbid"}
