"""FLAGSHIP `neoantigen_tcr_flags` regression suite (2026-08-16).

Root cause it locks down: phase-5 wired the deterministic clonotype loaders (`add_reactive_tcr` /
`add_tcr_track`), which correctly removed the hand-building of `tcr_clonotypes` -- and the per-
(neoantigen x patient) TCR flags, a BY-PRODUCT of that hand-work, vanished with it (Keskin 3->0,
Hu 5->0, Rojas 1->0) because `_tcr_gateway_gap` is satisfied by tcr_seq_methods OR data_depositions
and NOTHING enforced the flag lane.

Covers: (A) the symmetric `_tcr_flag_gap` finalize gate + its HARD override, (B) loader-side flag
DERIVATION and — the fabrication boundary — that a flag is NEVER minted onto an unresolvable or
ambiguous target, (C) the _RULES prompt line that stops the loaders from displacing the flag lane.
"""
import json
import pathlib

import agent_core
import schema
from schema import ExtractedPaper

META = {"pmid": "30568305", "journal": "Nature", "year": 2019, "title": "flag gate test",
        "cohort_size": 1, "indication_summary": "GBM"}
PROV = {"quoted_text": "q", "section_ref": "Supplementary Table 8"}
CLONO = {"paper_local_id": "TCR1", "patient_paper_id": "P1", **PROV}
METHOD = {"method_local_id": "M1", "platform": "10x_5p_vdj", **PROV}
FLAG = {"target_kind": "immunizing_peptide", "immunizing_peptide_paper_id": "IMP1",
        "patient_paper_id": "P1", "tcr_identified": "cd4", **PROV}


def _paper(**kw):
    base = dict(
        patients=[{"paper_local_id": "P1", "indication": "GBM", "n_peptides_synthesized": 1,
                   "n_peptides_immunogenic": 1, **PROV}],
        immunizing_peptides=[{"paper_local_id": "IMP1", "sequence": "PQVDGEIPLHRSDRVKVLSIGEGGF",
                              "gene_symbol": "SHANK2", "is_neoantigen": True, **PROV}],
    )
    base.update(kw)
    return ExtractedPaper(**META, **base)


# ---------------------------------------------------------------- A. the gate predicate

def test_flag_gap_silent_without_a_tcr_claim():
    assert agent_core._tcr_flag_gap(_paper()) == ""                          # status None
    assert agent_core._tcr_flag_gap(_paper(tcr_seq_status="none")) == ""

def test_flag_gap_silent_for_mentioned_only():
    # 'mentioned_only' = TCR-seq named with nothing extractable; that is the GATEWAY's business.
    assert agent_core._tcr_flag_gap(
        _paper(tcr_seq_status="mentioned_only", tcr_seq_methods=[METHOD])) == ""

def test_flag_gap_silent_when_no_tcr_content_at_all():
    # empty everything is the gateway gap, not the flag gap -- do not double-block.
    assert agent_core._tcr_flag_gap(_paper(tcr_seq_status="both")) == ""

def test_flag_gap_fires_on_clonotypes_without_flags():
    gap = agent_core._tcr_flag_gap(_paper(tcr_seq_status="single_cell", tcr_clonotypes=[CLONO]))
    assert gap and "neoantigen_tcr_flags is EMPTY" in gap

def test_flag_gap_fires_on_methods_without_flags():
    assert agent_core._tcr_flag_gap(_paper(tcr_seq_status="bulk", tcr_seq_methods=[METHOD]))

def test_flag_gap_silent_when_flags_present():
    assert agent_core._tcr_flag_gap(
        _paper(tcr_seq_status="single_cell", tcr_clonotypes=[CLONO],
               neoantigen_tcr_flags=[FLAG])) == ""

def test_flag_override_is_hard_routes_to_needs_review():
    assert "allow_missing_tcr_flags" in agent_core.HARD_OVERRIDES
    assert agent_core.overrides_are_soft_only(["allow_missing_tcr_flags"]) is False


# ---------------------------------------------------------------- A. the gate in finalize

def _init(tmp_path, n_pep=0, **meta_extra):
    out = tmp_path / "r.json"
    agent_core.init_partial(str(out), json.dumps({**META, **meta_extra}))
    ok, msg = agent_core.append_section(str(out), "patients", json.dumps([
        {"paper_local_id": "P1", "indication": "GBM", "n_peptides_synthesized": n_pep,
         "n_peptides_immunogenic": 0, **PROV}]))
    assert ok, msg
    ok, msg = agent_core.append_section(str(out), "tcr_seq_methods", json.dumps([METHOD]))
    assert ok, msg
    return out

def test_finalize_blocks_on_missing_flags_and_says_how(tmp_path):
    out = _init(tmp_path, tcr_seq_status="single_cell")
    ok, msg = agent_core.finalize_partial(str(out))
    assert not ok and "TCR flags" in msg
    # the message must teach the SHAPE ...
    assert "neoantigen_tcr_flags" in msg and "tcr_identified" in msg
    for v in ("cd4", "cd8", "both", "none", "not_assessed"):
        assert v in msg
    # ... and, critically, HOW to get the loader-minted ids (schema._check_target rejects a guess).
    assert "partial_status" in msg

def test_finalize_passes_once_a_flag_exists(tmp_path):
    out = _init(tmp_path, n_pep=1, tcr_seq_status="single_cell")
    agent_core.append_section(str(out), "immunizing_peptides", json.dumps([
        {"paper_local_id": "IMP1", "sequence": "PQVDGEIPLHRSDRVKVLSIGEGGF", "gene_symbol": "SHANK2",
         "is_neoantigen": True, **PROV}]))
    ok, msg = agent_core.append_section(str(out), "neoantigen_tcr_flags", json.dumps([FLAG]))
    assert ok, msg
    ok, msg = agent_core.finalize_partial(str(out))
    assert ok, msg
    assert "allow_missing_tcr_flags" not in json.loads(pathlib.Path(out).read_text()).get(
        "finalize_overrides_used", [])

def test_finalize_override_routes_to_needs_review(tmp_path):
    out = _init(tmp_path, tcr_seq_status="single_cell")
    ok, msg = agent_core.finalize_partial(str(out), allow_missing_tcr_flags=True)
    assert ok, msg
    used = json.loads(pathlib.Path(out).read_text()).get("finalize_overrides_used", [])
    assert "allow_missing_tcr_flags" in used
    assert agent_core.overrides_are_soft_only(used) is False

def test_flags_lane_is_appendable_never_locked(tmp_path):
    # explicit: the loaders lock tcr_clonotypes; the flag lane must stay open for text-derived flags
    # (Rojas's one flag came from Fig. 2d, which no supplementary sheet can supply).
    out = _init(tmp_path, tcr_seq_status="single_cell")
    agent_core._lock_lanes(str(out), "tcr_clonotypes")
    assert "neoantigen_tcr_flags" not in agent_core._read_locks(str(out))
    agent_core.append_section(str(out), "immunizing_peptides", json.dumps([
        {"paper_local_id": "IMP1", "sequence": "PQVDGEIPLHRSDRVKVLSIGEGGF", "gene_symbol": "SHANK2",
         "is_neoantigen": True, **PROV}]))
    ok, msg = agent_core.append_section(str(out), "neoantigen_tcr_flags", json.dumps([FLAG]))
    assert ok, msg


def test_partial_status_hands_over_the_real_target_ids(tmp_path):
    out = _init(tmp_path, tcr_seq_status="single_cell")
    agent_core.append_section(str(out), "immunizing_peptides", json.dumps([
        {"paper_local_id": "IMP_PQVDGEIPLHRS_25_ab12cd", "sequence": "PQVDGEIPLHRSDRVKVLSIGEGGF",
         "gene_symbol": "SHANK2", "is_neoantigen": True, **PROV}]))
    ok, msg = agent_core.partial_status(str(out))
    assert ok and "IMP_PQVDGEIPLHRS_25_ab12cd" in msg and "SHANK2" in msg
    # and it stays quiet once the lane is filled (no noise on every status call)
    agent_core.append_section(str(out), "neoantigen_tcr_flags", json.dumps([
        {**FLAG, "immunizing_peptide_paper_id": "IMP_PQVDGEIPLHRS_25_ab12cd"}]))
    ok, msg = agent_core.partial_status(str(out))
    assert ok and "IMP_PQVDGEIPLHRS_25_ab12cd" not in msg


# ---------------------------------------------------------------- B. loader-side derivation

REC = {
    "patients": [{"paper_local_id": "Pt8"}, {"paper_local_id": "Pt7"}],
    "immunizing_peptides": [
        {"paper_local_id": "IMP_A", "gene_symbol": "SVEP1", "sequence": "GPPAHVENAIARGIHYQYGDMITYS"},
        {"paper_local_id": "IMP_B", "gene_symbol": "SHANK2", "sequence": "PQVDGEIPLHRSDRVKVLSIGEGGF"},
        {"paper_local_id": "IMP_C", "gene_symbol": "SHANK2", "sequence": "KPYQPQVDGEIPLHRSDRVKVLSIGEG"},
    ],
    "epitopes": [{"paper_local_id": "EPI_A", "gene_symbol": "ARHGAP35", "sequence": "MVNTVAGAMK"}],
    "neoantigen_tcr_flags": [],
}

def _derive(antigen, sheet, rec=None, patient="Pt8"):
    return agent_core.derive_reactive_tcr_flags(
        json.loads(json.dumps(rec if rec is not None else REC)), patient_paper_id=patient,
        antigen=antigen, sheet=sheet, section_ref=f"Supplementary Table 8 ({sheet})", n_clonotypes=39)

def test_derives_a_flag_when_sheet_gives_subset_and_a_unique_target():
    flags, why = _derive("mutant SVEP1", "(b)Pt8 mutSVEP1-reactive CD4+ Tcells")
    assert len(flags) == 1 and not why
    f = flags[0]
    assert f["target_kind"] == "immunizing_peptide" and f["immunizing_peptide_paper_id"] == "IMP_A"
    assert f["patient_paper_id"] == "Pt8" and f["tcr_identified"] == "cd4"
    # house rule: real provenance, derived from the sheet actually read -- never an invented figure cite.
    assert "mutSVEP1" in f["section_ref"] and "39 clonotype" in f["quoted_text"]
    assert f["needs_review"] is False
    schema.NeoantigenTcrFlag(**f)                       # the derived row is schema-valid as written

def test_derives_cd8_from_the_sheet_name():
    flags, _ = _derive("mutant SVEP1", "(c)Pt8 mutSVEP1-react CD8")
    assert flags[0]["tcr_identified"] == "cd8"

def test_epitope_target_resolves_when_no_peptide_matches():
    flags, _ = _derive("ARHGAP35 neoantigen CD8", "(c)Pt7 Neoantigen-reactive CD8", patient="Pt7")
    assert flags[0]["target_kind"] == "epitope" and flags[0]["epitope_paper_id"] == "EPI_A"

# --- the fabrication boundary: a flag is NEVER minted onto a guessed/ambiguous/absent target ---

def test_no_flag_when_the_gene_matches_several_peptides():
    # Keskin's real shape: SHANK2 has TWO overlapping vaccine peptides -> which one carried the TCR is
    # NOT determined by the sheet. Refuse, and name the candidates so the agent can finish it.
    flags, why = _derive("mutant SHANK2", "(a)Pt8 mutSHANK2-react CD4 Tcells")
    assert flags == [] and "ambiguous" in why and "IMP_B" in why and "IMP_C" in why

def test_no_flag_when_the_antigen_names_no_entity():
    flags, why = _derive("pool C neoantigen (CD4)", "(b)Pt7 Neoantigen-reactive CD4+", patient="Pt7")
    assert flags == [] and "matches no" in why

def test_no_flag_when_the_sheet_does_not_name_a_subset():
    # Keskin's Pt8 sheets carry no CD4/CD8 column or label; TcrIdentified has no 'found, subset unknown'.
    flags, why = _derive("mutant SVEP1", "(b)Pt8 mutSVEP1-react Tcells")
    assert flags == [] and "CD4" in why and "CD8" in why

def test_no_flag_when_the_sheet_names_both_subsets():
    flags, why = _derive("mutant SVEP1", "Pt8 CD4 and CD8 sorted")
    assert flags == [] and "BOTH" in why

def test_no_flag_for_an_unknown_patient():
    flags, why = _derive("mutant SVEP1", "Pt9 mutSVEP1 CD4", patient="Pt9")
    assert flags == [] and "not in the record" in why

def test_no_flag_duplicates_an_existing_one():
    rec = json.loads(json.dumps(REC))
    rec["neoantigen_tcr_flags"] = [{"target_kind": "immunizing_peptide",
                                    "immunizing_peptide_paper_id": "IMP_A",
                                    "patient_paper_id": "Pt8", "tcr_identified": "cd8", **PROV}]
    flags, why = _derive("mutant SVEP1", "Pt8 mutSVEP1 CD4", rec=rec)
    assert flags == [] and "already exists" in why

def test_every_derived_target_resolves_in_the_record():
    """The invariant the schema enforces (`_check_target`), asserted at the derivation site: whatever the
    sheet says, a minted flag's target id is one the record actually holds."""
    ids = {r["paper_local_id"] for lane in ("immunizing_peptides", "epitopes") for r in REC[lane]}
    for antigen, sheet in (("mutant SVEP1", "Pt8 CD4"), ("SVEP1", "Pt8 CD8 sorted"),
                           ("ARHGAP35", "Pt8 CD4"), ("mutSVEP1", "Pt8 CD4"),
                           ("mutant SHANK2", "Pt8 CD4"), ("pool C", "Pt8 CD4"),
                           ("GPPAHVENAIARGIHYQYGDMITYS", "Pt8 CD4"),
                           ("totally unknown antigen", "Pt8 CD4"), ("", "Pt8 CD4")):
        flags, _ = _derive(antigen, sheet)
        for f in flags:
            tid = (f.get("immunizing_peptide_paper_id") or f.get("epitope_paper_id")
                   or f.get("candidate_paper_id"))
            assert tid in ids, (antigen, sheet, tid)


# ---------------------------------------------------------------- C. prompt guidance

def test_rules_block_says_loaders_do_not_fill_the_flag_lane():
    import prompt_render
    txt = prompt_render._RULES
    i = txt.find("add_tcr_track")
    assert i > 0
    tail = txt[i:]
    assert "neoantigen_tcr_flags" in tail and "do NOT populate" in tail
    assert "FLAGSHIP" in tail
