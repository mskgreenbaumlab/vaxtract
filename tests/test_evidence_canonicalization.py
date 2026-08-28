"""Evidence-lane determinism (stability lab, 2026-07-06): the three finalize passes that take Keskin
evidence stability 0.0 -> 1.0 by removing keying/grain noise (the runs recorded identical biology under
inconsistent patient ids + target grain). See agent_core.canonicalize_patient_ids / canonicalize_evidence
/ derive_pool_evidence.
"""
import agent_core as ac


# ---- Fix 1: patient id canonicalization ----

def test_canon_patient_id_numeric_aliases():
    for raw in ("7", "Patient 7", "patient 7", "P7", "P-7", "pt_7", "P.7", "#7", "Pt7"):
        assert ac.canon_patient_id(raw) == "Pt7", raw

def test_canon_patient_id_idempotent_and_leaves_nonnumeric():
    assert ac.canon_patient_id("Pt7") == "Pt7"
    assert ac.canon_patient_id("MEL-21B") == "MEL-21B"     # not a plain numeric alias -> untouched
    assert ac.canon_patient_id("healthy_donor_A") == "healthy_donor_A"
    assert ac.canon_patient_id(None) is None

def test_canonicalize_patient_ids_rewrites_defs_and_refs():
    rec = {
        "patients": [{"paper_local_id": "7"}, {"paper_local_id": "8"}],
        "evidence": [{"patient_paper_id": "7"}, {"patient_paper_id": "Patient 8"}],
        "pools": [{"patient_paper_id": "7"}],
        "immunizing_peptides": [{"patient_paper_id": None}],   # None ref is skipped, not crashed
    }
    n = ac.canonicalize_patient_ids(rec)
    assert n == 5                                          # 2 patient defs + 2 evidence refs + 1 pool ref
    assert [p["paper_local_id"] for p in rec["patients"]] == ["Pt7", "Pt8"]
    assert [e["patient_paper_id"] for e in rec["evidence"]] == ["Pt7", "Pt8"]
    assert rec["pools"][0]["patient_paper_id"] == "Pt7"
    assert ac.canonicalize_patient_ids(rec) == 0            # idempotent


# ---- Fix 2: evidence grain canonicalization ----

def _rec_with_targets():
    return {
        "immunizing_peptides": [
            {"paper_local_id": "IMP_A", "gene_symbol": "ARHGAP35", "sequence": "AAAAAAAAAA"},
            {"paper_local_id": "IMP_A2", "gene_symbol": "ARHGAP35", "sequence": "ZZZZZZZZZZ"},  # same gene, higher seq
        ],
        "candidates": [{"paper_local_id": "CAND_A", "gene_symbol": "ARHGAP35", "sequence": "CCCCCCCCCC"}],
        "epitopes": [{"paper_local_id": "EPI_X", "gene_symbol": "GPC1", "sequence": "GGGGGGGGGG"}],
        "evidence": [
            {"patient_paper_id": "Pt7", "target_kind": "candidate", "candidate_paper_id": "CAND_A",
             "assay": "elispot", "outcome": "immunogenic", "t_cell_subset": "cd8"},
            {"patient_paper_id": "Pt7", "target_kind": "epitope", "epitope_paper_id": "EPI_X",
             "assay": "elispot", "outcome": "immunogenic", "t_cell_subset": "cd4"},
        ],
    }

def test_canonicalize_evidence_repoints_candidate_to_deterministic_imp():
    rec = _rec_with_targets()
    ac.canonicalize_evidence(rec)
    cand_row = rec["evidence"][0]
    assert cand_row["target_kind"] == "immunizing_peptide"
    assert cand_row["immunizing_peptide_paper_id"] == "IMP_A"    # min-sequence rep, deterministic
    assert cand_row["candidate_paper_id"] is None
    # epitope-target row is LEFT ALONE (epitope-specific evidence is legitimate)
    assert rec["evidence"][1]["target_kind"] == "epitope"

def test_canonicalize_evidence_dedups_repoint_collisions_but_keeps_cd4_cd8_separate():
    rec = _rec_with_targets()
    # two rows, same patient/gene/assay/subset via different entities -> collapse to one after re-point
    rec["evidence"] = [
        {"patient_paper_id": "Pt7", "target_kind": "candidate", "candidate_paper_id": "CAND_A",
         "assay": "elispot", "outcome": "immunogenic", "t_cell_subset": "cd8"},
        {"patient_paper_id": "Pt7", "target_kind": "immunizing_peptide", "immunizing_peptide_paper_id": "IMP_A2",
         "assay": "elispot", "outcome": "immunogenic", "t_cell_subset": "cd8"},
        {"patient_paper_id": "Pt7", "target_kind": "candidate", "candidate_paper_id": "CAND_A",
         "assay": "elispot", "outcome": "immunogenic", "t_cell_subset": "cd4"},  # different subset -> kept
    ]
    ac.canonicalize_evidence(rec)
    subsets = sorted(e["t_cell_subset"] for e in rec["evidence"])
    assert subsets == ["cd4", "cd8"]                            # one per subset, collision merged
    assert all(e["immunizing_peptide_paper_id"] == "IMP_A" for e in rec["evidence"])

def test_canonicalize_evidence_ignores_non_immunogenic():
    rec = _rec_with_targets()
    rec["evidence"][0]["outcome"] = "non_immunogenic"
    ac.canonicalize_evidence(rec)
    assert rec["evidence"][0]["target_kind"] == "candidate"      # untouched


# ---- Fix 3: pool-evidence roll-up ----

def _rec_for_pool():
    return {
        "immunizing_peptides": [
            {"paper_local_id": "IMP_A", "gene_symbol": "ARHGAP35", "sequence": "AAAAAAAAAA"},
            {"paper_local_id": "IMP_B", "gene_symbol": "GPC1", "sequence": "BBBBBBBBBB"},
        ],
        "pools": [{"paper_local_id": "POOL-Pt7", "patient_paper_id": "Pt7",
                   "member_peptide_ids": ["IMP_A", "IMP_B"]}],
        "evidence": [
            {"patient_paper_id": "Pt7", "target_kind": "immunizing_peptide",
             "immunizing_peptide_paper_id": "IMP_A", "assay": "elispot", "outcome": "immunogenic"},
        ],
    }

def test_derive_pool_evidence_rolls_up_member_response():
    rec = _rec_for_pool()
    added = ac.derive_pool_evidence(rec)
    assert added == 1
    pool_rows = [e for e in rec["evidence"] if e["target_kind"] == "pool"]
    assert len(pool_rows) == 1
    row = pool_rows[0]
    assert row["patient_paper_id"] == "Pt7" and row["pool_paper_id"] == "POOL-Pt7"
    assert row["outcome"] == "immunogenic" and row["needs_review"] is True
    assert row["magnitude"]["raw"]                              # has a raw so the magnitude nudge stays quiet
    assert ac.derive_pool_evidence(rec) == 0                    # idempotent

def test_derive_pool_evidence_no_response_no_row():
    rec = _rec_for_pool()
    rec["evidence"] = []                                        # nobody responded -> nothing entailed
    assert ac.derive_pool_evidence(rec) == 0

def test_derive_pool_evidence_skips_when_response_not_in_this_pool():
    rec = _rec_for_pool()
    rec["pools"][0]["member_peptide_ids"] = ["IMP_B"]           # responder IMP_A is NOT a member
    assert ac.derive_pool_evidence(rec) == 0


# ---- candidates lane: drop administered re-lists (gold records zero) ----

def test_drop_candidate_that_duplicates_a_peptide_and_repoints_evidence():
    rec = {
        "immunizing_peptides": [{"paper_local_id": "IMP_A", "gene_symbol": "G", "sequence": "AAAAAAAAAA"}],
        "epitopes": [],
        "candidates": [{"paper_local_id": "CAND_A", "gene_symbol": "G", "sequence": "AAAAAAAAAA",
                        "candidate_status": "administered"}],
        "evidence": [{"patient_paper_id": "Pt7", "target_kind": "candidate", "candidate_paper_id": "CAND_A",
                      "assay": "elispot", "outcome": "presented"}],
    }
    dropped = ac.drop_administered_candidates(rec)
    assert dropped == 1 and rec["candidates"] == []
    row = rec["evidence"][0]                                    # re-pointed to the peptide, not dangling
    assert row["target_kind"] == "immunizing_peptide"
    assert row["immunizing_peptide_paper_id"] == "IMP_A" and row["candidate_paper_id"] is None

def test_drop_candidate_that_duplicates_an_epitope():
    rec = {
        "immunizing_peptides": [],
        "epitopes": [{"paper_local_id": "EPI_A", "gene_symbol": "G", "sequence": "SIINFEKL"}],
        "candidates": [{"paper_local_id": "CAND_A", "sequence": "SIINFEKL", "candidate_status": "administered"}],
        "evidence": [{"patient_paper_id": "Pt7", "target_kind": "candidate", "candidate_paper_id": "CAND_A",
                      "assay": "elispot", "outcome": "immunogenic"}],
    }
    assert ac.drop_administered_candidates(rec) == 1
    assert rec["evidence"][0]["target_kind"] == "epitope"
    assert rec["evidence"][0]["epitope_paper_id"] == "EPI_A"

def test_keep_genuinely_unselected_candidate():
    # a candidate whose sequence is in NEITHER lane (a real screened-out prediction) is preserved
    rec = {
        "immunizing_peptides": [{"paper_local_id": "IMP_A", "sequence": "AAAAAAAAAA"}],
        "epitopes": [],
        "candidates": [{"paper_local_id": "CAND_X", "sequence": "WWWWWWWWWW", "candidate_status": "predicted"}],
        "evidence": [],
    }
    assert ac.drop_administered_candidates(rec) == 0
    assert len(rec["candidates"]) == 1

def test_drop_candidates_repoints_tcr_clonotype_specificity():
    # a TcrClonotype points at the candidate via `specific_for_kind` (not target_kind) -> must re-point too
    rec = {
        "immunizing_peptides": [{"paper_local_id": "IMP_A", "sequence": "AAAAAAAAAA"}],
        "epitopes": [],
        "candidates": [{"paper_local_id": "CAND_A", "sequence": "AAAAAAAAAA", "candidate_status": "administered"}],
        "tcr_clonotypes": [{"paper_local_id": "TCR1", "specific_for_kind": "candidate",
                            "candidate_paper_id": "CAND_A"}],
    }
    assert ac.drop_administered_candidates(rec) == 1
    tcr = rec["tcr_clonotypes"][0]
    assert tcr["specific_for_kind"] == "immunizing_peptide"
    assert tcr["immunizing_peptide_paper_id"] == "IMP_A" and tcr["candidate_paper_id"] is None


# ---- Measurement identity: what makes two evidence rows "the same fact" ----
#
# A reactivity matrix measures the SAME patient/epitope/timepoint/stimulation several ways --
# peptide-pulsed, by minigene, against autologous tumour, and again with the response blocked by
# antibody. Those are distinct measurements. Before assay_antigen_format and assay_detail joined
# the key they collapsed onto one row (16 lost on Hu 33479501 'Suppl Dataset 4b' alone).

def _row(**over):
    base = {"patient_paper_id": "Pt3", "target_kind": "epitope", "epitope_paper_id": "ASP_3-ASP7",
            "assay": "elispot", "outcome": "immunogenic", "t_cell_subset": "cd4",
            "assay_stimulation": "in_vitro_expanded", "timepoint_phase": "post_vaccine"}
    base.update(over)
    return base


def test_antigen_format_and_detail_are_part_of_measurement_identity():
    pulsed = _row(assay_antigen_format="peptide_pulsed")
    minigene = _row(assay_antigen_format="minigene")
    blocked = _row(assay_antigen_format="minigene",
                   assay_detail="minigene response blocked by anti-DR antibody")
    keys = {ac._evidence_measurement_key(e) for e in (pulsed, minigene, blocked)}
    assert len(keys) == 3


def test_identical_measurements_still_share_one_key():
    a = _row(assay_antigen_format="peptide_pulsed", magnitude={"raw": "120 SFC"})
    b = _row(assay_antigen_format="peptide_pulsed", magnitude={"raw": "different"})
    # magnitude is NOT identity -- the same measurement reported twice stays one fact
    assert ac._evidence_measurement_key(a) == ac._evidence_measurement_key(b)


def test_canonicalize_evidence_keeps_distinct_assay_formats():
    rec = {"evidence": [_row(assay_antigen_format="peptide_pulsed"),
                        _row(assay_antigen_format="minigene"),
                        _row(assay_antigen_format="autologous_tumor"),
                        _row(assay_antigen_format="peptide_pulsed")]}   # a true duplicate
    ac.canonicalize_evidence(rec)
    assert len(rec["evidence"]) == 3
    assert {e["assay_antigen_format"] for e in rec["evidence"]} == {
        "peptide_pulsed", "minigene", "autologous_tumor"}
