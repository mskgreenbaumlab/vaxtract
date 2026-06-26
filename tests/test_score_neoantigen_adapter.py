"""Unit tests for the score_extraction neoantigen-identity adapter (2026-06-26).

Regression guard for the bug where agent neoantigens scored 0 true-positives because
ImmunizingPeptide carries no patient field — the patient must be recovered from the
NeoantigenCandidate that selected the peptide, or every (patient, gene, mutation)
match fails on the empty patient slot.
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "eval"))
import score_extraction as sx  # noqa: E402


def test_build_peptide_patient_map_links_candidate_to_peptide():
    agent = {"candidates": [{"selected_peptide_id": "IMP1", "patient_paper_id": "P1"}]}
    assert sx.build_peptide_patient_map(agent) == {"IMP1": "P1"}


def test_patient_recovered_so_agent_matches_gold_grain():
    # Peptide has gene+mutation but NO patient; the candidate supplies the patient.
    agent = {
        "candidates": [{"selected_peptide_id": "IMP1", "patient_paper_id": "P1",
                        "gene_symbol": "ATP10B", "mutation": "p.R821Q"}],
        "immunizing_peptides": [{"paper_local_id": "IMP1", "gene_symbol": "ATP10B",
                                 "mutation": "p.R821Q"}],  # no patient_paper_id
    }
    agent_set = sx.agent_neoantigen_set(agent)
    # gold side keys the same neoantigen as ("Pt1", "ATP10B", "R821Q")
    gold_id = sx.neo_identity("Pt1", "ATP10B", "R821Q")
    assert gold_id == ("1", "ATP10B", "R821Q")
    assert gold_id in agent_set, "patient recovery must make the agent neoantigen match the gold grain"


def test_without_recovery_patient_slot_would_be_empty():
    # Sanity: a peptide alone (no candidate) keeps an empty patient -> does NOT match gold.
    agent = {"immunizing_peptides": [{"paper_local_id": "IMP1", "gene_symbol": "ATP10B",
                                      "mutation": "p.R821Q"}]}
    agent_set = sx.agent_neoantigen_set(agent)
    assert ("", "ATP10B", "R821Q") in agent_set
    assert ("1", "ATP10B", "R821Q") not in agent_set


def test_degenerate_identity_without_gene_dropped():
    agent = {"immunizing_peptides": [{"paper_local_id": "X", "gene_symbol": None, "mutation": None}]}
    assert sx.agent_neoantigen_set(agent) == set()


def test_mutation_and_patient_normalisation():
    assert sx.neo_identity("Patient 3", "kras", "p.G12D") == ("3", "KRAS", "G12D")
