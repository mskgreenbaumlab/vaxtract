"""#5: pr_compare.qc_metrics — the per-paper batch QC (source-fact verification + optional gold P/R)."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "pr_compare"))
import pr_compare  # noqa: E402


def test_normalize_collapses_and_casefolds():
    assert pr_compare.normalize("  Hi\tThere… ") == "hi there..."


def test_qc_no_source_is_failsoft():
    rec = {"immunizing_peptides": [{"paper_local_id": "i1", "sequence": "SIINFEKL"}]}
    m = pr_compare.qc_metrics(rec)  # no paper_dir -> no source index
    assert m["fact_seq_rate"] is None
    assert m["evidence_precision"] is None and m["evidence_recall"] is None


def test_qc_fact_rate_catches_unsourced_sequence(tmp_path):
    (tmp_path / "src.txt").write_text("the vaccine peptide SIINFEKL was tested by ELISpot")
    rec = {"immunizing_peptides": [
        {"paper_local_id": "i1", "sequence": "SIINFEKL"},        # in source
        {"paper_local_id": "i2", "sequence": "MADEUPSEQ"},       # NOT in source -> unverified
    ]}
    m = pr_compare.qc_metrics(rec, paper_dir=str(tmp_path))
    assert m["fact_seq_rate"] == 0.5
    assert m["fact_unverified"] == 1


def test_qc_gold_pr_added_when_gold_given(tmp_path):
    (tmp_path / "src.txt").write_text("SIINFEKL")
    rec = {"immunizing_peptides": [{"paper_local_id": "i1", "sequence": "SIINFEKL"}],
           "evidence": [{"patient_paper_id": "P1", "immunizing_peptide_paper_id": "i1",
                         "assay": "elispot", "outcome": "immunogenic"}]}
    gold = {"immunizing_peptides": [{"paper_local_id": "X9", "sequence": "SIINFEKL"}],
            "evidence": [{"patient_paper_id": "P1", "immunizing_peptide_paper_id": "X9",
                          "assay": "elispot", "outcome": "immunogenic"}]}
    m = pr_compare.qc_metrics(rec, paper_dir=str(tmp_path), gold=gold)
    # the evidence rows bridge to the SAME sequence under different local ids -> they MATCH (P/R = 1.0)
    assert m["evidence_precision"] == 1.0 and m["evidence_recall"] == 1.0
