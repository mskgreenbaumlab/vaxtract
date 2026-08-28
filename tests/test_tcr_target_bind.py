"""Keskin clonotype target bind: unique entity / unique flag, never mixed-sheet guesses."""
from vaxtract.schema import ExtractedPaper, TcrClonotype
from vaxtract.tcr_target_bind import apply_reporter_refs, bind_clonotype, bind_record

PROV = {"quoted_text": "q", "section_ref": "S"}


def _imp(pid, seq, gene, patient="Pt8"):
    return {
        "paper_local_id": pid, "sequence": seq, "gene_symbol": gene,
        "is_neoantigen": True, "patient_paper_id": patient, **PROV,
    }


def _flag(patient, kind, tid, gene_in_quote, section_ref="Supplementary Table 8a"):
    return {
        "target_kind": kind,
        "immunizing_peptide_paper_id": tid if kind == "immunizing_peptide" else None,
        "epitope_paper_id": tid if kind == "epitope" else None,
        "candidate_paper_id": None,
        "patient_paper_id": patient,
        "tcr_identified": "cd4",
        "quoted_text": gene_in_quote,
        "section_ref": section_ref,
        "needs_review": False,
    }


def _clono(pid="TCR_Pt8_mutantSVEP1_1", patient="Pt8", antigen="mutant SVEP1",
           section_ref="Supplementary Table 8b"):
    return {
        "paper_local_id": pid,
        "patient_paper_id": patient,
        "receptor_modality": "single_cell_tcr_seq",
        "chain_pairing": "paired_ab",
        "chains": [
            {"locus": "TRA", "junction_aa": "CAVNDYKLSF", "cdr3_scheme": "imgt_junction_aa"},
            {"locus": "TRB", "junction_aa": "CASSLGQPQHF", "cdr3_scheme": "imgt_junction_aa"},
        ],
        "specificity_evidence": "not_established",
        "needs_review": False,
        "quoted_text": f"{patient} {antigen}-reactive single-cell clonotype (CDR3a=CAVNDYKLSF, CDR3b=CASSLGQPQHF; 2 cells)",
        "section_ref": section_ref,
    }


def _rec(**kwargs):
    rec = {
        "pmid": "30568305",
        "title": "t",
        "journal": "j",
        "year": 2018,
        "cohort_size": 2,
        "indication_summary": "melanoma",
        "patients": [{"paper_local_id": "Pt8", "indication": "melanoma",
                      "n_peptides_synthesized": 3, "n_peptides_immunogenic": 3, **PROV},
                     {"paper_local_id": "Pt7", "indication": "melanoma",
                      "n_peptides_synthesized": 2, "n_peptides_immunogenic": 2, **PROV}],
        "immunizing_peptides": [
            _imp("IMP_GPPAHVENAIAR_25_440275", "GPPAHVENAIARGIHYQYGDMITYS", "SVEP1"),
            _imp("IMP_PQVDGEIPLHRS_25_878e6c", "PQVDGEIPLHRSDRVKVLSIGEGGF", "SHANK2"),
            _imp("IMP_KPYQPQVDGEIP_27_b78284", "KPYQPQVDGEIPLHRSDRVKVLSIGEG", "SHANK2"),
            _imp("IMP_HNLDLAEKDFMV_20_8eee2c", "HNLDLAEKDFMVNTVAGAMK", "ARHGAP35", "Pt7"),
            _imp("IMP_ITFTTAATHREK_16_eadbf0", "ITFTTAATHREKLQGR", "SLX4", "Pt7"),
        ],
        "epitopes": [{
            "paper_local_id": "EPI_MVNTVAGAMK", "sequence": "MVNTVAGAMK",
            "gene_symbol": "ARHGAP35", "is_neoantigen": True, "mhc_class": "I",
            "predicted_affinity": {"unit": "unknown", "raw": "n", "tier": "predicted"},
            **PROV,
        }],
        "evidence": [{
            "patient_paper_id": "Pt7", "target_kind": "immunizing_peptide",
            "immunizing_peptide_paper_id": "IMP_HNLDLAEKDFMV_20_8eee2c",
            "assay": "tcr_reporter", "outcome": "positive",
            "evidence_local_id": "EV_TCR_H02",
            "quoted_text": "H02 and F10 specific for ARHGAP35MUT",
            "section_ref": "Fig. 4c",
        }, {
            "patient_paper_id": "Pt7", "target_kind": "epitope",
            "epitope_paper_id": "EPI_MVNTVAGAMK",
            "assay": "tcr_reporter", "outcome": "positive",
            "evidence_local_id": "EV_TCR_CD8_EPT12A",
            "quoted_text": "predominant CD8 clone specific to EPT12A",
            "section_ref": "Extended Data Fig. 7d",
        }],
        "neoantigen_tcr_flags": [
            _flag("Pt8", "immunizing_peptide", "IMP_PQVDGEIPLHRS_25_878e6c",
                  "stimulation with SHANK2MUT and SVEP1MUT", "Supplementary Table 8a"),
            _flag("Pt8", "immunizing_peptide", "IMP_GPPAHVENAIAR_25_440275",
                  "stimulation with SHANK2MUT and SVEP1MUT", "Supplementary Table 8b"),
        ],
        "tcr_clonotypes": [],
        "pools": [],
        "candidates": [],
    }
    rec.update(kwargs)
    return rec


def test_unique_gene_binds_ids_and_stays_not_established():
    rec = _rec()
    cl = _clono()
    why = bind_clonotype(rec, cl, antigen="mutant SVEP1",
                         section_ref="Supplementary Table 8b",
                         patient_paper_id="Pt8")
    assert why.startswith("bound")
    assert cl["specific_for_kind"] == "immunizing_peptide"
    assert cl["immunizing_peptide_paper_id"] == "IMP_GPPAHVENAIAR_25_440275"
    assert cl["epitope_paper_id"] is None
    assert cl["needs_review"] is False
    assert cl["specificity_evidence"] == "not_established"
    TcrClonotype.model_validate(cl)


def test_two_imps_same_gene_without_unique_flag_unbound():
    rec = _rec(neoantigen_tcr_flags=[])
    cl = _clono(pid="TCR_Pt8_mutantSHANK2_1", antigen="mutant SHANK2",
                section_ref="Supplementary Table 8a")
    why = bind_clonotype(rec, cl, antigen="mutant SHANK2",
                         section_ref="Supplementary Table 8a",
                         patient_paper_id="Pt8")
    assert why.startswith("skip")
    assert cl.get("specific_for_kind") is None
    assert cl.get("immunizing_peptide_paper_id") is None
    assert cl["specificity_evidence"] == "not_established"


def test_shank2_unique_flag_binds_flag_target():
    rec = _rec()
    cl = _clono(pid="TCR_Pt8_mutantSHANK2_1", antigen="mutant SHANK2",
                section_ref="Supplementary Table 8a")
    why = bind_clonotype(rec, cl, antigen="mutant SHANK2",
                         section_ref="Supplementary Table 8a",
                         patient_paper_id="Pt8")
    assert "IMP_PQVDGEIPLHRS_25_878e6c" in why
    assert cl["immunizing_peptide_paper_id"] == "IMP_PQVDGEIPLHRS_25_878e6c"
    assert cl["specificity_evidence"] == "not_established"
    assert cl["needs_review"] is False


def test_mixed_antigen_string_unbound():
    rec = _rec()
    cl = _clono(pid="TCR_Pt7_mixed_1", patient="Pt7",
                antigen="ARHGAP35 / SLX4 mutant EPT (CD8+)",
                section_ref="Supplementary Table 11c")
    why = bind_clonotype(
        rec, cl, antigen="ARHGAP35 / SLX4 mutant EPT (CD8+)",
        section_ref="Supplementary Table 11c", patient_paper_id="Pt7",
    )
    assert "mixed" in why
    assert cl.get("specific_for_kind") is None


def test_h02_quoted_text_sets_reporter_ref():
    rec = _rec()
    cl = _clono(pid="TCR_H02", patient="Pt7", antigen="ARHGAP35",
                section_ref="Fig. 4c")
    cl["quoted_text"] = "CD4+ TCR H02 expressed in Jurkat, specific for ARHGAP35MUT"
    apply_reporter_refs(rec, clonotypes=[cl])
    assert cl["validation_evidence_ref"] == "EV_TCR_H02"
    assert cl["specificity_evidence"] == "functional_tcr_clone"
    assert cl["immunizing_peptide_paper_id"] == "IMP_HNLDLAEKDFMV_20_8eee2c"
    TcrClonotype.model_validate(cl)


def test_reporter_skips_when_evidence_local_id_missing():
    rec = _rec()
    rec["evidence"] = [e for e in rec["evidence"] if e.get("evidence_local_id") != "EV_TCR_H02"]
    cl = _clono(pid="TCR_H02", patient="Pt7", antigen="ARHGAP35",
                section_ref="Fig. 4c")
    cl["quoted_text"] = "clone H02"
    report = apply_reporter_refs(rec, clonotypes=[cl])
    assert cl.get("validation_evidence_ref") is None
    assert report["n_linked"] == 0
    assert any("EV_TCR_H02" in (s.get("reason") or "") for s in report["skipped"])


def test_bind_record_schema_valid_and_mixed_counted():
    rec = _rec()
    rec["tcr_clonotypes"] = [
        _clono(),
        _clono(pid="TCR_Pt8_mutantSHANK2_1", antigen="mutant SHANK2",
               section_ref="Supplementary Table 8a"),
        _clono(pid="TCR_Pt7_mixed_1", patient="Pt7",
               antigen="ARHGAP35 / SLX4 mutant EPT (CD8+)",
               section_ref="Supplementary Table 11c"),
    ]
    report = bind_record(rec)
    assert report["n_bound"] == 2
    assert report["n_mixed"] == 1
    ExtractedPaper.model_validate(rec)
