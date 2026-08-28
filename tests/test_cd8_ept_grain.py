"""CD8 / predicted class-I (EPT) evidence must target the 8-11mer, not the long IMP.

Keskin Table S5: HAKLETALR is EPT13A nested in IMP13/IMP14. A quote that names
predicted class I epitopes (EPT) must point at the 9-mer, not the 25-mer.
"""
import json
import pathlib
from types import SimpleNamespace as NS

import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]
META = {"pmid": "12345678", "journal": "Test J", "year": 2024, "title": "CD8 EPT grain test",
        "cohort_size": 1, "indication_summary": "glioblastoma"}

IMP_SEQ = "SEMEENLANRSHAKLETALRDSSRVL"
MER9 = "HAKLETALR"
AFF = {"unit": "unknown", "raw": "n/a"}


def _patient():
    return {"quoted_text": "pt P1", "section_ref": "M", "paper_local_id": "P1",
            "indication": "glioblastoma", "n_peptides_synthesized": 1, "n_peptides_immunogenic": 1}


def _imp():
    return {"paper_local_id": "IMP_GPC1", "sequence": IMP_SEQ, "is_neoantigen": True,
            "gene_symbol": "GPC1", "quoted_text": "immunizing peptide GPC1", "section_ref": "Table S5"}


def _class_i_epitope():
    return {"paper_local_id": "EPI_HAKLETALR", "sequence": MER9, "is_neoantigen": True,
            "mhc_class": "I", "hla_allele": "HLA-A*68:01", "gene_symbol": "GPC1",
            "parent_peptide_ids": ["IMP_GPC1"],
            "predicted_affinity": AFF,
            "quoted_text": "predicted MHC class I neoepitope HAKLETALR", "section_ref": "Table S5"}


def _ev_cd8_ept():
    return {"patient_paper_id": "P1", "target_kind": "epitope",
            "epitope_paper_id": "EPI_HAKLETALR", "assay": "elispot", "outcome": "immunogenic",
            "quoted_text": "mutant-specific CD8+ T cell responses to predicted class I epitopes (EPT)",
            "section_ref": "Results", "t_cell_subset": "cd8", "mhc_class": "class_i",
            "magnitude": {"raw": "positive"}}


def _ev_cd8_imp_unnamed():
    # SLX4MUT shape: says EPT in the generic sense, does not name EPT9A / a 9-mer sequence
    return {"patient_paper_id": "P1", "target_kind": "immunizing_peptide",
            "immunizing_peptide_paper_id": "IMP_GPC1", "assay": "elispot", "outcome": "immunogenic",
            "quoted_text": "mutant-specific CD8+ T cell responses to predicted class I epitopes (EPT) that arose from mutated SLX4 (SLX4MUT)",
            "section_ref": "Fig. 2b", "t_cell_subset": "cd8",
            "assay_detail": "The paper does not state which of the four predicted SLX4 class I EPTs was recognized; recorded on the parent immunizing peptide",
            "magnitude": {"raw": "positive"}}


def _ev_cd8_imp_named_ept():
    return {"patient_paper_id": "P1", "target_kind": "immunizing_peptide",
            "immunizing_peptide_paper_id": "IMP_GPC1", "assay": "elispot", "outcome": "immunogenic",
            "quoted_text": "the predominant CD8+ clone was specific to predicted class I epitope EPT12A, sequence HAKLETALR",
            "section_ref": "Results", "t_cell_subset": "cd8", "mhc_class": "class_i",
            "magnitude": {"raw": "positive"}}


def _ev_cd8_imp_ept_id_only():
    # Names EPT12A but not the 9-mer sequence — H1 cannot uniquely re-point.
    return {"patient_paper_id": "P1", "target_kind": "immunizing_peptide",
            "immunizing_peptide_paper_id": "IMP_GPC1", "assay": "elispot", "outcome": "immunogenic",
            "quoted_text": "the predominant CD8+ clone was specific to predicted class I epitope EPT12A",
            "section_ref": "Results", "t_cell_subset": "cd8", "mhc_class": "class_i",
            "magnitude": {"raw": "positive"}}


def _ev_cd4_imp():
    return {"patient_paper_id": "P1", "target_kind": "immunizing_peptide",
            "immunizing_peptide_paper_id": "IMP_GPC1", "assay": "elispot", "outcome": "immunogenic",
            "quoted_text": "CD4+ T cell responses against mutated GPC1 neoepitope",
            "section_ref": "Figure 2b", "t_cell_subset": "cd4",
            "magnitude": {"raw": "positive"}}


def _assemble(out, *, epitopes, evidence):
    ok, msg = agent_core.init_partial(str(out), json.dumps(META))
    assert ok, msg
    ok, msg = agent_core.append_section(str(out), "patients", json.dumps([_patient()]))
    assert ok, msg
    ok, msg = agent_core.append_section(str(out), "immunizing_peptides", json.dumps([_imp()]))
    assert ok, msg
    if epitopes:
        ok, msg = agent_core.append_section(str(out), "epitopes", json.dumps(epitopes))
        assert ok, msg
    ok, msg = agent_core.append_section(str(out), "evidence", json.dumps(evidence))
    assert ok, msg


def test_gap_empty_for_cd8_on_ept():
    rec = NS(
        immunizing_peptides=[NS(paper_local_id="IMP_GPC1", sequence=IMP_SEQ)],
        epitopes=[NS(paper_local_id="EPI_HAKLETALR", sequence=MER9, mhc_class="I",
                     parent_peptide_ids=["IMP_GPC1"])],
        evidence=[NS(t_cell_subset="cd8", target_kind="epitope",
                     epitope_paper_id="EPI_HAKLETALR", immunizing_peptide_paper_id=None)],
    )
    assert agent_core._cd8_long_peptide_ept_gap(rec) == []


def test_gap_empty_for_cd8_on_imp_when_ept_not_named():
    rec = NS(
        immunizing_peptides=[NS(paper_local_id="IMP_GPC1", sequence=IMP_SEQ)],
        epitopes=[NS(paper_local_id="EPI_HAKLETALR", sequence=MER9, mhc_class="I",
                     parent_peptide_ids=["IMP_GPC1"])],
        evidence=[NS(t_cell_subset="cd8", target_kind="immunizing_peptide",
                     immunizing_peptide_paper_id="IMP_GPC1",
                     quoted_text="CD8+ responses to predicted class I epitopes (EPT) from SLX4MUT",
                     assay_detail="four predicted SLX4 class I EPTs were not deconvoluted")],
    )
    assert agent_core._cd8_long_peptide_ept_gap(rec) == []


def test_gap_empty_when_quote_has_long_imp_sequence():
    # The 9-mer is a substring of the 26-mer; quoting the IMP is not naming the bite.
    rec = NS(
        immunizing_peptides=[NS(paper_local_id="IMP_GPC1", sequence=IMP_SEQ)],
        epitopes=[NS(paper_local_id="EPI_HAKLETALR", sequence=MER9, mhc_class="I",
                     parent_peptide_ids=["IMP_GPC1"])],
        evidence=[NS(t_cell_subset="cd8", target_kind="immunizing_peptide",
                     immunizing_peptide_paper_id="IMP_GPC1",
                     quoted_text=f"CD8 response to immunizing peptide {IMP_SEQ}")],
    )
    assert agent_core._cd8_long_peptide_ept_gap(rec) == []


def test_gap_fires_for_cd8_on_imp_when_ept_id_named():
    rec = NS(
        immunizing_peptides=[NS(paper_local_id="IMP_GPC1", sequence=IMP_SEQ)],
        epitopes=[NS(paper_local_id="EPI_HAKLETALR", sequence=MER9, mhc_class="I",
                     parent_peptide_ids=["IMP_GPC1"],
                     quoted_text="predicted class I neoepitope HAKLETALR EPT12A")],
        evidence=[NS(t_cell_subset="cd8", target_kind="immunizing_peptide",
                     immunizing_peptide_paper_id="IMP_GPC1",
                     quoted_text="the clone was specific to EPT12A")],
    )
    msgs = agent_core._cd8_long_peptide_ept_gap(rec)
    assert msgs and "IMP_GPC1" in msgs[0]


def test_gap_empty_when_quote_names_three_epts():
    """COL22A1: several named EPTs, not deconvoluted → one CD8 row on the IMP, no gap."""
    rec = NS(
        immunizing_peptides=[NS(paper_local_id="IMP_COL", sequence="FPVVQSTEDVFPQGLPNEYAFVT")],
        epitopes=[
            NS(paper_local_id="EPI_FPQGLPNEY", sequence="FPQGLPNEY", mhc_class="I",
               quoted_text="COL22A1 FPQGLPNEY 5-EPT11A"),
            NS(paper_local_id="EPI_LPNEYAFVT", sequence="LPNEYAFVT", mhc_class="I",
               quoted_text="COL22A1 LPNEYAFVT 5-EPT11D"),
        ],
        evidence=[NS(t_cell_subset="cd8", target_kind="immunizing_peptide",
                     immunizing_peptide_paper_id="IMP_COL",
                     quoted_text="CD8+ T-cell response against COL22A1",
                     assay_detail="3 reactive EPT: 5-EPT11A FPQGLPNEY, 5-EPT11C LPNEYAFVTT, 5-EPT11D LPNEYAFVT")],
    )
    assert agent_core._cd8_long_peptide_ept_gap(rec) == []


def test_gap_empty_when_three_ept_ids_only_one_minted():
    """Stick-out EPTs quarantined: three named ids, one 9-mer still in the catalog.

    Without the ≥2 EPT-id guard, H1 would hang the whole CD8 finding on that one 9-mer.
    """
    rec = NS(
        immunizing_peptides=[NS(paper_local_id="IMP_COL", sequence="FPVVQSTEDVFPQGLPNEYAFVT")],
        epitopes=[
            NS(paper_local_id="EPI_FPQGLPNEY", sequence="FPQGLPNEY", mhc_class="I",
               quoted_text="COL22A1 FPQGLPNEY 5-EPT11A"),
        ],
        evidence=[NS(t_cell_subset="cd8", target_kind="immunizing_peptide",
                     immunizing_peptide_paper_id="IMP_COL",
                     quoted_text="CD8+ T-cell response against COL22A1",
                     assay_detail="3 reactive EPT: 5-EPT11A FPQGLPNEY, 5-EPT11C LPNEYAFVTT, 5-EPT11D LPNEYAFVT")],
    )
    assert agent_core._cd8_long_peptide_ept_gap(rec) == []
    d = {
        "immunizing_peptides": [{"paper_local_id": "IMP_COL",
                                 "sequence": "FPVVQSTEDVFPQGLPNEYAFVT"}],
        "epitopes": [{"paper_local_id": "EPI_FPQGLPNEY", "sequence": "FPQGLPNEY",
                      "mhc_class": "I", "quoted_text": "COL22A1 FPQGLPNEY 5-EPT11A"}],
        "evidence": [{"t_cell_subset": "cd8", "target_kind": "immunizing_peptide",
                      "immunizing_peptide_paper_id": "IMP_COL",
                      "quoted_text": "CD8+ T-cell response against COL22A1",
                      "assay_detail": ("3 reactive EPT: 5-EPT11A FPQGLPNEY, "
                                       "5-EPT11C LPNEYAFVTT, 5-EPT11D LPNEYAFVT")}],
    }
    assert agent_core._repoint_cd8_named_ept(d) == 0
    assert d["evidence"][0]["target_kind"] == "immunizing_peptide"


def test_gap_fires_for_cd8_on_imp_when_9mer_sequence_named():
    rec = NS(
        immunizing_peptides=[NS(paper_local_id="IMP_GPC1", sequence=IMP_SEQ)],
        epitopes=[NS(paper_local_id="EPI_HAKLETALR", sequence=MER9, mhc_class="I",
                     parent_peptide_ids=["IMP_GPC1"])],
        evidence=[NS(t_cell_subset="cd8", target_kind="immunizing_peptide",
                     immunizing_peptide_paper_id="IMP_GPC1",
                     quoted_text="CD8 cells recognized HAKLETALR")],
    )
    msgs = agent_core._cd8_long_peptide_ept_gap(rec)
    assert msgs and "EPI_HAKLETALR" in msgs[0]


def test_gap_fires_when_assay_detail_names_9mer():
    rec = NS(
        immunizing_peptides=[NS(paper_local_id="IMP_GPC1", sequence=IMP_SEQ)],
        epitopes=[NS(paper_local_id="EPI_HAKLETALR", sequence=MER9, mhc_class="I",
                     parent_peptide_ids=["IMP_GPC1"])],
        evidence=[NS(t_cell_subset="cd8", target_kind="immunizing_peptide",
                     immunizing_peptide_paper_id="IMP_GPC1",
                     quoted_text="CD8+ response to SLX4MUT",
                     assay_detail="the clone recognized HAKLETALR")],
    )
    msgs = agent_core._cd8_long_peptide_ept_gap(rec)
    assert msgs and "EPI_HAKLETALR" in msgs[0]


def test_gap_empty_for_cd8_on_imp_without_nested_ept():
    rec = NS(
        immunizing_peptides=[NS(paper_local_id="IMP_GPC1", sequence=IMP_SEQ)],
        epitopes=[],
        evidence=[NS(t_cell_subset="cd8", target_kind="immunizing_peptide",
                     immunizing_peptide_paper_id="IMP_GPC1")],
    )
    assert agent_core._cd8_long_peptide_ept_gap(rec) == []


def test_gap_empty_for_cd4_on_long_imp():
    rec = NS(
        immunizing_peptides=[NS(paper_local_id="IMP_GPC1", sequence=IMP_SEQ)],
        epitopes=[NS(paper_local_id="EPI_HAKLETALR", sequence=MER9, mhc_class="I",
                     parent_peptide_ids=["IMP_GPC1"])],
        evidence=[NS(t_cell_subset="cd4", target_kind="immunizing_peptide",
                     immunizing_peptide_paper_id="IMP_GPC1")],
    )
    assert agent_core._cd8_long_peptide_ept_gap(rec) == []


def test_finalize_cd8_on_ept_ok(tmp_path):
    out = tmp_path / "r.json"
    _assemble(out, epitopes=[_class_i_epitope()], evidence=[_ev_cd8_ept()])
    ok, msg = agent_core.finalize_partial(str(out))
    assert ok, msg
    rec = json.loads(out.read_text())
    assert rec["evidence"][0]["target_kind"] == "epitope"
    assert "allow_cd8_on_immunizing_peptide" not in (rec.get("finalize_overrides_used") or [])


def test_finalize_cd8_on_imp_unnamed_ok(tmp_path):
    """Decision A: SLX4-style CD8 on IMP, no named EPT, finalizes without override."""
    out = tmp_path / "r.json"
    _assemble(out, epitopes=[_class_i_epitope()], evidence=[_ev_cd8_imp_unnamed()])
    ok, msg = agent_core.finalize_partial(str(out))
    assert ok, msg
    rec = json.loads(out.read_text())
    assert rec["evidence"][0]["target_kind"] == "immunizing_peptide"
    assert "allow_cd8_on_immunizing_peptide" not in (rec.get("finalize_overrides_used") or [])


def test_repoint_cd8_on_candidate_named_sequence():
    rec = {
        "immunizing_peptides": [{"paper_local_id": "IMP_GPC1", "sequence": IMP_SEQ}],
        "epitopes": [{"paper_local_id": "EPI_HAKLETALR", "sequence": MER9, "mhc_class": "I"}],
        "evidence": [{"t_cell_subset": "cd8", "target_kind": "candidate",
                      "candidate_paper_id": "C1", "quoted_text": "CD8 clone recognized HAKLETALR",
                      "assay_detail": None, "section_ref": None}],
    }
    assert agent_core._repoint_cd8_named_ept(rec) == 1
    assert rec["evidence"][0]["target_kind"] == "epitope"
    assert rec["evidence"][0]["epitope_paper_id"] == "EPI_HAKLETALR"
    assert rec["evidence"][0].get("candidate_paper_id") in (None, "")


def test_finalize_repoints_cd8_named_sequence(tmp_path):
    """H1: quote names the nested 9-mer → finalize re-points to the epitope, no override."""
    out = tmp_path / "r.json"
    _assemble(out, epitopes=[_class_i_epitope()], evidence=[_ev_cd8_imp_named_ept()])
    ok, msg = agent_core.finalize_partial(str(out))
    assert ok, msg
    rec = json.loads(out.read_text())
    assert rec["evidence"][0]["target_kind"] == "epitope"
    assert rec["evidence"][0]["epitope_paper_id"] == "EPI_HAKLETALR"
    assert rec["evidence"][0].get("immunizing_peptide_paper_id") in (None, "")
    assert "allow_cd8_on_immunizing_peptide" not in (rec.get("finalize_overrides_used") or [])


def test_finalize_repoints_cd8_ept_id_via_epitope_quote(tmp_path):
    """Ott shape: evidence quote names EPT12A; epitope quoted_text carries that id."""
    epi = _class_i_epitope()
    epi["quoted_text"] = "predicted MHC class I neoepitope HAKLETALR EPT12A"
    out = tmp_path / "r.json"
    _assemble(out, epitopes=[epi], evidence=[_ev_cd8_imp_ept_id_only()])
    ok, msg = agent_core.finalize_partial(str(out))
    assert ok, msg
    rec = json.loads(out.read_text())
    assert rec["evidence"][0]["target_kind"] == "epitope"
    assert rec["evidence"][0]["epitope_paper_id"] == "EPI_HAKLETALR"


def test_finalize_unresolvable_ept_id_stays_on_imp(tmp_path):
    """EPT12A in the quote but no epitope carries that id → cannot re-point; stay on IMP, no gap."""
    out = tmp_path / "r.json"
    _assemble(out, epitopes=[_class_i_epitope()], evidence=[_ev_cd8_imp_ept_id_only()])
    ok, msg = agent_core.finalize_partial(str(out))
    assert ok, msg
    rec = json.loads(out.read_text())
    assert rec["evidence"][0]["target_kind"] == "immunizing_peptide"
    assert "allow_cd8_on_immunizing_peptide" not in (rec.get("finalize_overrides_used") or [])


def test_finalize_cd4_on_imp_still_ok_with_nested_ept(tmp_path):
    out = tmp_path / "r.json"
    _assemble(out, epitopes=[_class_i_epitope()], evidence=[_ev_cd4_imp()])
    ok, msg = agent_core.finalize_partial(str(out))
    assert ok, msg


def test_tool_schema_exposes_cd8_grain_override():
    from vaxtract.tool_registry import finalize_override_names, get_tool
    spec = get_tool("finalize")
    props = spec.input_schema["properties"]
    names = finalize_override_names()
    assert "allow_cd8_on_immunizing_peptide" in props
    assert "allow_cd8_on_immunizing_peptide" in names
    assert "allow_class_ii_asp_dump" in props
    assert "allow_class_ii_asp_dump" in names


def test_override_is_hard():
    assert "allow_cd8_on_immunizing_peptide" in agent_core.HARD_OVERRIDES
    assert not agent_core.overrides_are_soft_only(["allow_cd8_on_immunizing_peptide"])


def test_finalize_clears_stale_cd8_override_stamp(tmp_path):
    """Replay path: a seed/partial that still lists a HARD CD8 stamp must drop it
    when the gap no longer fires (unnamed SLX4-style CD8 on IMP is legal)."""
    out = tmp_path / "r.json"
    _assemble(out, epitopes=[_class_i_epitope()], evidence=[_ev_cd8_imp_unnamed()])
    partial = pathlib.Path(str(out) + ".partial.json")
    rec = json.loads(partial.read_text())
    rec["finalize_overrides_used"] = ["allow_cd8_on_immunizing_peptide"]
    partial.write_text(json.dumps(rec))
    ok, msg = agent_core.finalize_partial(str(out))
    assert ok, msg
    saved = json.loads(out.read_text())
    assert saved.get("finalize_overrides_used") == []
    assert saved["evidence"][0]["target_kind"] == "immunizing_peptide"


def test_refinalize_refuses_to_overwrite_source(tmp_path):
    p = tmp_path / "r.json"
    p.write_text("{}")
    ok, msg = agent_core.refinalize(str(p), str(p))
    assert not ok
    assert "overwrite" in msg


def test_refinalize_idempotent_on_clean_cd8_on_ept(tmp_path):
    src = tmp_path / "src.json"
    _assemble(src, epitopes=[_class_i_epitope()], evidence=[_ev_cd8_ept()])
    ok, msg = agent_core.finalize_partial(str(src))
    assert ok, msg
    before = json.loads(src.read_text())
    dest = tmp_path / "dest.json"
    ok, msg = agent_core.refinalize(str(src), str(dest))
    assert ok, msg
    after = json.loads(dest.read_text())
    assert after.get("finalize_overrides_used") == []
    assert len(after["evidence"]) == len(before["evidence"])
    assert after["evidence"][0]["target_kind"] == "epitope"
    assert after["evidence"][0]["epitope_paper_id"] == "EPI_HAKLETALR"
    assert src.read_text()  # source keep-file untouched as a path; dest is the replay


def test_refinalize_repoints_unique_named_and_clears_stamp(tmp_path):
    """Ott leftover shape: unique named CD8 still on IMP + unnamed CD8 + stale stamp.
    Replay re-points only the unique named rows and drops the HARD stamp."""
    src = tmp_path / "src.json"
    named = _ev_cd8_imp_named_ept()
    unnamed = _ev_cd8_imp_unnamed()
    _assemble(src, epitopes=[_class_i_epitope()], evidence=[named, unnamed])
    ok, msg = agent_core.finalize_partial(str(src))
    assert ok, msg
    # Simulate a keep-file saved before H1: unique named row hung on the IMP + fossil stamp.
    rec = json.loads(src.read_text())
    rec["evidence"][0]["target_kind"] = "immunizing_peptide"
    rec["evidence"][0]["immunizing_peptide_paper_id"] = "IMP_GPC1"
    rec["evidence"][0]["epitope_paper_id"] = None
    rec["finalize_overrides_used"] = ["allow_cd8_on_immunizing_peptide"]
    src.write_text(json.dumps(rec))
    dest = tmp_path / "dest.json"
    ok, msg = agent_core.refinalize(str(src), str(dest))
    assert ok, msg
    out = json.loads(dest.read_text())
    kinds = [(e.get("target_kind"), e.get("epitope_paper_id"), e.get("immunizing_peptide_paper_id"))
             for e in out["evidence"]]
    assert ("epitope", "EPI_HAKLETALR", None) in kinds
    assert any(k[0] == "immunizing_peptide" and k[2] == "IMP_GPC1" for k in kinds)
    assert out.get("finalize_overrides_used") == []
    assert json.loads(src.read_text())["finalize_overrides_used"] == [
        "allow_cd8_on_immunizing_peptide"
    ]
    assert not pathlib.Path(str(dest) + ".partial.json").exists()
