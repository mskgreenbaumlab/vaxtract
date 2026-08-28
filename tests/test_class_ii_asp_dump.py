"""O2: overlapping ASP 15-mers dumped as class-II epitopes (Ott SI5), not Rojas named-DR."""
import json
from types import SimpleNamespace as NS

import agent_core


AFF = {"value": 2.16, "unit": "pct_rank", "raw": "2.16", "method": "IEDB"}
CORE = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKL"  # 30 aa


def _imp(i, seq):
    return NS(paper_local_id=f"IMP{i}", sequence=seq)


def _ii(i, seq, parent, **kw):
    d = dict(paper_local_id=f"EPTII{i}", sequence=seq, mhc_class="II",
             parent_peptide_ids=[parent])
    d.update(kw)
    return NS(**d)


def _dump(n_imp=5, n_off=12, *, hla=False, aff=False):
    imps, epis = [], []
    n = 0
    for i in range(n_imp):
        imps.append(_imp(i, CORE))
        for off in range(n_off):
            kw = {}
            if hla:
                kw["hla_allele"] = "HLA-DPA1*01:03"
            if aff:
                kw["predicted_affinity"] = dict(AFF)
            epis.append(_ii(n, CORE[off:off + 15], f"IMP{i}", **kw))
            n += 1
    return imps, epis


def test_ott_shaped_dump_is_dropped_even_with_predicted_hla():
    """fix3 shape: SI5 stamps DP/DR on every 15-mer but stores no affinity."""
    imps, epis = _dump(hla=True, aff=False)
    rec = {
        "immunizing_peptides": [{"paper_local_id": p.paper_local_id, "sequence": p.sequence}
                                for p in imps],
        "epitopes": [vars(e) for e in epis],
        "evidence": [],
    }
    assert agent_core._drop_class_ii_asp_windows(rec) == 60
    assert rec["epitopes"] == []
    wrapped = NS(immunizing_peptides=imps, epitopes=[], evidence=[])
    assert agent_core._class_ii_asp_dump_gap(wrapped) == []


def test_rojas_named_dr_with_affinity_is_kept():
    iseq = "SEMEENLANRSHAKLETALRDSSRVL"
    rec = {
        "immunizing_peptides": [{"paper_local_id": "IMP0", "sequence": iseq}],
        "epitopes": [{
            "paper_local_id": "EPTII0", "sequence": iseq[:15], "mhc_class": "II",
            "hla_allele": "HLA-DRB1*04:01", "parent_peptide_ids": ["IMP0"],
            "predicted_affinity": dict(AFF),
            "quoted_text": "HLA-DRB1*04:01-restricted class II epitope",
        }],
        "evidence": [],
    }
    assert agent_core._drop_class_ii_asp_windows(rec) == 0
    assert rec["epitopes"][0]["paper_local_id"] == "EPTII0"
    wrapped = NS(
        immunizing_peptides=[_imp(0, iseq)],
        epitopes=[_ii(0, iseq[:15], "IMP0", hla_allele="HLA-DRB1*04:01",
                      predicted_affinity=dict(AFF))],
        evidence=[],
    )
    assert agent_core._class_ii_asp_dump_gap(wrapped) == []


def test_named_ii_with_evidence_not_dropped():
    imps, epis = _dump(hla=False, aff=False)
    rec = {
        "immunizing_peptides": [{"paper_local_id": p.paper_local_id, "sequence": p.sequence}
                                for p in imps],
        "epitopes": [vars(e) for e in epis],
        "evidence": [{"target_kind": "epitope", "epitope_paper_id": f"EPTII{k}"}
                     for k in range(len(epis))],
    }
    assert agent_core._drop_class_ii_asp_windows(rec) == 0
    assert len(rec["epitopes"]) == 60
    wrapped = NS(immunizing_peptides=imps, epitopes=epis,
                 evidence=[NS(target_kind="epitope", epitope_paper_id=f"EPTII{k}")
                           for k in range(60)])
    assert agent_core._class_ii_asp_dump_gap(wrapped) == []


def test_affinity_dump_leftover_still_warns():
    """Windows that have a score are not auto-dropped; a large leftover dump still nudges."""
    imps, epis = _dump(hla=True, aff=True)
    rec = {
        "immunizing_peptides": [{"paper_local_id": p.paper_local_id, "sequence": p.sequence}
                                for p in imps],
        "epitopes": [vars(e) for e in epis],
        "evidence": [],
    }
    assert agent_core._drop_class_ii_asp_windows(rec) == 0
    wrapped = NS(immunizing_peptides=imps, epitopes=epis, evidence=[])
    msgs = agent_core._class_ii_asp_dump_gap(wrapped)
    assert msgs and "survived auto-drop" in msgs[0]


def test_rojas_shaped_one_per_imp_does_not_fire():
    imps, epis = [], []
    core = "SEMEENLANRSHAKLETALRDSSRVL"
    for i in range(8):
        iseq = ("G" * i) + core
        imps.append(_imp(i, iseq))
        epis.append(_ii(i, iseq[i:i + 15], f"IMP{i}", predicted_affinity=dict(AFF),
                        hla_allele="HLA-DRB1*04:01"))
    rec = NS(immunizing_peptides=imps, epitopes=epis, evidence=[])
    assert agent_core._class_ii_asp_dump_gap(rec) == []
    drec = {
        "immunizing_peptides": [{"paper_local_id": p.paper_local_id, "sequence": p.sequence}
                                for p in imps],
        "epitopes": [vars(e) for e in epis],
        "evidence": [],
    }
    assert agent_core._drop_class_ii_asp_windows(drec) == 0


def test_override_is_hard():
    assert "allow_class_ii_asp_dump" in agent_core.HARD_OVERRIDES
    assert not agent_core.overrides_are_soft_only(["allow_class_ii_asp_dump"])


def test_finalize_drops_ott_dump_keeps_named_dr(tmp_path):
    out = tmp_path / "r.json"
    meta = {"pmid": "12345678", "journal": "Test J", "year": 2024, "title": "ASP drop",
            "cohort_size": 1, "indication_summary": "melanoma"}
    patient = {"quoted_text": "pt P1", "section_ref": "M", "paper_local_id": "P1",
               "indication": "melanoma", "n_peptides_synthesized": 5, "n_peptides_immunogenic": 1}
    imps, dump = [], []
    n = 0
    for i in range(5):
        core_i = ("W" * i) + CORE
        imps.append({"paper_local_id": f"IMP{i}", "sequence": core_i, "is_neoantigen": True,
                     "quoted_text": f"IMP{i}", "section_ref": "SI5"})
        for off in range(12):
            dump.append({
                "paper_local_id": f"EPTII{n}", "sequence": core_i[off:off + 15],
                "mhc_class": "II", "hla_allele": "HLA-DPA1*01:03",
                "parent_peptide_ids": [f"IMP{i}"], "is_neoantigen": True,
                "quoted_text": "predicted class II 15-mer", "section_ref": "SI5",
            })
            n += 1
    named = {
        "paper_local_id": "EPTII_NAMED", "sequence": CORE[:15], "mhc_class": "II",
        "hla_allele": "HLA-DRB1*04:01", "parent_peptide_ids": ["IMP0"],
        "predicted_affinity": dict(AFF), "is_neoantigen": True,
        "quoted_text": "HLA-DRB1*04:01-restricted class II epitope", "section_ref": "ED Table 2",
    }
    ev = {"patient_paper_id": "P1", "target_kind": "immunizing_peptide",
          "immunizing_peptide_paper_id": "IMP0", "assay": "elispot", "outcome": "immunogenic",
          "quoted_text": "CD4 response to IMP0", "section_ref": "Fig 1", "t_cell_subset": "cd4",
          "magnitude": {"raw": "positive"}}
    ok, msg = agent_core.init_partial(str(out), json.dumps(meta))
    assert ok, msg
    ok, msg = agent_core.append_section(str(out), "patients", json.dumps([patient]))
    assert ok, msg
    ok, msg = agent_core.append_section(str(out), "immunizing_peptides", json.dumps(imps))
    assert ok, msg
    ok, msg = agent_core.append_section(str(out), "epitopes", json.dumps(dump + [named]))
    assert ok, msg
    ok, msg = agent_core.append_section(str(out), "evidence", json.dumps([ev]))
    assert ok, msg
    ok, msg = agent_core.finalize_partial(str(out))
    assert ok, msg
    rec = json.loads(out.read_text())
    ids = {e["paper_local_id"] for e in rec["epitopes"]}
    assert "EPTII_NAMED" in ids
    assert not any(i.startswith("EPTII") and i != "EPTII_NAMED" for i in ids)
    assert "allow_class_ii_asp_dump" not in (rec.get("finalize_overrides_used") or [])
    assert "ASP-window class-II" in msg
