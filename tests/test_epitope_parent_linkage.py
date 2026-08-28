"""Regression: orphaned epitopes from the DETERMINISTIC manifest loader (tcr_phase5, 2026-08-16).

THE DEFECT. `add_epitope_manifest` -> `epitope_manifest_adapter` emits every manifest epitope with NO
`parent_peptide_ids`, and `load_epitope_manifest_into` then LOCKS the epitopes lane. Result: a 100%
orphaned epitope lane that `outer_guard` rejects and that the agent cannot repair -- `clear_entities`
is refused on a locked lane and `add_entities` only appends. Three of the four phase-5 live papers died
here (33479501 203/203 orphaned, 37165196 225 of the 682, 30568305 99/99); the one that never called
the loader finalized clean, and the whole pre-loader phase-4 generation has zero orphans.

THE REPAIR (agent_core.link_orphan_epitopes / quarantine_orphan_epitopes). A minimal epitope is a
SUBSTRING of the long peptide that was immunized, so the parent is DERIVABLE from the record -- the
same containment rule `build_crossreactivity_evidence` already uses. Nothing is invented: a parent must
be a loaded immunizing peptide that actually CONTAINS the epitope; anything left over is quarantined
out of the record, never given a made-up parent.
"""
import json
import pathlib
import tempfile

import openpyxl
import pytest

import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]
REF = json.loads((PKT / "reference_records" / "rojas_extracted.json").read_text())

META = json.dumps({"pmid": "30568305", "title": "t", "journal": "j", "year": 2020,
                   "cohort_size": 8, "indication_summary": "melanoma"})

# Keskin-shape manifest (two-row merged header): the epitope and the long peptide it came from sit on
# the SAME row, so the linkage the loader drops is recoverable by containment.
MANIFEST = [
    ["Patient ID", "Gene", "Protein change", "HLA allele", "Mutated peptide", None,
     "Immunizing peptide", None],
    [None, None, None, None, "Sequence", "Affinity (nM)", "Sequence", "ID"],
    [1, "ATP10B", "p.R821Q", "A01:01", "QTQKHLDLY", 29.87, "VPDINMEKKLRKIRAQTQKHLDLYARDG", "IMP03"],
    [1, "DOK7", "p.Q151K", "A02:01", "IPPAVTGKW", 40.0, "LARDIPPAVTGKWKLSDLRRYGAVPSG", "IMP04"],
]


def _xlsx(rows, sheet="S"):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = sheet
    for r in rows:
        ws.append(list(r))
    p = pathlib.Path(tempfile.mkdtemp()) / "m.xlsx"
    wb.save(p)
    return str(p)


def _new_partial():
    out = str(pathlib.Path(tempfile.mkdtemp()) / "o.json")
    ok, msg = agent_core.init_partial(out, META)
    assert ok, msg
    return out


def _orphans(rec):
    return [e["paper_local_id"] for e in rec["epitopes"] if not e.get("parent_peptide_ids")]


# ---------------------------------------------------------------------------------------------
# the live failure, end to end through the tool the agent actually calls
# ---------------------------------------------------------------------------------------------

def test_manifest_loader_leaves_no_orphaned_epitope():
    """THE regression. Before the fix this lane came back 2/2 orphaned."""
    xlsx = _xlsx(MANIFEST)
    out = _new_partial()
    ok, msg = agent_core.load_epitope_manifest_into(out, xlsx)
    assert ok, msg
    rec = json.loads(pathlib.Path(out + ".partial.json").read_text())
    assert len(rec["epitopes"]) == 2
    assert _orphans(rec) == [], "manifest epitopes came back orphaned"
    # each epitope points at the long peptide from its OWN row
    by_seq = {e["sequence"]: e for e in rec["epitopes"]}
    imp = {p["sequence"]: p["paper_local_id"] for p in rec["immunizing_peptides"]}
    assert by_seq["QTQKHLDLY"]["parent_peptide_ids"] == [imp["VPDINMEKKLRKIRAQTQKHLDLYARDG"]]
    assert by_seq["IPPAVTGKW"]["parent_peptide_ids"] == [imp["LARDIPPAVTGKWKLSDLRRYGAVPSG"]]


def test_loader_locks_the_lane_it_just_linked():
    """Why the linkage MUST happen inside the loader: the lane is locked immediately after, so an
    orphan minted here can never be cleared and re-added by the agent."""
    xlsx = _xlsx(MANIFEST)
    out = _new_partial()
    assert agent_core.load_epitope_manifest_into(out, xlsx)[0]
    ok, msg = agent_core.clear_section(out, "epitopes")
    assert not ok and "locked" in msg          # unfixable-by-agent, hence must arrive correct
    rec = json.loads(pathlib.Path(out + ".partial.json").read_text())
    assert _orphans(rec) == []


def test_loader_reports_linkage_in_its_message():
    out = _new_partial()
    ok, msg = agent_core.load_epitope_manifest_into(out, _xlsx(MANIFEST))
    assert ok and "parent linkage: derived for 2 epitope(s)" in msg


# ---------------------------------------------------------------------------------------------
# the derivation rule itself
# ---------------------------------------------------------------------------------------------

def _rec(epitopes, peptides):
    return {"epitopes": epitopes, "immunizing_peptides": peptides}


def test_containment_is_required_no_parent_is_invented():
    rec = _rec([{"paper_local_id": "E1", "sequence": "WLCSSFMGL", "parent_peptide_ids": []}],
               [{"paper_local_id": "P1", "sequence": "WLSSSFMGLSQNLLLRSPGFRQL"}])   # 1 aa apart
    assert agent_core.link_orphan_epitopes(rec) == 0
    assert rec["epitopes"][0]["parent_peptide_ids"] == []


def test_already_linked_epitopes_are_never_touched():
    rec = _rec([{"paper_local_id": "E1", "sequence": "QTQKHLDLY", "parent_peptide_ids": ["HAND"]}],
               [{"paper_local_id": "P1", "sequence": "IRAQTQKHLDLYARDG"}])
    assert agent_core.link_orphan_epitopes(rec) == 0
    assert rec["epitopes"][0]["parent_peptide_ids"] == ["HAND"]


def test_gene_symbol_disambiguates_a_shared_substring():
    rec = _rec([{"paper_local_id": "E1", "sequence": "AAAKKKLLL", "gene_symbol": "GENEB",
                 "parent_peptide_ids": []}],
               [{"paper_local_id": "PA", "sequence": "QQAAAKKKLLLQQ", "gene_symbol": "GENEA"},
                {"paper_local_id": "PB", "sequence": "WWAAAKKKLLLWW", "gene_symbol": "GENEB"}])
    assert agent_core.link_orphan_epitopes(rec) == 1
    assert rec["epitopes"][0]["parent_peptide_ids"] == ["PB"]


def test_tiling_epitope_keeps_every_containing_parent():
    """schema.MinimalEpitope.parent_peptide_ids is many-to-many by design: one minimal epitope tiling
    two overlapping long peptides is one epitope with two parents."""
    rec = _rec([{"paper_local_id": "E1", "sequence": "AAAKKKLLL", "parent_peptide_ids": []}],
               [{"paper_local_id": "P2", "sequence": "QQAAAKKKLLLQQ"},
                {"paper_local_id": "P1", "sequence": "WWAAAKKKLLLWW"}])
    assert agent_core.link_orphan_epitopes(rec) == 1
    assert rec["epitopes"][0]["parent_peptide_ids"] == ["P1", "P2"]      # sorted -> deterministic


def test_too_many_candidate_parents_stays_orphaned():
    """Past the schema's 12-parent cap the source has not pinned the link down -- leave it orphaned
    rather than truncate an arbitrary twelve."""
    peps = [{"paper_local_id": f"P{i:02d}", "sequence": f"{'AC'*i}AAAKKKLLLQQ"} for i in range(1, 15)]
    rec = _rec([{"paper_local_id": "E1", "sequence": "AAAKKKLLL", "parent_peptide_ids": []}], peps)
    assert agent_core.link_orphan_epitopes(rec) == 0
    assert rec["epitopes"][0]["parent_peptide_ids"] == []


def test_linking_is_idempotent():
    rec = _rec([{"paper_local_id": "E1", "sequence": "QTQKHLDLY", "parent_peptide_ids": []}],
               [{"paper_local_id": "P1", "sequence": "IRAQTQKHLDLYARDG"}])
    assert agent_core.link_orphan_epitopes(rec) == 1
    assert agent_core.link_orphan_epitopes(rec) == 0
    assert rec["epitopes"][0]["parent_peptide_ids"] == ["P1"]


# ---------------------------------------------------------------------------------------------
# quarantine over fabrication
# ---------------------------------------------------------------------------------------------

def _unlinkable_record():
    """The real Hu 33479501 residue: the paper's own 'Immunizing peptide' cell (…LPNEYAFVT) is one
    residue short of containing the epitope it is printed next to, so NO containment parent exists."""
    rec = json.loads(json.dumps(REF))
    parent = rec["immunizing_peptides"][0]["paper_local_id"]
    rec["epitopes"].append({
        "paper_local_id": "EPI_ORPHAN", "sequence": "LPNEYAFVTT", "gene_symbol": "COL22A1",
        "mhc_class": "I", "parent_peptide_ids": [],
        "quoted_text": "COL22A1 LPNEYAFVTT predicted MHC-I epitope",
        "section_ref": "Supplementary Dataset 4a",
    })
    assert parent      # sanity: the reference record does carry peptides
    return rec


def test_unlinkable_epitope_is_quarantined_not_fabricated(tmp_path):
    rec = _unlinkable_record()
    n_eps = len(rec["epitopes"])
    removed, err = agent_core.quarantine_orphan_epitopes(rec)
    assert err is None
    assert [e["paper_local_id"] for e in removed["epitopes"]] == ["EPI_ORPHAN"]
    assert len(rec["epitopes"]) == n_eps - 1
    assert _orphans(rec) == []


def test_orphaned_epitope_msg_spares_evidence_targeted_orphans():
    paper = type("P", (), {})()
    paper.epitopes = [type("E", (), {"paper_local_id": "EPI_ORPHAN", "parent_peptide_ids": []})()]
    paper.evidence = [type("V", (), {"epitope_paper_id": "EPI_ORPHAN"})()]
    assert agent_core._orphaned_epitope_msg(paper) is None
    paper.evidence = [type("V", (), {"epitope_paper_id": None})()]
    assert "EPI_ORPHAN" in agent_core._orphaned_epitope_msg(paper)


def test_evidence_targeted_orphan_is_kept():
    """Hu Supp11b: assayed 15-mers are not nested in the immunizing peptides. Empty parents
    are legal; quarantining them deletes the ELISPOT calls. An epitope evidence already
    targets is not dump noise."""
    rec = _unlinkable_record()
    ev = json.loads(json.dumps(rec["evidence"][0]))
    ev.update({"target_kind": "epitope", "epitope_paper_id": "EPI_ORPHAN",
               "immunizing_peptide_paper_id": None, "pool_paper_id": None, "candidate_paper_id": None,
               "evidence_local_id": None})
    rec["evidence"].append(ev)
    n_ev = len(rec["evidence"])
    n_eps = len(rec["epitopes"])
    removed, err = agent_core.quarantine_orphan_epitopes(rec)
    assert err is None
    assert "epitopes" not in removed
    assert len(rec["epitopes"]) == n_eps
    assert len(rec["evidence"]) == n_ev
    assert any(e.get("paper_local_id") == "EPI_ORPHAN" for e in rec["epitopes"])


def test_a_systematic_orphan_lane_is_refused_not_swept_away():
    """Every epitope unlinkable means the PEPTIDE lane is wrong (wrong sheet/column), not that the
    source is noisy -- refuse and say so instead of quarantining the whole lane."""
    rec = json.loads(json.dumps(REF))
    for e in rec["epitopes"]:
        e["parent_peptide_ids"] = []
    rec["immunizing_peptides"] = []                       # nothing to link against
    n_eps = len(rec["epitopes"])
    assert agent_core.link_orphan_epitopes(rec) == 0
    removed, err = agent_core.quarantine_orphan_epitopes(rec)
    assert removed == {} and err is not None
    assert "too large to be source noise" in err
    assert len(rec["epitopes"]) == n_eps                  # record untouched


def test_finalize_repairs_the_lane_and_writes_the_quarantine_sidecar(tmp_path):
    out = str(tmp_path / "o.json")
    rec = _unlinkable_record()
    for e in rec["epitopes"]:                             # simulate the loader's orphaned arrival
        e["parent_peptide_ids"] = []
    pathlib.Path(out + ".partial.json").write_text(json.dumps(rec))
    ok, msg = agent_core.finalize_partial(
        out, allow_missing_magnitudes=True, allow_missing_pools=True,
        allow_member_level_pool_evidence=True, allow_candidate_bridge_mismatch=True,
        allow_unknown_funnel_size=True, allow_regimen_divergence=True,
        allow_evidence_count_mismatch=True, allow_peptide_count_mismatch=True,
        allow_sparse_evidence=True, allow_missing_class_ii=True,
        allow_missing_minimal_epitopes=True, allow_ungrounded_safety_grade=True,
        allow_missing_tcr_gateway=True, allow_tcr_status_mismatch=True)
    assert ok, msg
    assert "derived parent_peptide_ids" in msg and "QUARANTINED 1" in msg
    saved = json.loads(pathlib.Path(out).read_text())
    assert _orphans(saved) == []
    assert not any(e["paper_local_id"] == "EPI_ORPHAN" for e in saved["epitopes"])
    side = json.loads(pathlib.Path(out + ".quarantine.json").read_text())
    assert [e["paper_local_id"] for e in side["removed"]["epitopes"]] == ["EPI_ORPHAN"]
    assert "not contained in any loaded immunizing peptide" in side["reason"]


def test_finalize_rejects_a_wholly_unlinkable_lane(tmp_path):
    out = str(tmp_path / "o.json")
    rec = json.loads(json.dumps(REF))
    for e in rec["epitopes"]:
        e["parent_peptide_ids"] = []
    rec["immunizing_peptides"] = []
    pathlib.Path(out + ".partial.json").write_text(json.dumps(rec))
    ok, msg = agent_core.finalize_partial(out)
    assert not ok and "too large to be source noise" in msg
    assert pathlib.Path(out + ".partial.json").exists()   # partial kept so the agent can fix it


@pytest.mark.parametrize("stale", [True, False])
def test_init_clears_a_stale_quarantine_sidecar(tmp_path, stale):
    out = str(tmp_path / "o.json")
    if stale:
        pathlib.Path(out + ".quarantine.json").write_text("{}")
    assert agent_core.init_partial(out, META)[0]
    assert not pathlib.Path(out + ".quarantine.json").exists()
