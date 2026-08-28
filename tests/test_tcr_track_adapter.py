"""Deterministic tet-spec track loader (vaxtract.tcr_track_adapter): unpivot a wide
CDR3b × timepoint sheet into TcrClonotype + nested observations, reproducibly. Second half is hermetic
coverage of the wired tool path (agent_core.load_tcr_track_into / MCP `add_tcr_track`): merge/idempotency,
the tcr_seq_status sticker, per-patient dedupe, and the not-a-track-sheet refusal."""
import hashlib
import json
import pathlib
import tempfile

import openpyxl

import agent_core as ac
from vaxtract import tcr_track_adapter as ad

# header row: col0=CDR3b, col1=tet-spec flag, col2..=timepoint labels
HEADER = ["CDR3β", "1=yes", "Pre-Vax", "Week 4", "Long-term"]
ROWS = [
    ["CASSLGNQPQHF", 1, 0.01, 0.05, None],       # tet-spec: 2 timepoints (Long-term blank -> skipped)
    ["AGQGFNQPQHF", None, 0.02, 0.03, 0.04],      # NOT tet-spec (flag blank) -> dropped entirely
    ["CASSDISYEQYF", 1, None, None, 0.007],       # tet-spec: 1 timepoint
]


def _parse():
    return ad.parse_tetraspec_track(HEADER, ROWS, patient_paper_id="2")


def test_keeps_only_tetramer_specific_rows():
    clones = _parse()
    assert len(clones) == 2                       # the flag-blank row is dropped
    seqs = {c.chains[0].junction_aa for c in clones}
    assert "CASSLGNQPQHF" in seqs and "AGQGFNQPQHF" not in seqs

def test_unpivots_timepoints_into_observations():
    clones = {c.chains[0].junction_aa: c for c in _parse()}
    assert len(clones["CASSLGNQPQHF"].observations) == 2   # Pre-Vax + Week 4 (Long-term blank skipped)
    assert len(clones["CASSDISYEQYF"].observations) == 1   # Long-term only
    obs = clones["CASSLGNQPQHF"].observations
    assert [(o.timepoint_label, o.frequency_value) for o in obs] == [("Pre-Vax", 0.01), ("Week 4", 0.05)]

def test_frequency_stored_as_reported_fraction():
    o = _parse()[0].observations[0]
    assert o.frequency_basis == "fraction" and o.frequency_value == 0.01

def test_timepoint_phase_mapped():
    obs = {o.timepoint_label: o.timepoint_phase for c in _parse() for o in c.observations}
    assert obs["Pre-Vax"] == "pre_vaccine" and obs["Week 4"] == "on_treatment" and obs["Long-term"] == "memory"

def test_identity_is_paper_explicit_no_review():
    c = _parse()[0]
    assert c.observations_identity == "paper_explicit" and c.needs_review is False

def test_junction_form_recognized():
    # C…F sequences are junctions, not cdr3-proper
    ch = _parse()[0].chains[0]
    assert ch.junction_aa == "CASSLGNQPQHF" and ch.cdr3_aa is None and ch.cdr3_scheme == "imgt_junction_aa"

def test_deterministic():
    def digest():
        return hashlib.sha1(json.dumps([c.model_dump() for c in _parse()], sort_keys=True,
                                       default=str).encode()).hexdigest()
    assert digest() == digest()   # pure function of the input -> identical every run

def test_outputs_are_schema_valid_clonotypes():
    # parse_tetraspec_track constructs schema.TcrClonotype instances, so any invalid row would raise.
    import schema
    assert all(isinstance(c, schema.TcrClonotype) for c in _parse())


# --- the shape variants that occur inside the ONE Hu supplement -----------------------------------------

# Pt3-Pt6 write the identity as the paper's comma-joined 'V,J,CDR3' token, on TCRa sheets too, with a
# patient-specific week set and a label row padded with blank trailing columns.
VJ_HEADER = ["TCR⍺", "1=yes", "Pre-Vax", "Week 3", "Week 20", None, None]
VJ_ROWS = [
    ["TRAV2,TRAJ39,CAVEDNNAGNMLTF", "1", None, 0.027359781, 0.003345824, None, None],
    ["TRAV22,TRAJ45,CAGNSGGGADGLTF", None, 0.01, None, None, None, None],   # unflagged -> dropped
]


def test_splits_the_vj_cdr3_identity_token():
    c = ad.parse_tetraspec_track(VJ_HEADER, VJ_ROWS, patient_paper_id="Pt3")[0]
    ch = c.chains[0]
    assert ch.junction_aa == "CAVEDNNAGNMLTF"          # the CDR3 is the LAST comma field
    assert (ch.v_call, ch.j_call) == ("TRAV2", "TRAJ39")
    assert ch.raw == "TRAV2,TRAJ39,CAVEDNNAGNMLTF"     # verbatim token kept (lossless)


def test_locus_inferred_from_the_identity_header():
    a = ad.parse_tetraspec_track(VJ_HEADER, VJ_ROWS, patient_paper_id="Pt3")[0]
    assert a.chains[0].locus == "TRA" and a.chain_pairing == "tra_only"
    b = _parse()[0]
    assert b.chains[0].locus == "TRB" and b.chain_pairing == "trb_only"


def test_ambiguous_header_raises_rather_than_guessing_the_locus():
    import pytest
    with pytest.raises(ValueError, match="cannot tell which chain"):
        ad.parse_tetraspec_track(["clone", "1=yes", "Pre-Vax"], [["CASSLGNQPQHF", 1, 0.01]],
                                 patient_paper_id="Pt3")
    # ...unless the caller states it
    assert ad.parse_tetraspec_track(["clone", "1=yes", "Pre-Vax"], [["CASSLGNQPQHF", 1, 0.01]],
                                    patient_paper_id="Pt3", locus="TRB")[0].chains[0].locus == "TRB"


def test_any_week_label_is_on_treatment():
    # the weeks DIFFER by patient (3/8/20 here vs 4/12/16/24 elsewhere) -- a rule, not an enumeration
    obs = {o.timepoint_label: o.timepoint_phase
           for c in ad.parse_tetraspec_track(VJ_HEADER, VJ_ROWS, patient_paper_id="Pt3")
           for o in c.observations}
    assert obs["Week 3"] == "on_treatment" and obs["Week 20"] == "on_treatment"


def test_blank_trailing_label_columns_are_ignored():
    c = ad.parse_tetraspec_track(VJ_HEADER, VJ_ROWS, patient_paper_id="Pt3")[0]
    assert [o.timepoint_label for o in c.observations] == ["Week 3", "Week 20"]


def test_row_whose_identity_is_not_a_cdr3_is_skipped():
    rows = [["not a sequence!!", 1, 0.01, None, None, None, None], *VJ_ROWS]
    assert len(ad.parse_tetraspec_track(VJ_HEADER, rows, patient_paper_id="Pt3")) == 1


def test_same_cdr3_on_both_loci_gets_distinct_ids():
    trb = _parse()[0]
    tra = ad.parse_tetraspec_track(["TCR⍺", "1=yes", "Pre-Vax"], [["CASSLGNQPQHF", 1, 0.01]],
                                   patient_paper_id="2")[0]
    assert trb.paper_local_id != tra.paper_local_id and trb.paper_local_id.endswith("_TRB_CASSLGNQPQHF")


# --- wired tool path (agent_core.load_tcr_track_into) ---------------------------------------------------

META = json.dumps({"pmid": "33479501", "title": "t", "journal": "j", "year": 2021,
                   "cohort_size": 8, "indication_summary": "melanoma"})


def _track_xlsx(rows=ROWS, title="Suppl 8 Pt2 Tet-spec TCRb track"):
    """The Hu-shaped sheet: a title row, then the timepoint-label header (row index 1), then data."""
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = title
    ws.append(["Pt2 tetramer-specific TCRb tracking"])          # title row (row 0)
    ws.append(HEADER)                                           # label row (row 1 = header_row_index)
    for r in rows:
        ws.append(r)
    p = pathlib.Path(tempfile.mkdtemp()) / "track.xlsx"
    wb.save(p)
    return str(p), ws.title


def _fresh_record():
    out = str(pathlib.Path(tempfile.mkdtemp()) / "o.json")
    ok, _ = ac.init_partial(out, META)
    assert ok
    return out


def _partial(out):
    return json.loads(pathlib.Path(out + ".partial.json").read_text())


def test_tool_path_loads_trajectories_and_is_idempotent():
    path, sheet = _track_xlsx()
    out = _fresh_record()
    ok, msg = ac.load_tcr_track_into(out, path, sheet, "Pt2")
    assert ok and "loaded 2" in msg and "3 timepoint observations" in msg   # 2 + 1 observations
    rec = _partial(out)
    assert len(rec["tcr_clonotypes"]) == 2
    assert sum(len(c["observations"]) for c in rec["tcr_clonotypes"]) == 3
    ok2, msg2 = ac.load_tcr_track_into(out, path, sheet, "Pt2")
    assert ok2 and "loaded 0" in msg2                           # (patient, chain) dedupe -> nothing new
    assert len(_partial(out)["tcr_clonotypes"]) == 2


def test_tool_path_locks_the_lane():
    path, sheet = _track_xlsx()
    out = _fresh_record()
    assert ac.load_tcr_track_into(out, path, sheet, "Pt2")[0]
    ok, msg = ac.clear_section(out, "tcr_clonotypes")           # sticky: a loader owns this lane
    assert not ok and "tcr_clonotypes" in msg
    assert len(_partial(out)["tcr_clonotypes"]) == 2


def test_status_sticker_bulk_and_upgrade_to_both():
    path, sheet = _track_xlsx()
    out = _fresh_record()
    ac.load_tcr_track_into(out, path, sheet, "Pt2")
    assert _partial(out)["tcr_seq_status"] == "bulk"            # a track sheet IS bulk repertoire

    out2 = _fresh_record()
    rec = _partial(out2)
    rec["tcr_seq_status"] = "single_cell"                       # single-cell clonotypes already loaded
    pathlib.Path(out2 + ".partial.json").write_text(json.dumps(rec))
    ac.load_tcr_track_into(out2, path, sheet, "Pt2")
    assert _partial(out2)["tcr_seq_status"] == "both"           # never downgraded


def test_same_cdr3_in_two_patients_is_not_deduped():
    path, sheet = _track_xlsx()
    out = _fresh_record()
    ac.load_tcr_track_into(out, path, sheet, "Pt2")
    ok, msg = ac.load_tcr_track_into(out, path, sheet, "Pt3")   # a public CDR3b recurs across patients
    assert ok and "loaded 2" in msg
    rec = _partial(out)
    assert len(rec["tcr_clonotypes"]) == 4
    assert {c["patient_paper_id"] for c in rec["tcr_clonotypes"]} == {"Pt2", "Pt3"}


def test_refuses_a_sheet_with_no_tetramer_specific_rows():
    unflagged = [[r[0], None, *r[2:]] for r in ROWS]            # same sheet, flag column blank
    path, sheet = _track_xlsx(rows=unflagged)
    out = _fresh_record()
    ok, msg = ac.load_tcr_track_into(out, path, sheet, "Pt2")
    assert not ok and "no tetramer-specific tracked clonotypes" in msg
    assert _partial(out)["tcr_clonotypes"] == []                # nothing added, record unchanged


def test_section_ref_override_reaches_the_rows():
    path, sheet = _track_xlsx()
    out = _fresh_record()
    ac.load_tcr_track_into(out, path, sheet, "Pt2", section_ref="Supplementary Dataset 8 (Pt2)")
    assert all(c["section_ref"] == "Supplementary Dataset 8 (Pt2)"
               for c in _partial(out)["tcr_clonotypes"])
