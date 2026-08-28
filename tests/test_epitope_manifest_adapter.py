"""H2 stability-lab work wired into the live extractor: the deterministic epitope-manifest loader.

Hermetic — synthesizes the two real manifest SHAPES in-memory (no corpus dependency):
  - Keskin 30568305 "Table S5": two-row MERGED header, class-I only, explicit affinity columns.
  - Rojas 37165196 "targets_with_elispot": single-row header, MHC-I AND MHC-II epitope columns, no
    printed affinity (the class-I affinity slot must still be minted as a lossless unit='unknown').
  - Li 33879241 "Table S2" BOUNDARY: only the long 21-mer is a column, so epitopes must come back EMPTY
    (the loader never invents a predicted minimal binder).

Covers the adapter, its schema-validity, and the agent_core tool path (merge-dedupe + idempotency).
"""
import json
import pathlib
import tempfile

import openpyxl
import pytest

import agent_core
from vaxtract import epitope_manifest_adapter as ema
from cancervac_packet import schema

META = json.dumps({"pmid": "30568305", "title": "t", "journal": "j", "year": 2020,
                   "cohort_size": 8, "indication_summary": "melanoma"})


def _xlsx(rows, sheet="S"):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = sheet
    for r in rows:
        ws.append(list(r))
    p = pathlib.Path(tempfile.mkdtemp()) / "m.xlsx"
    wb.save(p)
    return str(p)


KESKIN = [
    ["Supplementary Table 5. class I prediction related to the immunizing peptides"],
    [],
    ["Patient ID", "Immunizing pool", "Gene", "Protein change", "Peptide length", "HLA allele",
     "Mutated peptide", None, "Wild type peptide", None, "Immunizing peptide", None, "EPT peptide ID"],
    [None, None, None, None, None, None, "Sequence", "Affinity (nM)", "Sequence", "Affinity (nM)",
     "Sequence", "ID", None],
    [1, "A", "ATP10B", "p.R821Q", 9, "A01:01", "QTQKHLDLY", 29.87, "RTQKHLDLY", 69.02,
     "VPDINMEKKLRKIRAQTQKHLDLYARDG", "IMP03", "EPT1"],
    [1, "A", "GENE2", "p.G12D", 9, "B07:02", "SRALEEKKG", 12.5, "SRALEEKKV", 40.0,
     "LAALCPASRALEEKKGNYVVTDHGSCV", "IMP04", "EPT2"],
]

ROJAS = [
    ["Patient number", "Neoantigen number", "Gene", "RefSeq transcript", "Substitution",
     "Mutant Neoantigen Sequence", "WT Neoantigen Sequence", "mRNA",
     "MHC-I  Allele (Best Prediction)", "MHC-I Mutant Epitope (Best Prediction)", "MHC-I WT Epitope",
     "MHC-II  Allele (Best Prediction)", "MHC-II Mutant Epitope (Best Prediction)", "MHC-II WT Epitope"],
    [1, 3, "GPR75", "NM_006794", "D245N", "NAQVRKCPPVITVNASRPQPFMGVPVQ", "NAQVRKCPPVITVDASRPQPFMGVPVQ",
     "ACGT", "HLA-A*32:01", "ITVNASRPQPF", "ITVDASRPQPF", "HLA-DRB1*03:01", "VITVNASRPQPFMGV",
     "VITVDASRPQPFMGV"],
]

LI_S2 = [
    ["Table S2. Selective neoantigens"],
    ["Mutation", "MT 21-mer seq*", "MT score", "WT score", "Fold change", "H-2 allele",
     "Normal VAF", "Tumor VAF", "RNA VAF", "Gene FPKM", "ELISPOT"],
    ["Tmem101.G96V", "QLASTYTAYIVGYVHYGDWLK", 109, 11934, 109.5, "Kb", 0.78, 36.8, 45.78, 7.83, "-"],
]


def test_keskin_shape_two_row_header():
    m = ema.load_epitope_manifest(_xlsx(KESKIN))
    assert {e["sequence"] for e in m["epitopes"]} == {"QTQKHLDLY", "SRALEEKKG"}
    assert {p["sequence"] for p in m["immunizing_peptides"]} == {
        "VPDINMEKKLRKIRAQTQKHLDLYARDG", "LAALCPASRALEEKKGNYVVTDHGSCV"}
    e = next(e for e in m["epitopes"] if e["sequence"] == "QTQKHLDLY")
    assert e["mhc_class"] == "I" and e["hla_allele"] == "HLA-A*01:01"
    assert e["predicted_affinity"]["value"] == 29.87 and e["wild_type_sequence"] == "RTQKHLDLY"
    assert "EPT1" in e["quoted_text"]


def test_skips_epitope_not_contained_in_row_immunizing_peptide():
    rows = [
        ["Gene", "HLA allele", "Mutated peptide Sequence", "Immunizing peptide Sequence", "EPT peptide ID"],
        ["COL22A1", "B35:01", "LPNEYAFVTT", "FPVVQSTEDVFPQGLPNEYAFVT", "5-EPT11C"],
        ["COL22A1", "B35:01", "FPQGLPNEY", "FPVVQSTEDVFPQGLPNEYAFVT", "5-EPT11A"],
    ]
    m = ema.load_epitope_manifest(_xlsx(rows))
    seqs = {e["sequence"] for e in m["epitopes"]}
    assert "FPQGLPNEY" in seqs
    assert "LPNEYAFVTT" not in seqs
    e = next(e for e in m["epitopes"] if e["sequence"] == "FPQGLPNEY")
    assert "5-EPT11A" in e["quoted_text"]


def test_rojas_shape_single_row_dual_class():
    m = ema.load_epitope_manifest(_xlsx(ROJAS))
    seqs = {e["sequence"]: e for e in m["epitopes"]}
    assert set(seqs) == {"ITVNASRPQPF", "VITVNASRPQPFMGV"}          # one MHC-I + one MHC-II per row
    assert seqs["ITVNASRPQPF"]["mhc_class"] == "I"
    assert seqs["VITVNASRPQPFMGV"]["mhc_class"] == "II"
    # class-I with no printed affinity still gets the lossless slot the schema requires
    assert seqs["ITVNASRPQPF"]["predicted_affinity"]["unit"] == "unknown"
    assert {p["sequence"] for p in m["immunizing_peptides"]} == {"NAQVRKCPPVITVNASRPQPFMGVPVQ"}


def test_li_boundary_peptides_only_no_invented_epitopes():
    m = ema.load_epitope_manifest(_xlsx(LI_S2, sheet="Table S2"), sheet="Table S2")
    assert m["epitopes"] == []                                     # minimal epitope is not a column
    assert {p["sequence"] for p in m["immunizing_peptides"]} == {"QLASTYTAYIVGYVHYGDWLK"}


# ---- immunizing-peptide id collisions (live Keskin 30568305 finalize hard-fail, 2026-08-16) ----
#
# Supplementary manifests carry CONTAINMENT VARIANTS: the same neoantigen synthesized at two
# lengths. Under the old `IMP_{seq[:12]}` scheme both rows minted the same paper_local_id and
# `finalize` hard-failed with `duplicate paper_local_id in immunizing_peptides` (6 pairs in Keskin;
# no allow-flag overrides that check, so the paper produced NO record at all).

# a short peptide and its length-extended synthesis variant (real Keskin pair, Pt1 DOK7)
CONTAINMENT = [
    ["Patient ID", "Gene", "Protein change", "HLA allele", "Mutated peptide", None,
     "Immunizing peptide", None],
    [None, None, None, None, "Sequence", "Affinity (nM)", "Sequence", "ID"],
    [1, "DOK7", "p.Q151K", "A01:01", "RDIPPAVTG", 30.0, "LARDIPPAVTGKWKLSDLRRYGA", "IMP01"],
    [1, "DOK7", "p.Q151K", "A02:01", "IPPAVTGKW", 40.0, "LARDIPPAVTGKWKLSDLRRYGAVPSG", "IMP02"],
]

# harder case: same first 12 aa AND the same length — a mutant/WT long-peptide pair whose
# substitution sits past position 12, which a length suffix alone would NOT separate
SAME_LEN_VARIANTS = [
    ["Patient ID", "Gene", "Protein change", "HLA allele", "Mutated peptide", None,
     "Immunizing peptide", None],
    [None, None, None, None, "Sequence", "Affinity (nM)", "Sequence", "ID"],
    [2, "BTD", "p.V39L", "A01:01", "GSYLVALGA", 30.0, "GSYLVALGAHTGEESVADW", "IMP01"],
    [2, "BTD", "p.V39L", "A01:01", "GSYLVALGT", 31.0, "GSYLVALGAHTGEESVADY", "IMP02"],
]


def _imp_ids(rows):
    return [p["paper_local_id"] for p in ema.load_epitope_manifest(_xlsx(rows))["immunizing_peptides"]]


def test_containment_variants_get_distinct_ids():
    ids = _imp_ids(CONTAINMENT)
    assert len(ids) == 2 and len(set(ids)) == 2, f"containment variants collided: {ids}"
    assert all(i.startswith("IMP_LARDIPPAVTGK") for i in ids)     # readable prefix survives


def test_same_prefix_same_length_variants_get_distinct_ids():
    ids = _imp_ids(SAME_LEN_VARIANTS)
    assert len(ids) == 2 and len(set(ids)) == 2, f"same-length variants collided: {ids}"


@pytest.mark.parametrize("rows", [KESKIN, ROJAS, LI_S2, CONTAINMENT, SAME_LEN_VARIANTS])
def test_immunizing_peptide_ids_are_unique_per_manifest(rows):
    ids = _imp_ids(rows)
    assert len(set(ids)) == len(ids)


def test_imp_id_is_a_pure_function_of_the_sequence():
    # stable across runs, across sheets and across papers: same sequence -> same id, and it is
    # NOT a function of row order / position (no counters, no PYTHONHASHSEED dependence)
    seq = "LARDIPPAVTGKWKLSDLRRYGA"
    first = _imp_ids(CONTAINMENT)[0]
    assert first == ema._imp_id(seq) == _imp_ids(CONTAINMENT)[0]
    reordered = [CONTAINMENT[0], CONTAINMENT[1], CONTAINMENT[3], CONTAINMENT[2]]
    assert set(_imp_ids(reordered)) == set(_imp_ids(CONTAINMENT))
    assert _imp_ids(reordered)[1] == first                        # same seq, different row -> same id


def test_containment_variant_record_survives_the_uniqueness_validator():
    # the exact check that hard-failed live: ExtractedPaper's duplicate-paper_local_id guard
    m = ema.load_epitope_manifest(_xlsx(CONTAINMENT))
    rec = json.loads(META)
    rec["immunizing_peptides"] = m["immunizing_peptides"]
    ok, msg = agent_core.validate_record(json.dumps(rec))
    assert ok, msg


@pytest.mark.parametrize("rows,sheet", [(KESKIN, "S"), (ROJAS, "S"), (LI_S2, "Table S2")])
def test_adapter_output_validates_against_schema(rows, sheet):
    m = ema.load_epitope_manifest(_xlsx(rows, sheet=sheet), sheet=sheet)
    for e in m["epitopes"]:
        schema.MinimalEpitope(**e)
    for p in m["immunizing_peptides"]:
        schema.ImmunizingPeptide(**p)


def _partial(xlsx, sheet=None):
    out = str(pathlib.Path(tempfile.mkdtemp()) / "o.json")
    ok, msg = agent_core.init_partial(out, META)
    assert ok, msg
    return out


def test_tool_path_loads_and_is_idempotent():
    xlsx = _xlsx(KESKIN)
    out = _partial(xlsx)
    ok, msg = agent_core.load_epitope_manifest_into(out, xlsx)
    assert ok and "epitopes +2" in msg and "immunizing_peptides +2" in msg
    ok2, msg2 = agent_core.load_epitope_manifest_into(out, xlsx)         # second call adds nothing
    assert ok2 and "epitopes +0" in msg2 and "immunizing_peptides +0" in msg2
    rec = json.loads(pathlib.Path(out + ".partial.json").read_text())
    assert len(rec["epitopes"]) == 2 and len(rec["immunizing_peptides"]) == 2


def test_tool_rejects_non_manifest_sheet():
    # a plain 2-column list (no allele/epitope columns) is not a manifest -> clear failure, no mutation
    xlsx = _xlsx([["Epitope ID", "Amino acid sequence"], ["pp65", "GILARNLVPMVATVQGQNLK"]])
    out = _partial(xlsx)
    ok, msg = agent_core.load_epitope_manifest_into(out, xlsx)
    assert not ok and "no epitope/immunizing-peptide manifest recognized" in msg
    rec = json.loads(pathlib.Path(out + ".partial.json").read_text())
    assert rec["epitopes"] == [] and rec["immunizing_peptides"] == []


# ---- sticky lanes: a loaded lane is locked against clear_entities (the live tcr-wired run2 collapse) ----

def test_loaded_epitope_lane_is_locked_against_clear():
    out = _partial(_xlsx(KESKIN))
    agent_core.load_epitope_manifest_into(out, _xlsx(KESKIN))
    assert agent_core._read_locks(out) >= {"epitopes", "immunizing_peptides"}
    n = len(json.loads(pathlib.Path(out + ".partial.json").read_text())["epitopes"])
    ok, msg = agent_core.clear_section(out, "epitopes")         # the run2 attack
    assert not ok and "locked" in msg.lower()
    assert len(json.loads(pathlib.Path(out + ".partial.json").read_text())["epitopes"]) == n

def test_locked_lane_still_accepts_appends_and_unlocked_lane_clears():
    out = _partial(_xlsx(KESKIN))
    agent_core.load_epitope_manifest_into(out, _xlsx(KESKIN))
    before = len(json.loads(pathlib.Path(out + ".partial.json").read_text())["epitopes"])
    ok, _ = agent_core.append_section(out, "epitopes", json.dumps([{
        "paper_local_id": "EPI_extra", "sequence": "SIINFEKLA", "is_neoantigen": True, "mhc_class": "II",
        "hla_allele": "HLA-DRB1*01:01", "quoted_text": "a class-II epitope the manifest lacked",
        "section_ref": "Fig 3"}]))
    assert ok and len(json.loads(pathlib.Path(out + ".partial.json").read_text())["epitopes"]) == before + 1
    assert agent_core.clear_section(out, "evidence")[0] is True   # an unlocked lane still clears

def test_reinit_refused_when_a_loaded_lane_is_locked():
    # live run3 bypass: blocked from clear_entities, the agent re-init'd to wipe the record + lock.
    out = _partial(_xlsx(KESKIN))
    agent_core.load_epitope_manifest_into(out, _xlsx(KESKIN))
    n = len(json.loads(pathlib.Path(out + ".partial.json").read_text())["epitopes"])
    ok, msg = agent_core.init_partial(out, META)
    assert not ok and "re-init" in msg.lower()
    assert len(json.loads(pathlib.Path(out + ".partial.json").read_text())["epitopes"]) == n  # intact

def test_init_clears_a_stale_lock_with_no_partial():
    # a lock file left by an interrupted run (no partial) is not a live record -> init proceeds + clears it
    out = str(pathlib.Path(tempfile.mkdtemp()) / "o.json")
    agent_core._lock_lanes(out, "epitopes")                      # stale lock, no partial written
    assert agent_core._read_locks(out) == {"epitopes"}
    ok, _ = agent_core.init_partial(out, META)
    assert ok and agent_core._read_locks(out) == set()


# ---- the DOUBLE-WRITE: a loader-owned lane must not hold the same entity twice ----------------
# Hu 33479501 phase-6 run1 holds 125 immunizing peptides TWICE: 125 loader-minted
# `IMP_<seq>_<len>_<sha>` plus 125 the agent hand-transcribed under the paper's own ids
# ('1-IMP04', ...). The lock stopped clear_entities but nothing stopped the second WRITE. Its
# evidence then pointed at the hand-minted ids, which is why an id-keyed comparison of run1 against
# the keep-file reads Jaccard 0.000 while the biology agrees ~80%.

def _handwritten_copy_of_the_manifest(out):
    """What the agent transcribed by hand: the same peptides, under the paper's own ids."""
    rec = json.loads(pathlib.Path(out + ".partial.json").read_text())
    return [{"paper_local_id": f"1-IMP{i:02d}", "sequence": p["sequence"],
             "gene_symbol": p.get("gene_symbol"), "is_neoantigen": True,
             "quoted_text": "transcribed from Supplementary Table 5", "section_ref": "Supp Table 5"}
            for i, p in enumerate(rec["immunizing_peptides"])]


def test_hand_appending_the_loaders_own_peptides_does_not_double_the_lane():
    out = _partial(_xlsx(KESKIN))
    agent_core.load_epitope_manifest_into(out, _xlsx(KESKIN))
    n = len(json.loads(pathlib.Path(out + ".partial.json").read_text())["immunizing_peptides"])
    ok, msg = agent_core.append_section(out, "immunizing_peptides",
                                        json.dumps(_handwritten_copy_of_the_manifest(out)))
    assert ok, msg
    after = json.loads(pathlib.Path(out + ".partial.json").read_text())["immunizing_peptides"]
    assert len(after) == n, f"lane doubled: {len(after)} vs {n}"
    assert "MERGED AWAY" in msg and "1-IMP00" in msg


def test_the_merged_away_id_is_aliased_so_a_reference_to_it_still_resolves():
    """Without the alias the fix would trade a silent double-write for a dangling reference."""
    out = _partial(_xlsx(KESKIN))
    agent_core.load_epitope_manifest_into(out, _xlsx(KESKIN))
    agent_core.append_section(out, "immunizing_peptides",
                              json.dumps(_handwritten_copy_of_the_manifest(out)))
    aliases = agent_core.read_id_aliases(out)
    assert "1-IMP00" in aliases and aliases["1-IMP00"].startswith("IMP_")
    # an evidence row written against the id the agent THOUGHT it created
    rec = {"evidence": [{"target_kind": "immunizing_peptide",
                         "immunizing_peptide_paper_id": "1-IMP00", "patient_paper_id": "Pt1"}],
           "epitopes": [{"parent_peptide_ids": ["1-IMP00", "1-IMP01"]}]}
    n = agent_core.apply_id_aliases(rec, aliases)
    assert n == 3
    assert rec["evidence"][0]["immunizing_peptide_paper_id"] == aliases["1-IMP00"]
    assert rec["epitopes"][0]["parent_peptide_ids"] == [aliases["1-IMP00"], aliases["1-IMP01"]]


def test_a_genuinely_new_peptide_still_appends_to_a_locked_lane():
    out = _partial(_xlsx(KESKIN))
    agent_core.load_epitope_manifest_into(out, _xlsx(KESKIN))
    n = len(json.loads(pathlib.Path(out + ".partial.json").read_text())["immunizing_peptides"])
    ok, msg = agent_core.append_section(out, "immunizing_peptides", json.dumps([{
        "paper_local_id": "IMP_new", "sequence": "MKWVTFISLLLLFSSAYSRGV", "is_neoantigen": True,
        "quoted_text": "a peptide the manifest sheet did not carry", "section_ref": "Fig 1"}]))
    assert ok and "MERGED AWAY" not in msg
    after = json.loads(pathlib.Path(out + ".partial.json").read_text())["immunizing_peptides"]
    assert len(after) == n + 1


def test_same_sequence_different_allele_epitopes_stay_distinct():
    """Epitope identity is (sequence, hla_allele, mhc_class) -- one peptide presented by two alleles
    is two facts, and collapsing them would destroy real data."""
    out = _partial(_xlsx(KESKIN))
    agent_core.load_epitope_manifest_into(out, _xlsx(KESKIN))
    before = len(json.loads(pathlib.Path(out + ".partial.json").read_text())["epitopes"])
    aff = {"unit": "unknown", "raw": "not reported", "method": "NetMHCpan", "tier": "predicted"}
    rows = [{"paper_local_id": f"E{i}", "sequence": "QTQKHLDLY", "is_neoantigen": True,
             "mhc_class": "I", "hla_allele": a, "predicted_affinity": aff,
             "quoted_text": "q", "section_ref": "s"}
            for i, a in enumerate(("HLA-B*07:02", "HLA-C*07:01"))]
    ok, msg = agent_core.append_section(out, "epitopes", json.dumps(rows))
    assert ok and "MERGED AWAY" not in msg
    after = json.loads(pathlib.Path(out + ".partial.json").read_text())["epitopes"]
    assert len(after) == before + 2


def test_an_unlocked_lane_is_never_deduped():
    """PMID 33064988 mints PER-PATIENT peptide ids for one sequence ('IMP-B1-AAA', 'IMP-M2-AAA') and
    its pools reference them patient by patient; the Rojas gold likewise holds 232 peptides over 223
    distinct sequences. Dedup must not reach a lane no deterministic loader owns."""
    out = _partial(_xlsx(KESKIN))                       # NO loader call -> no lock
    rows = [{"paper_local_id": pid, "sequence": "AAAVGVGKSAL", "is_neoantigen": True,
             "quoted_text": "q", "section_ref": "s"} for pid in ("IMP-B1-AAA", "IMP-M2-AAA")]
    ok, msg = agent_core.append_section(out, "immunizing_peptides", json.dumps(rows))
    assert ok and "MERGED AWAY" not in msg
    rec = json.loads(pathlib.Path(out + ".partial.json").read_text())
    assert len(rec["immunizing_peptides"]) == 2


def test_init_record_clears_stale_aliases():
    out = _partial(_xlsx(KESKIN))
    agent_core.load_epitope_manifest_into(out, _xlsx(KESKIN))
    agent_core.append_section(out, "immunizing_peptides",
                              json.dumps(_handwritten_copy_of_the_manifest(out)))
    assert agent_core.read_id_aliases(out)
    pathlib.Path(out + ".partial.json").unlink()        # a genuine restart, not a lock bypass
    pathlib.Path(out + ".locked.json").unlink()
    agent_core.init_partial(out, META)
    assert agent_core.read_id_aliases(out) == {}
