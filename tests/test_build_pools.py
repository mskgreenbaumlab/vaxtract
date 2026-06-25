"""Deterministic per-patient pool builder (`build_patient_pools`).

33064988's per-patient pools were hand-built via add_entities and swung 0<->34 run-to-run
(live 2026-06-09). build_patient_pools groups the per-patient peptide-ASSIGNMENT sheet
('Patient Alias' | 'Peptide Sequence') by patient and matches each sequence to an
already-loaded immunizing_peptide -> one pool per patient, identical every run.
"""
import json
import pathlib

import openpyxl

import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]
REF = json.loads((PKT / "reference_records" / "rojas_extracted.json").read_text())
META = {k: v for k, v in REF.items() if not isinstance(v, list)}


def _partial(out):
    return json.loads(pathlib.Path(str(out) + ".partial.json").read_text())


def _init(out, peptides):
    """A partial record at `out` with the given immunizing_peptides already loaded."""
    pathlib.Path(out).parent.mkdir(parents=True, exist_ok=True)
    ok, msg = agent_core.init_partial(str(out), json.dumps(META))
    assert ok, msg
    ok, msg = agent_core.append_section(str(out), "immunizing_peptides", json.dumps(peptides))
    assert ok, msg


# short tokens -> valid 10-mer sequences (schema requires 8-55 AA over ACDEFGHIKLMNPQRSTVWY)
SEQ = {"AAA": "AAAAAAAAAA", "BBB": "CCCCCCCCCC", "CCC": "DDDDDDDDDD", "ZZZ": "EEEEEEEEEE"}


def _pep(pid, tok):
    return {"paper_local_id": pid, "sequence": SEQ[tok], "is_neoantigen": True,
            "quoted_text": f"{pid} synthesized", "section_ref": "Vaccine peptides"}


def _assignment_sheet(path, rows, header=("Patient Alias", "Diagnosis", "Peptide Sequence")):
    wb = openpyxl.Workbook(); ws = wb.active; ws.title = "Vaccine peptides"
    ws.append(list(header))
    for r in rows:
        ws.append([SEQ.get(c, c) for c in r])   # expand peptide tokens to real sequences
    wb.save(path)


def test_builds_one_pool_per_patient_matched_by_sequence(tmp_path):
    # ids embed the patient token (IMP-<patient>-<seq>) — the 33064988 scheme
    peps = [_pep("IMP-B1-AAA", "AAA"), _pep("IMP-B1-BBB", "BBB"),
            _pep("IMP-M2-AAA", "AAA"), _pep("IMP-M2-CCC", "CCC")]
    out = tmp_path / "rec.json"
    _init(out, peps)
    xlsx = tmp_path / "mmc2.xlsx"
    _assignment_sheet(xlsx, [("B1", "Bladder", "AAA"), ("B1", "Bladder", "BBB"),
                             ("M2", "Melanoma", "AAA"), ("M2", "Melanoma", "CCC")])
    ok, msg = agent_core.build_patient_pools(
        str(out), str(xlsx), "Patient Alias", "Peptide Sequence", section_ref="Supp Table 8")
    assert ok, msg
    pools = {p["patient_paper_id"]: p for p in _partial(out)["pools"]}
    assert set(pools) == {"B1", "M2"}
    # patient-token disambiguation: shared sequence 'AAA' resolves to the right patient's id
    assert set(pools["B1"]["member_peptide_ids"]) == {"IMP-B1-AAA", "IMP-B1-BBB"}
    assert set(pools["M2"]["member_peptide_ids"]) == {"IMP-M2-AAA", "IMP-M2-CCC"}
    assert pools["B1"]["paper_local_id"] == "POOL-B1"
    assert pools["B1"]["section_ref"] == "Supp Table 8"


def test_deterministic_across_two_calls(tmp_path):
    peps = [_pep("IMP-B1-AAA", "AAA"), _pep("IMP-B1-BBB", "BBB")]
    rows = [("B1", "Bladder", "AAA"), ("B1", "Bladder", "BBB")]
    seen = []
    for i in range(2):
        out = tmp_path / f"r{i}" / "rec.json"
        _init(out, peps)
        xlsx = tmp_path / f"r{i}" / "m.xlsx"
        _assignment_sheet(xlsx, rows)
        ok, _ = agent_core.build_patient_pools(str(out), str(xlsx), "Patient Alias", "Peptide Sequence")
        assert ok
        p = _partial(out)["pools"]
        seen.append([(x["paper_local_id"], sorted(x["member_peptide_ids"])) for x in p])
    assert seen[0] == seen[1] == [("POOL-B1", ["IMP-B1-AAA", "IMP-B1-BBB"])]


def test_column_by_index(tmp_path):
    peps = [_pep("IMP-B1-AAA", "AAA")]
    out = tmp_path / "rec.json"
    _init(out, peps)
    xlsx = tmp_path / "m.xlsx"
    _assignment_sheet(xlsx, [("B1", "Bladder", "AAA")])
    ok, msg = agent_core.build_patient_pools(str(out), str(xlsx), 0, 2)  # positional cols
    assert ok, msg
    assert _partial(out)["pools"][0]["member_peptide_ids"] == ["IMP-B1-AAA"]


def test_unmatched_sequences_reported_not_fatal(tmp_path):
    peps = [_pep("IMP-B1-AAA", "AAA")]
    out = tmp_path / "rec.json"
    _init(out, peps)
    xlsx = tmp_path / "m.xlsx"
    _assignment_sheet(xlsx, [("B1", "Bladder", "AAA"), ("B1", "Bladder", "ZZZ")])  # ZZZ not loaded
    ok, msg = agent_core.build_patient_pools(str(out), str(xlsx), "Patient Alias", "Peptide Sequence")
    assert ok, msg
    assert "unmatched" in msg
    assert _partial(out)["pools"][0]["member_peptide_ids"] == ["IMP-B1-AAA"]


def test_errors_when_no_peptides_loaded(tmp_path):
    out = tmp_path / "rec.json"
    _init(out, [])  # no peptides
    xlsx = tmp_path / "m.xlsx"
    _assignment_sheet(xlsx, [("B1", "Bladder", "AAA")])
    ok, msg = agent_core.build_patient_pools(str(out), str(xlsx), "Patient Alias", "Peptide Sequence")
    assert not ok
    assert "immunizing_peptides" in msg


def test_bad_column_name_errors(tmp_path):
    peps = [_pep("IMP-B1-AAA", "AAA")]
    out = tmp_path / "rec.json"
    _init(out, peps)
    xlsx = tmp_path / "m.xlsx"
    _assignment_sheet(xlsx, [("B1", "Bladder", "AAA")])
    ok, msg = agent_core.build_patient_pools(str(out), str(xlsx), "NopeCol", "Peptide Sequence")
    assert not ok
    assert "patient_col" in msg
