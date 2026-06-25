"""Deterministic mutant-vs-WT cross-reactivity loader (`build_crossreactivity_evidence`).

33064988 Supp Table 5 (mmc1.pdf) is a readable table — Peptide ID | Mutant | WT | Cross-reactive
to WT — but the agent hand-built its epitopes/evidence via add_entities and dropped them
run-to-run. This loads them deterministically: one MinimalEpitope + one epitope/immunogenic
evidence row per listed peptide, parent linked by sequence containment, mutation_specific from
the cross-reactivity column. FAITHFULNESS: a row whose parent peptide isn't loaded is skipped
(no orphan epitope invented).
"""
import json
import pathlib

import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]
REF = json.loads((PKT / "reference_records" / "rojas_extracted.json").read_text())
META = {k: v for k, v in REF.items() if not isinstance(v, list)}

# (peptide_id, mutant, wt, cross_reactive_to_WT)
ROWS = [
    ("M1-IM24-EPT0058", "IPYFQTINI", "IPYFQTKNI", False),          # 9mer -> class I, mutant-specific
    ("M1-IM01-ASP0037", "RVYPETYEWARKMAV", "RVHPETYEWARKMAV", True),  # 15mer -> class II, cross-reactive
    ("M4-IM16-ASP0028", "VAFTSQWHSTSSPVSM", "VAFTPQWHSTSSPVSM", False),
]


def _partial(out):
    return json.loads(pathlib.Path(str(out) + ".partial.json").read_text())


def _init(out, peptides):
    pathlib.Path(out).parent.mkdir(parents=True, exist_ok=True)
    ok, msg = agent_core.init_partial(str(out), json.dumps(META))
    assert ok, msg
    ok, msg = agent_core.append_section(str(out), "immunizing_peptides", json.dumps(peptides))
    assert ok, msg


def _pep(pid, seq):
    return {"paper_local_id": pid, "sequence": seq, "is_neoantigen": True,
            "quoted_text": "x", "section_ref": "Vaccine peptides"}


def _run(out, rows, monkeypatch):
    monkeypatch.setattr(agent_core, "_parse_crossreactivity_pdf", lambda _p: rows)
    return agent_core.build_crossreactivity_evidence(str(out), "ignored.pdf")


def test_builds_epitope_and_evidence_per_row(tmp_path, monkeypatch):
    out = tmp_path / "rec.json"
    # parents: the 9mer sits INSIDE a longer loaded peptide; the 15/16mers are loaded as-is
    _init(out, [_pep("IMP-M1-FQVKDIPYFQTINI", "FQVKDIPYFQTINI"),
                _pep("IMP-M1-RVYPETYEWARKMAV", "RVYPETYEWARKMAV"),
                _pep("IMP-M4-VAFTSQWHSTSSPVSM", "VAFTSQWHSTSSPVSM")])
    ok, msg = _run(out, ROWS, monkeypatch)
    assert ok, msg
    rec = _partial(out)
    eps = {e["paper_local_id"]: e for e in rec["epitopes"]}
    ev = {e["epitope_paper_id"]: e for e in rec["evidence"] if e["target_kind"] == "epitope"}
    assert len(eps) == 3 and len(ev) == 3
    # 9mer -> class I (+ lossless predicted_affinity); long -> class II
    assert eps["EP-M1-IM24-EPT0058"]["mhc_class"] == "I"
    assert eps["EP-M1-IM24-EPT0058"]["predicted_affinity"]["unit"] == "unknown"
    assert eps["EP-M1-IM01-ASP0037"]["mhc_class"] == "II"
    # parent linked by sequence containment (9mer -> the longer peptide that contains it)
    assert eps["EP-M1-IM24-EPT0058"]["parent_peptide_ids"] == ["IMP-M1-FQVKDIPYFQTINI"]
    # wild_type captured
    assert eps["EP-M1-IM24-EPT0058"]["wild_type_sequence"] == "IPYFQTKNI"
    # mutation_specific = NOT cross-reactive to WT
    assert ev["EP-M1-IM24-EPT0058"]["mutation_specific"] is True
    assert ev["EP-M1-IM01-ASP0037"]["mutation_specific"] is False
    # evidence shape
    e = ev["EP-M1-IM24-EPT0058"]
    assert e["outcome"] == "immunogenic" and e["assay"] == "elispot"
    assert e["patient_paper_id"] == "M1"
    assert e["provenance"][0]["kind"] == "table"
    assert e["magnitude"]["raw"]


def test_skips_row_without_a_loaded_parent(tmp_path, monkeypatch):
    out = tmp_path / "rec.json"
    _init(out, [_pep("IMP-M1-RVYPETYEWARKMAV", "RVYPETYEWARKMAV")])  # only the M1 ASP0037 parent
    ok, msg = _run(out, ROWS, monkeypatch)
    assert ok, msg
    eps = [e["paper_local_id"] for e in _partial(out)["epitopes"]]
    assert eps == ["EP-M1-IM01-ASP0037"]   # the 9mer + M4 row have no parent -> skipped
    assert "skipped" in msg


def test_deterministic_across_two_calls(tmp_path, monkeypatch):
    peps = [_pep("IMP-M1-FQVKDIPYFQTINI", "FQVKDIPYFQTINI"),
            _pep("IMP-M1-RVYPETYEWARKMAV", "RVYPETYEWARKMAV"),
            _pep("IMP-M4-VAFTSQWHSTSSPVSM", "VAFTSQWHSTSSPVSM")]
    seen = []
    for i in range(2):
        out = tmp_path / f"r{i}" / "rec.json"
        _init(out, peps)
        _run(out, ROWS, monkeypatch)
        seen.append(sorted((e["epitope_paper_id"], e["mutation_specific"])
                           for e in _partial(out)["evidence"] if e["target_kind"] == "epitope"))
    assert seen[0] == seen[1]


def test_errors_without_peptides(tmp_path, monkeypatch):
    out = tmp_path / "rec.json"
    pathlib.Path(out).parent.mkdir(parents=True, exist_ok=True)
    agent_core.init_partial(str(out), json.dumps(META))
    ok, msg = _run(out, ROWS, monkeypatch)
    assert not ok and "peptides" in msg


def test_real_pdf_parses_thirteen_rows():
    """Guarded smoke on the real mmc1.pdf (skipped if the corpus isn't present)."""
    import os
    corpus = os.environ.get("VAXTRACT_CORPUS_DIR")
    if not corpus:
        import pytest
        pytest.skip("set VAXTRACT_CORPUS_DIR to the corpus root to run this smoke test")
    pdf = pathlib.Path(corpus) / "33064988" / "supps" / "mmc1.pdf"
    if not pdf.exists():
        import pytest
        pytest.skip("33064988 mmc1.pdf not present")
    rows = agent_core._parse_crossreactivity_pdf(str(pdf))
    assert len(rows) == 13
    ids = {r[0] for r in rows}
    assert "M1-IM24-EPT0058" in ids
    # cross-reactivity parsed: only ASP0037 is cross-reactive to WT
    cross = {r[0]: r[3] for r in rows}
    assert cross["M1-IM01-ASP0037"] is True
    assert cross["M1-IM24-EPT0058"] is False
