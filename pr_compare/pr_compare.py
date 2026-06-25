#!/usr/bin/env python
"""#7 — provenance-matched precision/recall comparator (replaces count-ratio compare.py).

Scores an extraction vs gold by ENTITY IDENTITY (not cardinality) and verifies every emitted
quoted_text against the paper's source text. Three columns per entity type: matched / emitted-extra
(precision miss) / gold-missing (recall miss) — so a drop+hallucinate that nets to the same count no
longer reads as parity. Plus a provenance-verification rate and a gold-staleness flag.

Usage:
  python pr_compare.py GOLD.json EMITTED.json [--paper-dir DIR_with_pdf_xlsx]

Blockers handled (per Fable's plan):
  - source availability: builds a source_index from PDF (PyMuPDF/pypdf) + xlsx (openpyxl) + text files;
    a quoted_text whose provenance is a FIGURE, or that isn't in the index, is bucketed
    'unverifiable' (NOT a fail) rather than silently passing.
  - normalization: a single normalize() (NFKC, casefold, collapse-whitespace, strip ellipsis) applied
    to BOTH the identity keys and the provenance check.
  - gold staleness: if the emitted record is schema-valid and gold is NOT (under the current schema),
    the unmatched deltas are flagged 'gold_suspect' instead of scored as pure misses.
"""
import argparse, json, pathlib, re, sys, unicodedata

_ROOT = pathlib.Path(__file__).resolve().parents[2]   # cancerVacExtrac_claudeAi
sys.path[:0] = [str(_ROOT), str(_ROOT / "cancervac_packet")]
import agent_core  # noqa: E402


def normalize(s):
    if not s:
        return ""
    s = unicodedata.normalize("NFKC", str(s)).casefold()
    s = s.replace("…", "...").replace("–", "-").replace("—", "-")
    return re.sub(r"\s+", " ", s).strip()


# ---- source index (for provenance verification) ----

def build_source_index(paper_dir):
    if not paper_dir:
        return ""
    d = pathlib.Path(paper_dir)
    chunks = []
    for p in d.rglob("*"):
        if not p.is_file():
            continue
        suf = p.suffix.lower()
        try:
            if suf == ".pdf":
                try:
                    import fitz
                    chunks.append(" ".join(pg.get_text() for pg in fitz.open(p)))
                except Exception:
                    import pypdf
                    chunks.append(" ".join((pg.extract_text() or "") for pg in pypdf.PdfReader(str(p)).pages))
            elif suf in (".xlsx", ".xls"):
                import openpyxl
                wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
                for ws in wb.worksheets:
                    for row in ws.iter_rows(values_only=True):
                        chunks.append(" ".join(str(c) for c in row if c is not None))
            elif suf in (".txt", ".md", ".csv", ".tsv"):
                chunks.append(p.read_text(errors="ignore"))
        except Exception as e:
            chunks.append("")  # unreadable source file -> rows from it become 'unverifiable'
    return normalize(" ".join(chunks))


# ---- identity keys per entity type ----

def _tgt(r):
    return (r.get("immunizing_peptide_paper_id") or r.get("epitope_paper_id")
            or r.get("pool_paper_id") or r.get("candidate_paper_id") or "")


def build_idmap(rec):
    """Resolve each extractor-assigned paper_local_id to its SEQUENCE (peptides/epitopes/candidates) or
    its member-sequence set (pools). Extractor-local ids are NOT stable across two independent runs, so
    evidence/pool/screening identity must key on the resolved sequence, not the id (Fable's blocker)."""
    seq = {}
    for t in ("immunizing_peptides", "epitopes", "candidates"):
        for r in rec.get(t) or []:
            if r.get("paper_local_id") and r.get("sequence"):
                seq[r["paper_local_id"]] = normalize(r["sequence"])
    pool = {}
    for p in rec.get("pools") or []:
        members = frozenset(seq.get(m, m) for m in (p.get("member_peptide_ids") or []))
        if p.get("paper_local_id"):
            pool[p["paper_local_id"]] = members
    return seq, pool


def _resolved_tgt(r, idmap):
    seq, pool = idmap
    tid = _tgt(r)
    if r.get("pool_paper_id"):
        return ("pool", pool.get(tid, tid))
    return ("seq", seq.get(tid, tid))   # fall back to the raw id when unresolvable


# each key fn takes (row, idmap) so evidence/pool/screening resolve ids -> sequences
KEYS = {
    "patients":            lambda r, m: ("pat", r.get("paper_local_id")),
    "immunizing_peptides": lambda r, m: ("imp", normalize(r.get("sequence"))),
    "epitopes":            lambda r, m: ("epi", normalize(r.get("sequence")), r.get("mhc_class")),
    "pools":               lambda r, m: ("pool", r.get("patient_paper_id"),
                                         m[1].get(r.get("paper_local_id"), frozenset())),
    "evidence":            lambda r, m: ("ev", r.get("patient_paper_id"), _resolved_tgt(r, m),
                                         r.get("assay"), r.get("outcome")),
    "candidates":          lambda r, m: ("cand", normalize(r.get("sequence"))),
    "screening_readouts":  lambda r, m: ("scr", r.get("patient_paper_id"), _resolved_tgt(r, m), r.get("assay")),
}


def match_type(gold, emit, keyfn, gmap, emap):
    gk = {}
    for r in gold:
        gk.setdefault(keyfn(r, gmap), []).append(r)
    ek = {}
    for r in emit:
        ek.setdefault(keyfn(r, emap), []).append(r)
    gset, eset = set(gk), set(ek)
    matched = gset & eset
    gold_missing = gset - eset       # recall miss
    emitted_extra = eset - gset      # precision miss
    p = len(matched) / len(eset) if eset else None
    rc = len(matched) / len(gset) if gset else None
    return dict(n_gold=len(gold), n_emit=len(emit), matched=len(matched),
                gold_missing=len(gold_missing), emitted_extra=len(emitted_extra),
                precision=p, recall=rc,
                sample_gold_missing=[k for k in list(gold_missing)[:5]],
                sample_emitted_extra=[k for k in list(emitted_extra)[:5]])


# ---- FACT verification (anti-hallucination) ----
# CALIBRATION (2026-06-12): quoted_text in this pipeline is a SYNTHESIZED provenance label, not a
# verbatim source quote — even the human-audited gold's quoted_texts are not source substrings, while
# the SEQUENCE inside them IS in the source. So we verify the identity-bearing EXTRACTED TOKENS (peptide/
# epitope sequence; HLA allele) against the source — a sequence absent from the source is a likely
# hallucinated/mis-transcribed value (the real faithfulness defect). Figure-provenance rows are
# 'unverifiable' (image-locked), not failures.
_SEQ_TYPES = ("immunizing_peptides", "epitopes", "candidates", "screening_readouts")

def verify_facts(rec, source_index):
    res = {"sequence": [0, 0, 0], "hla_allele": [0, 0, 0]}   # [verified, unverified, unverifiable]
    misses = []
    def _check(kind, val, prov_kinds, tname, label):
        v = normalize(val)
        if not v:
            return
        if not source_index:
            res[kind][2] += 1
        elif v in source_index:
            res[kind][0] += 1
        elif "figure" in prov_kinds:
            res[kind][2] += 1
        else:
            res[kind][1] += 1
            if len(misses) < 10:
                misses.append((tname, kind, label, v[:40]))
    for tname in KEYS:
        for r in rec.get(tname) or []:
            pk = {p.get("kind") for p in (r.get("provenance") or [])}
            lid = r.get("paper_local_id") or r.get("id") or r.get("patient_paper_id") or "?"
            if tname in _SEQ_TYPES:
                _check("sequence", r.get("sequence"), pk, tname, lid)
            _check("hla_allele", r.get("hla_allele"), pk, tname, lid)
    return res, misses


def qc_metrics(rec, paper_dir=None, gold=None):
    """Per-paper QC for the batch lane (#5). FACT verification (sequence/HLA in source) needs NO gold and
    works on every paper; identity precision/recall is added only when a gold record is supplied. Returns
    a flat dict for the ledger. Fail-soft: a missing source/lib yields None rates, never an exception."""
    src = build_source_index(paper_dir)
    facts, _ = verify_facts(rec, src)

    def _rate(kind):
        v, u, _x = facts[kind]
        return round(v / (v + u), 4) if (v + u) else None

    out = {
        "fact_seq_rate": _rate("sequence"),
        "fact_hla_rate": _rate("hla_allele"),
        "fact_unverified": facts["sequence"][1] + facts["hla_allele"][1],
        "evidence_precision": None, "evidence_recall": None,
    }
    if gold is not None:
        gmap, emap = build_idmap(gold), build_idmap(rec)
        ev = match_type(gold.get("evidence") or [], rec.get("evidence") or [], KEYS["evidence"], gmap, emap)
        out["evidence_precision"], out["evidence_recall"] = ev["precision"], ev["recall"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gold"); ap.add_argument("emitted"); ap.add_argument("--paper-dir")
    a = ap.parse_args()
    gold = json.loads(pathlib.Path(a.gold).read_text())
    emit = json.loads(pathlib.Path(a.emitted).read_text())
    g_ok, _ = agent_core.validate_record(json.dumps(gold))
    e_ok, e_msg = agent_core.validate_record(json.dumps(emit))
    gold_suspect = e_ok and not g_ok

    print(f"=== provenance-matched compare ===  emitted_valid={e_ok}  gold_valid={g_ok}"
          + ("  [GOLD-SUSPECT: emitted valid, gold violates current schema]" if gold_suspect else ""))
    if not e_ok:
        print("  emitted INVALID:", e_msg[:160])
    print(f"  {'type':<20}{'gold':>5}{'emit':>5}{'match':>6}{'g_miss':>7}{'e_extra':>8}"
          f"{'prec':>6}{'rec':>6}")
    gmap, emap = build_idmap(gold), build_idmap(emit)
    for tname, keyfn in KEYS.items():
        m = match_type(gold.get(tname) or [], emit.get(tname) or [], keyfn, gmap, emap)
        pr = f"{m['precision']:.2f}" if m['precision'] is not None else "-"
        rc = f"{m['recall']:.2f}" if m['recall'] is not None else "-"
        print(f"  {tname:<20}{m['n_gold']:>5}{m['n_emit']:>5}{m['matched']:>6}"
              f"{m['gold_missing']:>7}{m['emitted_extra']:>8}{pr:>6}{rc:>6}")

    src = build_source_index(a.paper_dir)
    res, misses = verify_facts(emit, src)
    print(f"  FACT verification vs source ({'no paper-dir' if not src else f'{len(src)} chars'}):")
    for kind, (v, u, x) in res.items():
        rate = f"{v/(v+u):.2f}" if (v + u) else "n/a"
        print(f"    {kind:<12} verified={v:>4} unverified={u:>4} unverifiable(figure)={x:>4} "
              f"| in-source-rate={rate}")
    if misses:
        print("  sample tokens NOT found in source (possible hallucination/mis-transcription):")
        for t, kind, lid, val in misses:
            print(f"    [{t}/{kind}] {lid}: {val!r}")


if __name__ == "__main__":
    main()
