#!/usr/bin/env python3
"""Per-paper batch wrapper for HPC-scale extraction (subscription auth).

One PMID in -> runs `extraction_agent.py`, validates the output, routes it, appends a ledger
row, and writes a terminal marker. The exit-code contract is what makes the Snakemake batch
RESUMABLE and quota-safe on a subscription:

  outcome        output written            marker?  exit  Snakemake effect
  -------------   -----------------------   -------  ----  ------------------------------------
  ok             results/extracted/PMID    yes      0     done; never re-run
  invalid        results/needs_review/PMID  yes     0     done; human reviews (NOT auto-retried)
  usage_limit    (none)                    no       1     retried on the next run (after window)
  error          (none)                    no       1     retried on the next run (transient/net)

The marker (results/markers/PMID.done) is the Snakemake rule output: it is written ONLY for
terminal outcomes (ok / invalid), so a quota stop or a network blip leaves NO marker and the
paper re-runs next time, while a genuinely-bad paper is NOT re-extracted on every invocation
(it sits in needs_review for a human). The ledger CSV is the scaled RUNS.md.

Pure helpers (parse_summary / usage_limit_hit / classify) are unit-tested without any network
or subprocess; main() is the I/O shell.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import pathlib
import re
import subprocess
import sys

# --- the runner's own summary line, e.g.
#   "[info] agent loop ok: subtype=success turns=27 cost=$4.08 errors=None"
_SUMMARY_RE = re.compile(
    r"agent loop (?:ok|ENDED IN ERROR):\s*subtype=(?P<subtype>\S+)\s+"
    r"turns=(?P<turns>\d+)\s+cost=\$(?P<cost>[\d.]+)")
# transient / quota signals that mean "stop and retry later", NOT "the data is bad"
_USAGE_LIMIT_RE = re.compile(
    r"usage limit|rate.?limit|\b429\b|limit reached|resets? at|overloaded|"
    r"too many requests|quota", re.I)


def parse_summary(log_text: str) -> dict:
    """Extract (subtype, turns, cost) from the runner's summary line. Empty dict if absent
    (e.g. the process died before printing it — a transient/abort case)."""
    m = _SUMMARY_RE.search(log_text or "")
    if not m:
        return {}
    return {"subtype": m.group("subtype"),
            "turns": int(m.group("turns")),
            "cost_usd": float(m.group("cost"))}


def usage_limit_hit(log_text: str) -> bool:
    """True if the log shows a subscription/rate/quota limit — the signal to retry later."""
    return bool(_USAGE_LIMIT_RE.search(log_text or ""))


def classify(out_exists: bool, valid: bool, log_text: str) -> tuple[str, bool]:
    """Map (did the runner write an output?, does it validate?, log) -> (status, retryable).

    retryable=True means: do NOT write a marker; let the batch re-run this PMID later (a quota
    stop or a transient failure). retryable=False means terminal (ok, or a data problem a human
    must look at) -> write the marker so quota isn't burned re-running it.
    """
    if out_exists and valid:
        return "ok", False
    if out_exists and not valid:
        return "invalid", False          # extractor finished but the record fails the schema
    if usage_limit_hit(log_text):
        return "usage_limit", True       # rolling-window / rate cap -> resume after it resets
    return "error", True                 # no output + no quota signal -> transient (net/crash)


# #5: PMID -> signed-off reference record (the only papers with a gold to score identity P/R against).
_GOLD_REF = {"37165196": "rojas", "30568305": "keskin", "33879241": "li"}
# Below this source-verification rate, the record has likely-hallucinated sequences -> route to review.
_FACT_RATE_FLOOR = 0.98


def _run_qc(repo: pathlib.Path, record_path: pathlib.Path, paper_dir: str, pmid: str) -> dict:
    """#5 per-paper QC gate. FACT verification (every extracted sequence/HLA must be in the source)
    runs on EVERY paper (no gold needed); identity precision/recall is added for the gold papers.
    Fail-soft: any error returns an empty dict so QC never breaks the extraction."""
    try:
        sys.path[:0] = [str(repo / "pr_compare")]
        import pr_compare
        rec = json.loads(record_path.read_text())
        gold = None
        name = _GOLD_REF.get(pmid)
        if name:
            gref = repo / "reference_records" / f"{name}_extracted.json"
            if gref.exists():
                gold = json.loads(gref.read_text())
        return pr_compare.qc_metrics(rec, paper_dir=paper_dir, gold=gold)
    except Exception as e:  # QC is advisory; never let it fail the run
        return {"qc_error": str(e)[:120]}


def _append_ledger(ledger_path: pathlib.Path, row: dict) -> None:
    fields = ["timestamp", "pmid", "status", "subtype", "turns", "cost_usd", "valid", "out_path",
              "fact_seq_rate", "fact_hla_rate", "evidence_precision", "evidence_recall", "note"]
    new = not ledger_path.exists()
    ledger_path.parent.mkdir(parents=True, exist_ok=True)
    with ledger_path.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        if new:
            w.writeheader()
        w.writerow({k: row.get(k, "") for k in fields})


def _now() -> str:
    # real HPC runtime (not the workflow sandbox) -> datetime is fine here
    from datetime import datetime, timezone
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="batch-extract one paper (PMID)")
    ap.add_argument("--pmid", required=True)
    ap.add_argument("--paper-dir", required=True, help="dir with the paper's pdf/xlsx")
    ap.add_argument("--repo", required=True, help="path to the extractor repo (holds extraction_agent.py)")
    ap.add_argument("--results", required=True, help="results root (extracted/ needs_review/ markers/ logs/)")
    ap.add_argument("--subscription", action="store_true", default=True,
                    help="use the claude subscription login (default; drops ANTHROPIC_API_KEY)")
    ap.add_argument("--api-key", dest="subscription", action="store_false",
                    help="use ANTHROPIC_API_KEY instead (for the future paid lane)")
    a = ap.parse_args(argv)

    repo = pathlib.Path(a.repo).resolve()
    results = pathlib.Path(a.results).resolve()
    extracted = results / "extracted" / f"{a.pmid}.json"
    review = results / "needs_review" / f"{a.pmid}.json"
    marker = results / "markers" / f"{a.pmid}.done"
    log_path = results / "logs" / f"{a.pmid}.log"
    ledger = results / "ledger.csv"
    for p in (extracted.parent, review.parent, marker.parent, log_path.parent):
        p.mkdir(parents=True, exist_ok=True)

    # the runner writes its final JSON straight to the extracted/ target; we relocate on invalid.
    cmd = [sys.executable, str(repo / "extraction_agent.py")]
    if a.subscription:
        cmd.append("--subscription")
    cmd += [str(pathlib.Path(a.paper_dir).resolve()), str(extracted)]

    proc = subprocess.run(cmd, cwd=str(repo), capture_output=True, text=True)
    log_text = (proc.stdout or "") + "\n" + (proc.stderr or "")
    log_path.write_text(log_text)

    valid = False
    if extracted.exists():
        sys.path[:0] = [str(repo), str(repo / "cancervac_packet")]
        import agent_core  # noqa: import after path injection
        valid, vmsg = agent_core.validate_record(extracted.read_text())
    else:
        vmsg = "no output file produced"

    status, retryable = classify(extracted.exists(), valid, log_text)
    summ = parse_summary(log_text)

    # v2.11.3: a VALID record that finalized by OVERRIDING a soft completeness/consistency guard is not
    # clean-silver — re-route it from 'ok' to 'needs_review' (terminal, like 'invalid': a human confirms
    # the override). The override list is persisted on the record by finalize_partial.
    overrides = []
    if status == "ok" and extracted.exists():
        try:
            overrides = json.loads(extracted.read_text()).get("finalize_overrides_used") or []
        except Exception:
            overrides = []

    # v2.11.5: TIER the overrides. A record whose overrides are ALL soft (funnel-size / regimen
    # divergence — expected variance / metadata, not data-correctness gaps) stays in the clean lane with
    # the overrides still recorded for audit; only a HARD override (recall/correctness gap, or any
    # unknown override) sends it to needs_review. Tiered after the 2026-06-09 batch sent all 7 papers to
    # needs_review, draining the clean lane of meaning.
    soft_only = status == "ok" and agent_core.overrides_are_soft_only(overrides)

    # #5: per-paper QC gate. Source-FACT verification (every extracted sequence/HLA must appear in the
    # paper) runs on EVERY produced record; identity P/R is added for the 3 gold papers. A record whose
    # sequences fall below the source-verification floor has likely HALLUCINATED data -> needs_review.
    qc = _run_qc(repo, extracted, str(pathlib.Path(a.paper_dir).resolve()), a.pmid) if extracted.exists() else {}
    fsr = qc.get("fact_seq_rate")
    hallucinated = status == "ok" and fsr is not None and fsr < _FACT_RATE_FLOOR

    note = ""
    if status == "invalid":
        review.parent.mkdir(parents=True, exist_ok=True)
        extracted.replace(review)        # move the bad record out of the clean lane
        note = f"validation: {vmsg[:160]}"
    elif hallucinated:
        status = "needs_review"
        review.parent.mkdir(parents=True, exist_ok=True)
        extracted.replace(review)        # sequences not found in source -> human checks for hallucination
        note = f"QC: {qc.get('fact_unverified')} extracted tokens not in source (seq-rate {fsr})"
    elif status == "ok" and overrides and not soft_only:
        status = "needs_review"
        review.parent.mkdir(parents=True, exist_ok=True)
        extracted.replace(review)        # HARD override -> out of the clean lane for a human
        note = f"finalize overrides used (HARD): {overrides}"
    elif status == "ok" and overrides:   # soft-only -> stays in extracted/, but noted
        note = f"clean; soft overrides only: {overrides}"
    elif status == "usage_limit":
        note = "subscription/rate limit — re-run after the window resets"
    elif status == "error":
        note = "no output + no quota signal — transient/network/crash"

    _append_ledger(ledger, {
        "timestamp": _now(), "pmid": a.pmid, "status": status,
        "subtype": summ.get("subtype", ""), "turns": summ.get("turns", ""),
        "cost_usd": summ.get("cost_usd", ""), "valid": valid,
        "out_path": str(review if status in ("invalid", "needs_review")
                        else extracted if status == "ok" else ""),
        "fact_seq_rate": qc.get("fact_seq_rate", ""), "fact_hla_rate": qc.get("fact_hla_rate", ""),
        "evidence_precision": qc.get("evidence_precision", ""),
        "evidence_recall": qc.get("evidence_recall", ""),
        "note": note})

    if retryable:
        # leave NO marker -> Snakemake re-runs this PMID on the next invocation
        print(f"[{a.pmid}] {status} (retryable): {note}", file=sys.stderr)
        return 1
    # terminal -> write the marker so the batch never re-extracts this paper
    marker.write_text(f"{status}\t{summ.get('subtype','')}\t{summ.get('cost_usd','')}\t{_now()}\n")
    print(f"[{a.pmid}] {status} -> marker written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
