#!/usr/bin/env python
"""lite_extract — assemble the extraction prompt for the lighter architecture.

Lighter than the MCP harness: schema digest + explicit cross-field RULES + NATIVE tools
(Read incl. PDF-image rendering, Bash+python, Grep, Write) + a validate-against-the-real-schema
repair loop as the "thin finalize". No custom MCP tools.

Usage:
    python lite_extract/build_prompt.py --pmid 38538867 \
        --paper-dir /abs/path/to/paper_dir --out /abs/path/out.json

Prints the full prompt to stdout (feed it to a native-tool agent: claude -p, the Agent tool, etc.).
"""
import argparse, pathlib, sys

HERE = pathlib.Path(__file__).resolve().parent
PKG = HERE.parent  # cancerVacExtrac_claudeAi

DIGEST = PKG / "outputs" / "vanilla_compare_2026-06-11" / "SCHEMA_DIGEST.txt"
RULES = HERE / "RULES.md"


def build(pmid: str, paper_dir: str, out: str) -> str:
    digest = DIGEST.read_text() if DIGEST.exists() else "(schema digest missing — regenerate with schema_digest.build_schema_digest)"
    rules = RULES.read_text()
    return f"""\
You are extracting one cancer-vaccine paper into a single schema-valid JSON record. Use ONLY native
tools: Read, Bash, Grep, Glob, Write. There are NO custom extraction tools — you build the record yourself.

PMID = {pmid}
PAPER SOURCES = {paper_dir}

=== HOW TO READ THE SOURCES (native, image-first) ===
- Main text / PDFs: the **Read tool renders PDF pages as images** — use it. Many manifests and tables
  are IMAGE-LOCKED (printed in a figure or an "Extended Data Table" image, NOT in any spreadsheet).
  Read the relevant PDF pages as images and transcribe the table. THIS IS A FIRST-CLASS PATH, not a
  fallback — the most common recall miss is a peptide/epitope table that lives only in an image.
- .xlsx/.xls/.csv supplements: parse with Bash via {sys.executable} (openpyxl + pandas available).
  List every sheet of every workbook FIRST, then read the ones with the manifest / per-patient assays.
- .docx: parse with python (python-docx or unzip + XML).
- Inventory the whole directory (Glob/ls) before deciding what holds the manifest and the immunogenicity.

=== TARGET SCHEMA (v2.12.0) ===
{digest}

=== CROSS-FIELD RULES (these are what make the record VALID — read carefully) ===
{rules}

=== COMPLETENESS ===
Capture every patient, every immunizing peptide, every epitope, every immunogenicity measurement
(positive AND reported-negative), pools, survival. For big tables, WRITE A PYTHON SCRIPT to convert
rows → JSON entities; do not transcribe by hand. Do not consult any reference/gold answer.

=== THIN FINALIZE (the validator IS the contract) ===
After writing {out}, VALIDATE it against the real schema and REPAIR until it passes:
    {sys.executable} - <<'PY'
    import sys; sys.path[:0]=['{PKG}', '{PKG}/cancervac_packet']
    import agent_core
    ok, msg = agent_core.validate_record(open('{out}').read())
    print('PASS' if ok else msg)
    PY
Loop: read the errors, fix the JSON, re-validate. Do not finish until it prints PASS (or until you can
justify in curator_notes why a record legitimately stays minimal, e.g. a companion-deferred manifest).

FINAL MESSAGE: one line `COUNTS pat=.. IMP=.. epi=.. ev=.. pool=.. nm=..`, whether the validator PASSES,
and 2-3 sentences on where the manifest lived (xlsx vs image-locked vs companion-deferred).
"""


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmid", required=True)
    ap.add_argument("--paper-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    print(build(a.pmid, a.paper_dir, a.out))


if __name__ == "__main__":
    main()
