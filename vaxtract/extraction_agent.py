#!/usr/bin/env python3
"""extraction_agent.py — antVacDB paper -> validated JSON, as a Claude Agent SDK agent.

WHAT MAKES THIS AN AGENT: the model is given TOOLS and a LOOP. It reads the
tables/text itself, drafts candidate JSON, calls validate(), reads the schema's
errors, and FIXES its own output — repeating until the record passes the packet
schema or it gives up (-> outer-guard quarantine).

Pure logic lives in agent_core.py / prompt_render.py (SDK-free, unit-tested).
This file is the SDK shell: it wraps that logic as in-process MCP tools, builds
the options, and runs the loop.

Setup (the agent lives in the ``[agent]`` extra; core install is schema-only):
    python -m pip install 'vaxtract[agent]'
    export ANTHROPIC_API_KEY=...        # the SDK's primary auth
Usage:
    vaxtract <paper_dir> [out_path]   # console entry point (see cli.py)

Author: Samuel Ahuno (ekwame001@gmail.com)
"""
from __future__ import annotations

import os
import pathlib
import sys

from . import agent_core
from .prompt_render import build_system_prompt, field_guidance_only
from .schema_digest import build_schema_digest

try:
    from claude_agent_sdk import (
        ClaudeAgentOptions,
        ClaudeSDKClient,
        PermissionResultAllow,
        PermissionResultDeny,
        create_sdk_mcp_server,
        tool,
    )
except ModuleNotFoundError as exc:  # core-only install (no `[agent]` extra)
    # Narrow to ModuleNotFoundError on purpose: a plain ImportError (e.g. a future
    # pre-1.0 SDK renames one of these names) should surface its real message, not
    # this install hint — the SDK is present, so `pip install [agent]` wouldn't help.
    raise ImportError(
        "vaxtract's extraction agent requires the optional 'agent' dependencies "
        "(the Claude Agent SDK + file readers). Install them with:\n"
        "    pip install 'vaxtract[agent]'\n"
        "The schema/vocab data contract (`import vaxtract.schema`) works without them."
    ) from exc

# 1M-context variant (B, 2026-06-04): Keskin overflowed the 200K window -> 2 compactions
# -> re-orientation churn (extra turns/thinking). The [1m] model removes compaction entirely.
# Trade-off: 1M-tier pricing applies only ABOVE 200K, so small papers (Rojas, ~64K/turn) are
# unaffected; only a paper big enough to have compacted pays the premium -- exactly when it
# helps. (The SDK `betas=['context-1m-2025-08-07']` flag is Sonnet-only and does NOT enable
# 1M on Opus; the context window is selected via this model id.)
MODEL = "claude-opus-4-8[1m]"
# 120: add_table made the bulk fast, so the remaining budget is the judgment/reconciliation
# tail (evidence, patient counts, survival, finalize). Run 4 hit the old cap of 60 (RUNS.md).
MAX_TURNS = 120
# Runaway-cost backstop (run4 hit $13.34). Generous enough for a normal completion; caps a
# pathological run rather than letting it blow past unbounded. The on-disk partial survives.
MAX_BUDGET_USD = 18.0


# ---- TOOLS: thin wrappers over agent_core ----------------------------------
@tool("read_table",
      "Parse an .xlsx supplementary/source-data file into rows (Layer 1, highest fidelity). "
      "Returns up to max_rows rows and reports the total, so you know if a table was truncated. "
      "To INSPECT or SLICE a big table (count/find rows by a column value), use row_filter and/or "
      "columns HERE -- do NOT try to grep or read the tool's spilled result file (host file tools "
      "are unavailable). A huge sheet is byte-capped (you still get the TRUE total) -- narrow with "
      "row_filter/columns/max_rows rather than chasing a spill. Leading 'Table Sn.' title rows are "
      "auto-skipped for column projection. Set underline=true to reveal minimal-epitope substrings "
      "marked by underlining (returned wrapped in <u>…</u>). To BULK-ADD rows to the record, use add_table.",
      {"type": "object",
       "properties": {
           "path": {"type": "string", "description": "path to the .xlsx file"},
           "sheet": {"type": "string", "description": "sheet name; omit to use the first sheet"},
           "max_rows": {"type": "integer",
                        "description": "max rows to return (default 500); raise it to read a larger table fully"},
           "row_filter": {"type": "object",
                          "description": "keep only matching rows: {\"col\":H | \"col_idx\":N | \"col_letter\":\"L\", "
                                         "\"in\"|\"equals\"|\"not_empty\":...}; reports matched/total counts. Use "
                                         "col_idx/col_letter for blank/duplicated-header columns names can't reach."},
           "columns": {"type": "array", "items": {"type": ["string", "integer"]},
                       "description": "project to just these columns -- a header NAME (string) or a 0-based "
                                      "POSITION (integer) for nameless columns; output is keyed by name or [idx]"},
           "underline": {"type": "boolean",
                         "description": "reveal underlined sub-sequences (e.g. minimal epitopes marked by "
                                        "underlining inside a longer peptide) wrapped in <u>…</u>"},
           "header_row": {"type": "integer",
                          "description": "0-based header row; omit to auto-skip leading title rows"},
       },
       "required": ["path"]})
async def read_table(args):
    try:
        text = agent_core.read_table_rows(
            args["path"], args.get("sheet") or None, int(args.get("max_rows") or 500),
            row_filter=args.get("row_filter") or None, columns=args.get("columns") or None,
            underline=bool(args.get("underline")),
            header_row=args.get("header_row") if args.get("header_row") is not None else None)
    except Exception as e:
        text = f"ERROR reading table {args.get('path')!r}: {e}"
    return {"content": [{"type": "text", "text": text}]}


@tool("read_docx",
      "Parse a .docx supplement — common for BMC/Molecular-Cancer & Frontiers (Data_Sheet_*.docx) and "
      "often holding per-entity data the xlsx/pdf readers miss. THREE modes: (a) omit table_index -> a "
      "SUMMARY listing each table's rows×cols + its caption (the preceding 'Supplementary Table Sn…' "
      "paragraph); (b) table_index=N -> that table (0-based) as a capped preview, with the SAME "
      "row_filter/columns/underline/byte-cap behaviour as read_table -> bulk-load it with add_table just "
      "like an xlsx; (c) text_offset=0 -> the paragraph PROSE in a paged window (next-offset reported). "
      "Do NOT read a spilled result file (host file tools are unavailable).",
      {"type": "object",
       "properties": {
           "path": {"type": "string", "description": "path to the .docx file"},
           "table_index": {"type": "integer", "description": "0-based table to read; omit for the table summary"},
           "max_rows": {"type": "integer", "description": "max rows to return (default 500)"},
           "row_filter": {"type": "object",
                          "description": "keep only matching rows: {\"col\":H | \"col_idx\":N | \"col_letter\":\"L\", "
                                         "\"in\"|\"equals\"|\"not_empty\":...} (same DSL as read_table)"},
           "columns": {"type": "array", "items": {"type": ["string", "integer"]},
                       "description": "project to these columns -- header NAME or 0-based index"},
           "underline": {"type": "boolean", "description": "reveal underlined sub-sequences wrapped in <u>…</u>"},
           "text_offset": {"type": "integer", "description": "read PROSE from this char offset instead of a table"},
           "max_chars": {"type": "integer", "description": "prose window size (default 40000)"},
       },
       "required": ["path"]})
async def read_docx(args):
    try:
        text = agent_core.read_docx_from(
            args["path"],
            table_index=args.get("table_index") if args.get("table_index") is not None else None,
            max_rows=int(args.get("max_rows") or 500),
            row_filter=args.get("row_filter") or None, columns=args.get("columns") or None,
            underline=bool(args.get("underline")),
            text_offset=args.get("text_offset") if args.get("text_offset") is not None else None,
            max_chars=int(args.get("max_chars") or 40_000))
    except Exception as e:
        text = f"ERROR reading docx {args.get('path')!r}: {e}"
    return {"content": [{"type": "text", "text": text}]}


@tool("read_pdf_text",
      "Extract selectable text from the main or a supplementary PDF (methods, legends, results). "
      "Returns a paging window (default 40000 chars from offset 0); if more text remains the result "
      "reports the next offset -- call again with that offset to continue. Do NOT read a spilled "
      "result file (host file tools are unavailable).",
      {"type": "object",
       "properties": {
           "path": {"type": "string"},
           "offset": {"type": "integer", "description": "start char (default 0); use the next-offset the result reports"},
           "max_chars": {"type": "integer", "description": "window size (default 40000)"},
       },
       "required": ["path"]})
async def read_pdf_text(args):
    try:
        text = agent_core.read_pdf_text_from(
            args["path"], int(args.get("max_chars") or 40_000), int(args.get("offset") or 0))
    except Exception as e:
        text = f"ERROR reading pdf {args.get('path')!r}: {e}"
    return {"content": [{"type": "text", "text": text}]}


@tool("survey_sources",
      "INVENTORY every supplement in ONE call so you can LOCATE where data lives before opening files. "
      "Pass the paper directory (recursed) or a single file. Returns, per .xlsx/.pdf/.docx: xlsx -> each "
      "sheet's name, rows x cols, and HEADER ROW (so a peptide/immunogenicity table is recognizable by "
      "its columns even when the file/sheet name is unhelpful); pdf -> page count + first text "
      "(flags image-only figure PDFs); docx -> table count + caption + first table header. "
      "CALL THIS FIRST on a paper with many supplements (don't guess which file holds the manifest). "
      "Byte-capped; dropped files are listed (re-call on a subdir/file to see them).",
      {"type": "object",
       "properties": {
           "path": {"type": "string", "description": "paper directory (recursed) or a single supplement file"},
           "max_chars": {"type": "integer", "description": "byte cap for the digest (default 14000)"},
       },
       "required": ["path"]})
async def survey_sources(args):
    try:
        text = agent_core.survey_sources(args["path"], int(args.get("max_chars") or 14000))
    except Exception as e:
        text = f"ERROR surveying {args.get('path')!r}: {e}"
    return {"content": [{"type": "text", "text": text}]}


@tool("read_figure",
      "Read numbers off a figure when they exist in NO table/source-data/text. TWO STEPS: "
      "(1) call with path+page (no region) to SEE the page image and locate the panel; "
      "(2) call again with region=[x0,y0,x1,y1] (fractions of the page, top-left origin) "
      "to get a zoomed, legible crop, then read the values. RECORD conservatively: put the "
      "read estimate in Measurement.raw, set value=null (or a number only for a clean simple "
      "chart), tier='reported', confidence<=2; attach Provenance(kind='figure', "
      "needs_review=true); set the row's quoted_text to a VERBATIM fragment of the figure "
      "CAPTION (get it via read_pdf_text) and section_ref to the panel (e.g. 'Figure 1G'). "
      "NEVER fabricate: if a value is unreadable, leave value=null.",
      {"type": "object",
       "properties": {"path": {"type": "string"},
                      "page": {"type": "integer"},
                      "what": {"type": "string"},
                      "region": {"type": "array", "items": {"type": "number"},
                                 "description": "optional [x0,y0,x1,y1] fractions of the page"}},
       "required": ["path", "page", "what"]})
async def read_figure(args):
    try:
        region = args.get("region")
        region = tuple(region) if region else None
        b64, w, h = agent_core.render_figure_image(args["path"], args["page"], region=region)
    except Exception as e:
        return {"content": [{"type": "text", "text":
                f"ERROR rendering figure {args.get('path')!r} p{args.get('page')}: {e}"}]}
    note = (f"Rendered {'region ' + str(region) if region else 'full page'} "
            f"of {args['path']} page {args['page']} ({w}x{h}). Read the values for "
            f"{args.get('what')!r}. If this is the full page, locate the panel and call "
            f"read_figure again with region=[x0,y0,x1,y1] for a legible crop. RECORD: "
            f"value=null (or a number only for a clean simple chart) + estimate in raw, "
            f"tier='reported', confidence<=2, Provenance(kind='figure', needs_review=true), "
            f"quoted_text = a verbatim figure-caption fragment.")
    return {"content": [
        {"type": "image", "data": b64, "mimeType": "image/png"},
        {"type": "text", "text": note},
    ]}


@tool("validate",
      "(DEPRECATED — small records only; prefer init_record/add_entities/finalize) Check a candidate ExtractedPaper JSON against the schema. Returns 'VALID' or the errors to fix.",
      {"candidate_json": str})
async def validate(args):
    _, msg = agent_core.validate_record(args["candidate_json"])
    return {"content": [{"type": "text", "text": msg}]}


@tool("save_extraction",
      "(DEPRECATED — small records only; prefer init_record/add_entities/finalize) Validate AND write the final JSON in one shot. Fails if it exceeds the output limit on large papers.",
      {"candidate_json": str, "out_path": str})
async def save_extraction(args):
    _, msg = agent_core.save_record(args["candidate_json"], args["out_path"])
    return {"content": [{"type": "text", "text": msg}]}


@tool("init_record",
      "Start a new record from paper-level fields (pmid, journal, year, title, cohort_size, "
      "indication_summary, + optional doi/pmcid/nct_id/n_enrolled). Entity lists are added later "
      "with add_entities. Creates an on-disk partial that survives across turns.",
      {"type": "object",
       "properties": {"out_path": {"type": "string"},
                      "paper_meta_json": {"type": "string", "description": "JSON object of paper-level fields"}},
       "required": ["out_path", "paper_meta_json"]})
async def init_record(args):
    _, msg = agent_core.init_partial(args["out_path"], args["paper_meta_json"])
    return {"content": [{"type": "text", "text": msg}]}


@tool("add_entities",
      "Append a batch (<=~50) of entities to one section of the in-progress record. section is one "
      "of: patients, immunizing_peptides, epitopes, pools, evidence, survival_outcomes. Each item is "
      "validated against its schema; if any item is invalid the WHOLE batch is rejected and the "
      "partial is unchanged. Call repeatedly to add all rows (do NOT omit non-immunogenic peptides). "
      "For evidence, add one row per REPORTED response, not one per peptide. For bulk peptide tables "
      "prefer add_table and omit patient_paper_id.",
      {"type": "object",
       "properties": {"out_path": {"type": "string"},
                      "section": {"type": "string"},
                      "items_json": {"type": "string", "description": "JSON list of entity objects"}},
       "required": ["out_path", "section", "items_json"]})
async def add_entities(args):
    _, msg = agent_core.append_section(args["out_path"], args["section"], args["items_json"])
    return {"content": [{"type": "text", "text": msg}]}


@tool("set_safety_summary",
      "Set the PAPER-LEVEL safety_summary (CTCAE-grade headline toxicity facts). This is the ONLY "
      "way to record safety -- it is a scalar, not an entity list, so add_entities does NOT take it. "
      "Fields: max_related_grade (1-5, highest treatment-RELATED CTCAE grade -- NOT the immunogenicity "
      "response grade), any_grade3plus_related (bool), n_patients_with_related_ae (int), irae_present "
      "(bool), raw (verbatim safety sentence). Omit any field the paper doesn't state. Re-callable "
      "(overwrites). MOST vaccine trials report safety -- read the safety paragraph / AE table and set "
      "this before finalize. SERIOUSNESS != GRADE: a 'serious adverse event'/SAE is a REGULATORY "
      "category (death/hospitalization/life-threatening/disability), NOT a CTCAE grade -- never set "
      "any_grade3plus_related from 'serious'/'SAE'; set it true ONLY for a grade>=3 (or 'severe'/'life-"
      "threatening') AE the paper attributes to the TREATMENT, NOT to disease progression. The `raw` "
      "you pass must itself contain the grade + relatedness you assert (finalize cross-checks it). "
      "Pass safety_json='null' to clear it.",
      {"type": "object",
       "properties": {"out_path": {"type": "string"},
                      "safety_json": {"type": "string", "description": "JSON object of SafetySummary fields"}},
       "required": ["out_path", "safety_json"]})
async def set_safety_summary(args):
    _, msg = agent_core.set_safety_summary(args["out_path"], args["safety_json"])
    return {"content": [{"type": "text", "text": msg}]}


@tool("clear_entities",
      "Reset one section to empty (to correct a mistake), then re-add it with add_entities.",
      {"type": "object",
       "properties": {"out_path": {"type": "string"}, "section": {"type": "string"}},
       "required": ["out_path", "section"]})
async def clear_entities(args):
    _, msg = agent_core.clear_section(args["out_path"], args["section"])
    return {"content": [{"type": "text", "text": msg}]}


@tool("partial_status",
      "Report current per-section counts + the paper metadata. Use to check progress, or to resume "
      "after a context summary (the on-disk partial is the source of truth, not your memory).",
      {"type": "object", "properties": {"out_path": {"type": "string"}}, "required": ["out_path"]})
async def partial_status(args):
    _, msg = agent_core.partial_status(args["out_path"])
    return {"content": [{"type": "text", "text": msg}]}


@tool("finalize",
      "Validate the assembled record against the FULL schema and write it. If it reports errors "
      "(e.g. count reconciliation, orphaned epitopes, unknown patient ref), fix with "
      "clear_entities/add_entities and call finalize again. It also blocks ONCE if immunogenic "
      "ELISpot/ICS/etc. responses have no magnitude -- add magnitudes from the figure source-data, "
      "or pass allow_missing_magnitudes=true to proceed if none are reported. It also blocks ONCE if "
      "a patient has pooled evidence but no peptide-pool entity -- add the pool, or pass "
      "allow_missing_pools=true. It also blocks ONCE if a pooled response was encoded as per-member "
      "rows instead of ONE target_kind='pool' row -- collapse them to one pool row, or pass "
      "allow_member_level_pool_evidence=true if the source deconvolutes every member. For the "
      "candidate funnel it blocks ONCE if candidates exist "
      "but n_predicted_reported is unset -- set n_predicted_reported/n_selected_reported to "
      "the counts the paper states, or pass allow_unknown_funnel_size=true; and ONCE if a "
      "candidate's selected_peptide_id bridges to an IMP with a different sequence -- fix the "
      "bridge, or pass allow_candidate_bridge_mismatch=true. NEVER clear_entities the "
      "candidates to get past these -- that destroys the funnel; set the field or pass the flag. "
      "Terminal step.",
      {"type": "object",
       "properties": {"out_path": {"type": "string"},
                      "allow_missing_magnitudes": {"type": "boolean",
                          "description": "set true to proceed when some responses have no reported magnitude"},
                      "allow_missing_pools": {"type": "boolean",
                          "description": "set true to proceed when a patient has pooled evidence but pool membership is unresolvable"},
                      "allow_member_level_pool_evidence": {"type": "boolean",
                          "description": "set true to proceed when a pooled response is kept as per-member rows because the source deconvolutes every member"},
                      "allow_unknown_funnel_size": {"type": "boolean",
                          "description": "set true to proceed when candidates exist but the paper reports no predicted total"},
                      "allow_candidate_bridge_mismatch": {"type": "boolean",
                          "description": "set true to proceed when a candidate's selected_peptide_id IMP has a different sequence"},
                      "allow_regimen_divergence": {"type": "boolean",
                          "description": "set true to proceed when patients of one arm have different delivery regimens (e.g. dose escalation)"},
                      "allow_evidence_count_mismatch": {"type": "boolean",
                          "description": "set true to proceed when recorded evidence differs from the paper's stated immunogenic/negative counts because they are at a different grain (e.g. pooled)"},
                      "allow_peptide_count_mismatch": {"type": "boolean",
                          "description": "set true to proceed when recorded peptides fall short of n_selected_reported because the paper's count includes peptides not individually listed"},
                      "allow_sparse_evidence": {"type": "boolean",
                          "description": "set true to proceed when immune evidence covers only a subset of vaccinated patients because the trial immune-monitored only that subset"},
                      "allow_missing_class_ii": {"type": "boolean",
                          "description": "proceed despite a CD4/class-II cue with no class-II record minted"},
                      "allow_missing_minimal_epitopes": {"type": "boolean",
                          "description": "proceed despite a mislabeled long-peptide epitope or a dropped predicted minimal-epitope layer"},
                      "allow_ungrounded_safety_grade": {"type": "boolean",
                          "description": "proceed when a grade>=3 treatment-related safety claim cannot be grounded in a verbatim grade>=3 raw quote (routes to needs_review)"}},
       "required": ["out_path"]})
async def finalize(args):
    _, msg = agent_core.finalize_partial(
        args["out_path"], allow_missing_magnitudes=bool(args.get("allow_missing_magnitudes")),
        allow_missing_pools=bool(args.get("allow_missing_pools")),
        allow_member_level_pool_evidence=bool(args.get("allow_member_level_pool_evidence")),
        allow_unknown_funnel_size=bool(args.get("allow_unknown_funnel_size")),
        allow_candidate_bridge_mismatch=bool(args.get("allow_candidate_bridge_mismatch")),
        allow_regimen_divergence=bool(args.get("allow_regimen_divergence")),
        allow_evidence_count_mismatch=bool(args.get("allow_evidence_count_mismatch")),
        allow_peptide_count_mismatch=bool(args.get("allow_peptide_count_mismatch")),
        allow_sparse_evidence=bool(args.get("allow_sparse_evidence")),
        allow_missing_class_ii=bool(args.get("allow_missing_class_ii")),
        allow_missing_minimal_epitopes=bool(args.get("allow_missing_minimal_epitopes")),
        allow_ungrounded_safety_grade=bool(args.get("allow_ungrounded_safety_grade")))
    return {"content": [{"type": "text", "text": msg}]}


@tool("add_table",
      "Bulk-add entities to a section by mapping xlsx columns to schema fields - reads ALL rows in "
      "one deterministic call. PREFER this over add_entities for table-derived sections (peptides, "
      "epitopes). A column may be addressed by NAME or by 0-based POSITION (position reaches "
      "merged/blank/duplicated-header columns a name can't). "
      "mapping_json = {\"filter\"?: {\"col\":H|\"col_idx\":N|\"col_letter\":\"L\", \"in\"|\"equals\"|\"not_empty\":...}, "
      "\"fields\": {field: {\"col\":H | \"col_idx\":N | \"col_letter\":\"L\" | \"const\":v | "
      "\"template\":\"..{Col}..{#N}..{@L}..\" | \"template_list\":\"..\"}}} "
      "(template tokens: {Header} by name, {#N} by 0-based index, {@L} by Excel letter). "
      "Any field rule may add \"extract\":\"regex(group)\" to post-process its value (e.g. strip a prefix). "
      "MULTI-SHEET: pass \"sheets\":[name,...] to apply the SAME mapping across many per-entity sheets in "
      "ONE call (the cheap way to load papers that split data into one sheet per patient, e.g. ~34 "
      "'IAP-<patient>' immunogenicity tabs). Each row then carries a reserved \"__sheet__\" column "
      "(the sheet name), so derive per-sheet fields from it, e.g. "
      "{\"patient_paper_id\":{\"col\":\"__sheet__\",\"extract\":\"IAP-(.+)\"}}. "
      "Atomic: if any generated row is invalid, nothing is added and the bad rows are reported. "
      "Omit patient_paper_id when bulk-adding peptides to avoid per-patient count reconciliation; "
      "link patients via evidence.",
      {"type": "object",
       "properties": {"out_path": {"type": "string"},
                      "section": {"type": "string"},
                      "path": {"type": "string", "description": "path to the .xlsx file"},
                      "mapping_json": {"type": "string", "description": "the column->field mapping (see description)"},
                      "sheet": {"type": "string", "description": "single sheet name; omit for the first sheet"},
                      "sheets": {"type": "array", "items": {"type": "string"},
                                 "description": "MULTI-SHEET: list of sheet names; the same mapping is applied to each (each row gets a reserved __sheet__ column). Overrides `sheet`."},
                      "header_row": {"type": "integer",
                                     "description": "0-based header row; omit to auto-skip leading title rows"}},
       "required": ["out_path", "section", "path", "mapping_json"]})
async def add_table(args):
    sheets = args.get("sheets") or None
    _, msg = agent_core.table_to_entities(
        args["out_path"], args["section"], args["path"], args["mapping_json"], args.get("sheet") or None,
        header_row=args.get("header_row") if args.get("header_row") is not None else None,
        sheets=sheets)
    # OBSERVABILITY (2026-06-09): make the bulk-vs-per-sheet choice auditable in the per-paper log.
    # mode=multi(N) => the agent used add_table sheets=[...] across N tabs (the cheap path for
    # per-patient-sheet papers like 33064988); mode=single => one sheet per call. grep '[add_table]'.
    mode = f"multi({len(sheets)})" if sheets else "single"
    print(f"[add_table] section={args['section']} mode={mode} -> {msg[:240]}")
    return {"content": [{"type": "text", "text": msg}]}


@tool("build_pools",
      "DETERMINISTICALLY build one per-patient ExtractedPeptidePool from a per-patient peptide-ASSIGNMENT "
      "sheet (one row per patient x peptide, e.g. 33064988 'Vaccine peptides': 'Patient Alias' | "
      "'Peptide Sequence'). It groups the sheet by patient_col and, for each patient, sets "
      "member_peptide_ids to the ALREADY-LOADED immunizing_peptides matched by sequence. Call this ONCE "
      "(after loading the peptides with add_table) INSTEAD of hand-building pools with add_entities -- the "
      "group-by is identical every run, so it removes the per-patient pool variance. patient_col/peptide_col "
      "may be a header NAME or 0-based index.",
      {"type": "object",
       "properties": {"out_path": {"type": "string"},
                      "path": {"type": "string", "description": "path to the .xlsx file with the assignment table"},
                      "patient_col": {"type": ["string", "integer"], "description": "patient column (header name or 0-based index)"},
                      "peptide_col": {"type": ["string", "integer"], "description": "peptide-sequence column (header name or 0-based index)"},
                      "sheet": {"type": "string", "description": "sheet name; omit for the first sheet"},
                      "section_ref": {"type": "string", "description": "provenance section_ref for the pools (e.g. 'Supp Table 8; Figure 4A')"},
                      "quoted_text_template": {"type": "string", "description": "optional; may use {patient} and {n} (member count)"},
                      "header_row": {"type": "integer", "description": "0-based header row; omit to auto-skip title rows"}},
       "required": ["out_path", "path", "patient_col", "peptide_col"]})
async def build_pools(args):
    ok, msg = agent_core.build_patient_pools(
        args["out_path"], args["path"], args["patient_col"], args["peptide_col"],
        sheet=args.get("sheet") or None, section_ref=args.get("section_ref") or "",
        quoted_text_template=args.get("quoted_text_template") or None,
        header_row=args.get("header_row") if args.get("header_row") is not None else None)
    print(f"[build_pools] -> {msg[:240]}")
    return {"content": [{"type": "text", "text": msg}]}


@tool("build_pool_evidence",
      "DETERMINISTICALLY emit ONE pool-target immunogenic evidence row per MONITORED patient. Use for a "
      "paper whose per-patient pool immunogenicity is stated uniformly in TEXT ('de novo responses in all "
      "patients') but quantified only in per-patient FIGURES with no backing table (e.g. 33064988 Fig 4A / "
      "Supp Fig 3) -- the agent kept building these by hand and the count swung 0<->N run-to-run. Call this "
      "ONCE, AFTER build_pools. The MONITORED set = the patients the paper INDIVIDUALLY shows (NOT every "
      "vaccinated patient): pass patients=[...] OR sheets_path+sheet_pattern to derive them from the "
      "per-patient tab names (e.g. sheet_pattern='IAP-(.+)'). Each row references the patient's existing "
      "pool, magnitude=null + needs_review=true (figure-magnitude backfilled out-of-band), provenance "
      "kind='figure'. Faithfulness: it will NOT invent a row for a monitored patient that has no pool.",
      {"type": "object",
       "properties": {"out_path": {"type": "string"},
                      "patients": {"type": "array", "items": {"type": "string"},
                                   "description": "explicit monitored-patient ids; OR use sheets_path+sheet_pattern"},
                      "sheets_path": {"type": "string", "description": "xlsx whose per-patient tab names enumerate the monitored set"},
                      "sheet_pattern": {"type": "string", "description": "regex with one group over sheet names, e.g. 'IAP-(.+)'"},
                      "assay": {"type": "string", "description": "assay (default 'elispot')"},
                      "outcome": {"type": "string", "description": "outcome (default 'immunogenic')"},
                      "section_ref": {"type": "string", "description": "e.g. 'Figure 4A; Supplemental Figure 3'"},
                      "provenance_locator": {"type": "string", "description": "figure locator (default same as section_ref)"},
                      "quoted_text_template": {"type": "string", "description": "optional; may use {patient}"}},
       "required": ["out_path"]})
async def build_pool_evidence(args):
    ok, msg = agent_core.build_pool_evidence(
        args["out_path"], patients=args.get("patients") or None,
        sheets_path=args.get("sheets_path") or None, sheet_pattern=args.get("sheet_pattern") or None,
        assay=args.get("assay") or "elispot", outcome=args.get("outcome") or "immunogenic",
        section_ref=args.get("section_ref") or "Figure 4A; Supplemental Figure 3",
        provenance_locator=args.get("provenance_locator") or args.get("section_ref") or "Figure 4A; Supplemental Figure 3",
        quoted_text_template=args.get("quoted_text_template") or None)
    print(f"[build_pool_evidence] -> {msg[:240]}")
    return {"content": [{"type": "text", "text": msg}]}


@tool("build_crossreactivity_evidence",
      "DETERMINISTICALLY load a mutant-vs-WT cross-reactivity TABLE (33064988 Supplemental Table 5, in "
      "mmc1.pdf: 'Peptide ID | Mutant seq | WT seq | Cross reactive to WT') into one MinimalEpitope + one "
      "epitope-target immunogenic evidence row per listed peptide. The table is READABLE text, but the "
      "agent kept hand-building these via add_entities and dropping them run-to-run. Call this ONCE, AFTER "
      "the immunizing peptides are loaded (it links each epitope to its parent peptide by sequence "
      "containment). mutation_specific is set from the 'cross reactive to WT' column (No => mutant-specific); "
      "mhc_class is inferred from length (<=11 -> I else II) and the epitope is flagged needs_review. A row "
      "whose parent peptide isn't loaded is skipped (no orphan invented).",
      {"type": "object",
       "properties": {"out_path": {"type": "string"},
                      "pdf_path": {"type": "string", "description": "PDF holding the cross-reactivity table (e.g. mmc1.pdf)"},
                      "section_ref": {"type": "string", "description": "default 'Supplemental Table 5; Figure 4B'"},
                      "assay": {"type": "string", "description": "assay (default 'elispot')"},
                      "provenance_locator": {"type": "string", "description": "default 'Supplemental Table 5'"}},
       "required": ["out_path", "pdf_path"]})
async def build_crossreactivity_evidence(args):
    ok, msg = agent_core.build_crossreactivity_evidence(
        args["out_path"], args["pdf_path"],
        section_ref=args.get("section_ref") or "Supplemental Table 5; Figure 4B",
        assay=args.get("assay") or "elispot",
        provenance_locator=args.get("provenance_locator") or "Supplemental Table 5")
    print(f"[build_crossreactivity_evidence] -> {msg[:240]}")
    return {"content": [{"type": "text", "text": msg}]}


server = create_sdk_mcp_server(
    name="antvac",
    version="1.0.0",
    tools=[read_table, read_docx, read_pdf_text, survey_sources, read_figure, validate, save_extraction,
           init_record, add_entities, set_safety_summary, clear_entities, partial_status, finalize, add_table,
           build_pools, build_pool_evidence, build_crossreactivity_evidence],
)

SYSTEM_PROMPT = build_system_prompt(
    field_guidance_only(agent_core.DELTAS_TEXT),
    build_schema_digest(agent_core.schema, agent_core.vocab),
    agent_core.vocab,
)

_ANTVAC_TOOLS = ["mcp__antvac__read_table", "mcp__antvac__read_pdf_text",
                 "mcp__antvac__survey_sources",
                 "mcp__antvac__read_figure", "mcp__antvac__validate",
                 "mcp__antvac__save_extraction", "mcp__antvac__init_record",
                 "mcp__antvac__add_entities", "mcp__antvac__set_safety_summary",
                 "mcp__antvac__clear_entities",
                 "mcp__antvac__partial_status", "mcp__antvac__finalize",
                 "mcp__antvac__add_table", "mcp__antvac__build_pools",
                 "mcp__antvac__build_pool_evidence", "mcp__antvac__build_crossreactivity_evidence"]

# can_use_tool is NOT consulted for read-only host tools (Grep/Read/Glob) in default
# permission mode, so allowed_tools alone let the agent grep/read the SDK's *spilled*
# read_table result files to slice the big neoantigen table (run 4: 31 Greps). A hard
# DENYLIST does block them -- forcing the agent onto read_table(row_filter=...) / add_table,
# which read the xlsx deterministically. Bash/Write/Edit stay blocked for safety too.
# NOTE (2026-06-04): `tools=[]` now stops these from being LOADED at all (see options),
# so they no longer cost context tokens; the denylist + can_use_tool remain as
# defense-in-depth (deny anything that somehow appears).
_DENY_HOST_TOOLS = ["Bash", "Read", "Write", "Edit", "MultiEdit", "NotebookEdit",
                    "Grep", "Glob", "Task", "WebFetch", "WebSearch"]


async def _only_antvac(tool_name, tool_input, context):
    """Hard-restrict the agent to the antvac MCP tools.

    The SDK otherwise exposes the full host toolset (Bash/Read/Write/...). Denying
    everything else keeps the agent inside the validate→save loop and headless-safe
    (this callback IS the permission decision, so nothing blocks waiting on a human).
    """
    if tool_name.startswith("mcp__antvac__"):
        return PermissionResultAllow()
    return PermissionResultDeny(
        message=(f"Tool '{tool_name}' is not available. Use ONLY the antvac tools: "
                 "read_table, read_docx, read_pdf_text, read_figure, init_record, add_entities, "
                 "add_table, clear_entities, partial_status, finalize."))


async def extract_paper(paper_dir: str, out_path: str) -> None:
    files = [str(p) for p in pathlib.Path(paper_dir).rglob("*")
             if p.suffix.lower() in (".pdf", ".xlsx", ".docx")]
    if not files:
        print(f"[abort] no .pdf/.xlsx/.docx files found in {paper_dir}")
        sys.exit(1)
    options = ClaudeAgentOptions(
        model=MODEL,
        system_prompt=SYSTEM_PROMPT,
        mcp_servers={"antvac": server},
        # CONTEXT DIET (2026-06-04; count updated 2026-06-17): only the 16 antvac tools
        # in _ANTVAC_TOOLS should occupy context.
        # `tools=[]` disables ALL built-in tools (Bash/Read/Task/Cron*/Skill/AskUserQuestion/
        # ToolSearch/...; 22 of them) so their schemas are never loaded. strict_mcp_config
        # ignores every MCP server except the ones in mcp_servers, dropping the ~45 claude.ai
        # account servers (Vercel/PubMed/HuggingFace/Synthesize/Google/...). Net: 78 -> 16
        # tool schemas. cacheRead = turns x per-turn-context, so this trims the dominant cost
        # for every paper (and removes distractor tools). antvac arrives via mcp_servers, not
        # the built-in set, so tools=[] does not touch it.
        tools=[],
        strict_mcp_config=True,
        allowed_tools=_ANTVAC_TOOLS,
        disallowed_tools=_DENY_HOST_TOOLS,  # defense-in-depth (tools=[] already unloads these)
        can_use_tool=_only_antvac,    # belt-and-suspenders: deny every non-antvac tool
        setting_sources=[],           # do NOT load host user/project/local settings: no hooks/skills/CLAUDE.md
        permission_mode="default",    # so can_use_tool is consulted (bypassPermissions would skip it)
        max_turns=MAX_TURNS,
        max_budget_usd=MAX_BUDGET_USD,  # runaway-cost backstop
        # read_figure returns a base64 page image; one such message can exceed the SDK's
        # default 1 MB stdout NDJSON buffer and hard-crash the transport (run 9). 32 MB
        # leaves generous headroom for a full-page PNG.
        max_buffer_size=32 * 1024 * 1024,
    )
    prompt = (f"Extract this paper into {out_path}. Files available:\n"
              + "\n".join(files) +
              "\nStart by reading the neoantigen/ELISpot supplementary tables.")
    result = None
    async with ClaudeSDKClient(options=options) as client:
        await client.query(prompt)
        async for msg in client.receive_response():
            print(msg)  # stream tool calls / text to the log
            if type(msg).__name__ == "ResultMessage":
                result = msg

    if result is not None:
        tag = "[warn] agent loop ENDED IN ERROR" if getattr(result, "is_error", False) else "[info] agent loop ok"
        print(f"{tag}: subtype={getattr(result, 'subtype', '?')} "
              f"turns={getattr(result, 'num_turns', '?')} "
              f"cost=${(getattr(result, 'total_cost_usd', 0) or 0):.2f} "
              f"errors={getattr(result, 'errors', None)}")

    ok, gate_msg = agent_core.outer_guard(out_path)
    print(gate_msg)
    if not ok:
        sys.exit(1)


def _parse_cli(argv):
    """Parse argv (excluding the program name) → (paper_dir, out_path, subscription).

    `--subscription` is an optional, position-independent flag.
    """
    args = [a for a in argv if a != "--subscription"]
    subscription = "--subscription" in argv
    if not args:
        raise SystemExit(
            "usage: vaxtract [--subscription] <paper_dir> [out.json]")
    paper_dir = args[0]
    out = args[1] if len(args) > 1 else "newpaper_extracted.json"
    return paper_dir, out, subscription


def _apply_auth_mode(subscription):
    """Decide how the spawned `claude` CLI authenticates and return a label.

    The Agent SDK runs `claude` as a subprocess (it inherits this process's env).
    With `--subscription` we drop ANTHROPIC_API_KEY so the CLI falls back to the
    user's Claude subscription (plan quota) instead of API pay-per-token billing.
    """
    if subscription:
        os.environ.pop("ANTHROPIC_API_KEY", None)
        return "subscription (ANTHROPIC_API_KEY unset; claude CLI uses your plan quota)"
    return ("api-key (ANTHROPIC_API_KEY present)"
            if os.environ.get("ANTHROPIC_API_KEY")
            else "subscription (no ANTHROPIC_API_KEY in env)")


