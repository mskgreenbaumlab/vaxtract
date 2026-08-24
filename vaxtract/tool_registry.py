"""Provider-neutral tool catalog for the extraction agent.

Handlers call agent_core (SDK-free). Backends adapt ToolSpec to MCP or
Chat Completions function tools. Adding a finalize override: put it on
``finalize_partial`` — the schema and kwargs are derived from that signature.
"""
from __future__ import annotations

import inspect
from typing import Any

from . import agent_core
from .backends.base import ToolResult, ToolSpec

MCP_PREFIX = "mcp__antvac__"
DENY_HOST_TOOLS = ["Bash", "Read", "Write", "Edit", "MultiEdit", "NotebookEdit",
                   "Grep", "Glob", "Task", "WebFetch", "WebSearch"]


def mcp_tool_name(name: str) -> str:
    return f"{MCP_PREFIX}{name}"


def _ok(text: str) -> ToolResult:
    return ToolResult(text=text)


def _err(text: str) -> ToolResult:
    return ToolResult(text=text, is_error=True)


def _schema(properties: dict, required: list[str] | None = None) -> dict:
    spec: dict[str, Any] = {"type": "object", "properties": properties}
    if required:
        spec["required"] = required
    return spec


# ---------------------------------------------------------------------------
# handlers
# ---------------------------------------------------------------------------

async def handle_read_table(args: dict) -> ToolResult:
    try:
        text = agent_core.read_table_rows(
            args["path"], args.get("sheet") or None, int(args.get("max_rows") or 500),
            row_filter=args.get("row_filter") or None, columns=args.get("columns") or None,
            underline=bool(args.get("underline")),
            header_row=args.get("header_row") if args.get("header_row") is not None else None)
    except Exception as e:
        return _err(f"ERROR reading table {args.get('path')!r}: {e}")
    return _ok(text)


async def handle_read_docx(args: dict) -> ToolResult:
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
        return _err(f"ERROR reading docx {args.get('path')!r}: {e}")
    return _ok(text)


async def handle_read_pdf_text(args: dict) -> ToolResult:
    try:
        text = agent_core.read_pdf_text_from(
            args["path"], int(args.get("max_chars") or 40_000), int(args.get("offset") or 0))
    except Exception as e:
        return _err(f"ERROR reading pdf {args.get('path')!r}: {e}")
    return _ok(text)


async def handle_survey_sources(args: dict) -> ToolResult:
    try:
        text = agent_core.survey_sources(args["path"], int(args.get("max_chars") or 14000))
    except Exception as e:
        return _err(f"ERROR surveying {args.get('path')!r}: {e}")
    return _ok(text)


async def handle_read_figure(args: dict) -> ToolResult:
    try:
        region = args.get("region")
        region = tuple(region) if region else None
        b64, w, h = agent_core.render_figure_image(args["path"], args["page"], region=region)
    except Exception as e:
        return _err(f"ERROR rendering figure {args.get('path')!r} p{args.get('page')}: {e}")
    note = (f"Rendered {'region ' + str(region) if region else 'full page'} "
            f"of {args['path']} page {args['page']} ({w}x{h}). Read the values for "
            f"{args.get('what')!r}. If this is the full page, locate the panel and call "
            f"read_figure again with region=[x0,y0,x1,y1] for a legible crop. RECORD: "
            f"value=null (or a number only for a clean simple chart) + estimate in raw, "
            f"tier='reported', confidence<=2, Provenance(kind='figure', needs_review=true), "
            f"quoted_text = a verbatim figure-caption fragment.")
    return ToolResult(text=note, image_png_b64=b64)


async def handle_validate(args: dict) -> ToolResult:
    _, msg = agent_core.validate_record(args["candidate_json"])
    return _ok(msg)


async def handle_save_extraction(args: dict) -> ToolResult:
    _, msg = agent_core.save_record(args["candidate_json"], args["out_path"])
    return _ok(msg)


async def handle_init_record(args: dict) -> ToolResult:
    _, msg = agent_core.init_partial(args["out_path"], args["paper_meta_json"])
    return _ok(msg)


async def handle_add_entities(args: dict) -> ToolResult:
    _, msg = agent_core.append_section(args["out_path"], args["section"], args["items_json"])
    return _ok(msg)


async def handle_set_safety_summary(args: dict) -> ToolResult:
    _, msg = agent_core.set_safety_summary(args["out_path"], args["safety_json"])
    return _ok(msg)


async def handle_clear_entities(args: dict) -> ToolResult:
    _, msg = agent_core.clear_section(args["out_path"], args["section"])
    return _ok(msg)


async def handle_partial_status(args: dict) -> ToolResult:
    _, msg = agent_core.partial_status(args["out_path"])
    return _ok(msg)


_FINALIZE_OVERRIDE_DOCS: dict[str, str] = {
    "allow_missing_magnitudes":
        "set true to proceed when some responses have no reported magnitude",
    "allow_missing_pools":
        "set true to proceed when a patient has pooled evidence but pool membership is unresolvable",
    "allow_member_level_pool_evidence":
        "set true to proceed when a pooled response is kept as per-member rows because the source deconvolutes every member",
    "allow_unknown_funnel_size":
        "set true to proceed when candidates exist but the paper reports no predicted total",
    "allow_candidate_bridge_mismatch":
        "set true to proceed when a candidate's selected_peptide_id IMP has a different sequence",
    "allow_regimen_divergence":
        "set true to proceed when patients of one arm have different delivery regimens (e.g. dose escalation)",
    "allow_evidence_count_mismatch":
        "set true to proceed when recorded evidence differs from the paper's stated immunogenic/negative counts because they are at a different grain (e.g. pooled)",
    "allow_peptide_count_mismatch":
        "set true to proceed when recorded peptides fall short of n_selected_reported because the paper's count includes peptides not individually listed",
    "allow_sparse_evidence":
        "set true to proceed when immune evidence covers only a subset of vaccinated patients because the trial immune-monitored only that subset",
    "allow_missing_class_ii":
        "proceed despite a quoted class-II restriction (DR/DP/DQ, I-A/I-E, 'class II-restricted') with no class-II record minted; CD4/helper alone is not this gap",
    "allow_missing_minimal_epitopes":
        "proceed despite a mislabeled long-peptide epitope or a dropped predicted minimal-epitope layer",
    "allow_ungrounded_safety_grade":
        "proceed when a grade>=3 treatment-related safety claim cannot be grounded in a verbatim grade>=3 raw quote (routes to needs_review)",
    "allow_missing_tcr_gateway":
        "proceed when tcr_seq_status asserts TCR-seq was done but no TcrSeqMethod/DataDeposition was extractable (set tcr_seq_status='mentioned_only'; routes to needs_review)",
    "allow_tcr_status_mismatch":
        "proceed when tcr_seq_status does not match the extracted TCR content (normally just set the field correctly; routes to needs_review)",
    "allow_missing_tcr_flags":
        "proceed when TCR content was extracted but no per-(neoantigen x patient) neoantigen_tcr_flags row could be written (routes to needs_review)",
    "allow_cd4_on_class_i_minimal":
        "proceed when a CD4 evidence row targets a class-I 8-11mer epitope instead of the parent immunizing peptide — only if the paper quotes that exact 9-mer as the CD4 target (routes to needs_review)",
    "allow_cd8_on_immunizing_peptide":
        "proceed when a CD8/EPT evidence row targets a long immunizing peptide instead of the nested class-I 8-11mer — only if the paper assayed the LONG peptide for CD8 (routes to needs_review)",
    "allow_class_ii_asp_dump":
        "proceed when many class-II 12-18mer epitopes are overlapping ASP windows inside immunizing peptides rather than named minimal class-II determinants (routes to needs_review)",
}


def finalize_override_names() -> list[str]:
    return [p for p in inspect.signature(agent_core.finalize_partial).parameters
            if p.startswith("allow_")]


def finalize_input_schema() -> dict:
    props: dict[str, Any] = {"out_path": {"type": "string"}}
    for name in finalize_override_names():
        props[name] = {
            "type": "boolean",
            "description": _FINALIZE_OVERRIDE_DOCS.get(name, f"override flag {name}"),
        }
    return _schema(props, ["out_path"])


async def handle_finalize(args: dict) -> ToolResult:
    kwargs = {name: bool(args.get(name)) for name in finalize_override_names()}
    _, msg = agent_core.finalize_partial(args["out_path"], **kwargs)
    return _ok(msg)


async def handle_add_table(args: dict) -> ToolResult:
    sheets = args.get("sheets") or None
    _, msg = agent_core.table_to_entities(
        args["out_path"], args["section"], args["path"], args["mapping_json"], args.get("sheet") or None,
        header_row=args.get("header_row") if args.get("header_row") is not None else None,
        sheets=sheets)
    mode = f"multi({len(sheets)})" if sheets else "single"
    print(f"[add_table] section={args['section']} mode={mode} -> {msg[:240]}")
    return _ok(msg)


async def handle_build_pools(args: dict) -> ToolResult:
    ok, msg = agent_core.build_patient_pools(
        args["out_path"], args["path"], args["patient_col"], args["peptide_col"],
        sheet=args.get("sheet") or None, section_ref=args.get("section_ref") or "",
        quoted_text_template=args.get("quoted_text_template") or None,
        header_row=args.get("header_row") if args.get("header_row") is not None else None)
    print(f"[build_pools] -> {msg[:240]}")
    return _ok(msg) if ok else _err(msg)


async def handle_build_pool_evidence(args: dict) -> ToolResult:
    ok, msg = agent_core.build_pool_evidence(
        args["out_path"], patients=args.get("patients") or None,
        sheets_path=args.get("sheets_path") or None, sheet_pattern=args.get("sheet_pattern") or None,
        assay=args.get("assay") or "elispot", outcome=args.get("outcome") or "immunogenic",
        section_ref=args.get("section_ref") or "Figure 4A; Supplemental Figure 3",
        provenance_locator=args.get("provenance_locator") or args.get("section_ref") or "Figure 4A; Supplemental Figure 3",
        quoted_text_template=args.get("quoted_text_template") or None)
    print(f"[build_pool_evidence] -> {msg[:240]}")
    return _ok(msg) if ok else _err(msg)


async def handle_build_crossreactivity_evidence(args: dict) -> ToolResult:
    ok, msg = agent_core.build_crossreactivity_evidence(
        args["out_path"], args["pdf_path"],
        section_ref=args.get("section_ref") or "Supplemental Table 5; Figure 4B",
        assay=args.get("assay") or "elispot",
        provenance_locator=args.get("provenance_locator") or "Supplemental Table 5")
    print(f"[build_crossreactivity_evidence] -> {msg[:240]}")
    return _ok(msg) if ok else _err(msg)


# ---------------------------------------------------------------------------
# catalog
# ---------------------------------------------------------------------------

_ROW_FILTER = {
    "type": "object",
    "description": (
        'keep only matching rows: {"col":H | "col_idx":N | "col_letter":"L", '
        '"in"|"equals"|"not_empty":...}; reports matched/total counts. Use '
        "col_idx/col_letter for blank/duplicated-header columns names can't reach."
    ),
}
_COLUMNS = {
    "type": "array",
    "items": {"type": ["string", "integer"]},
    "description": (
        "project to just these columns -- a header NAME (string) or a 0-based "
        "POSITION (integer) for nameless columns; output is keyed by name or [idx]"
    ),
}

TOOLS: list[ToolSpec] = [
    ToolSpec("read_table",
             "Parse an .xlsx supplementary/source-data file into rows (Layer 1, highest fidelity). "
             "Returns up to max_rows rows and reports the total, so you know if a table was truncated. "
             "To INSPECT or SLICE a big table (count/find rows by a column value), use row_filter and/or "
             "columns HERE -- do NOT try to grep or read the tool's spilled result file (host file tools "
             "are unavailable). A huge sheet is byte-capped (you still get the TRUE total) -- narrow with "
             "row_filter/columns/max_rows rather than chasing a spill. Leading 'Table Sn.' title rows are "
             "auto-skipped for column projection. Set underline=true to reveal minimal-epitope substrings "
             "marked by underlining (returned wrapped in <u>…</u>). To BULK-ADD rows to the record, use add_table.",
             _schema({
                 "path": {"type": "string", "description": "path to the .xlsx file"},
                 "sheet": {"type": "string", "description": "sheet name; omit to use the first sheet"},
                 "max_rows": {"type": "integer",
                              "description": "max rows to return (default 500); raise it to read a larger table fully"},
                 "row_filter": _ROW_FILTER,
                 "columns": _COLUMNS,
                 "underline": {"type": "boolean",
                               "description": "reveal underlined sub-sequences (e.g. minimal epitopes marked by "
                                              "underlining inside a longer peptide) wrapped in <u>…</u>"},
                 "header_row": {"type": "integer",
                                "description": "0-based header row; omit to auto-skip leading title rows"},
             }, ["path"]),
             handle_read_table),
    ToolSpec("read_docx",
             "Parse a .docx supplement — common for BMC/Molecular-Cancer & Frontiers (Data_Sheet_*.docx) and "
             "often holding per-entity data the xlsx/pdf readers miss. THREE modes: (a) omit table_index -> a "
             "SUMMARY listing each table's rows×cols + its caption (the preceding 'Supplementary Table Sn…' "
             "paragraph); (b) table_index=N -> that table (0-based) as a capped preview, with the SAME "
             "row_filter/columns/underline/byte-cap behaviour as read_table -> bulk-load it with add_table just "
             "like an xlsx; (c) text_offset=0 -> the paragraph PROSE in a paged window (next-offset reported). "
             "Do NOT read a spilled result file (host file tools are unavailable).",
             _schema({
                 "path": {"type": "string", "description": "path to the .docx file"},
                 "table_index": {"type": "integer", "description": "0-based table to read; omit for the table summary"},
                 "max_rows": {"type": "integer", "description": "max rows to return (default 500)"},
                 "row_filter": _ROW_FILTER,
                 "columns": _COLUMNS,
                 "underline": {"type": "boolean", "description": "reveal underlined sub-sequences wrapped in <u>…</u>"},
                 "text_offset": {"type": "integer", "description": "read PROSE from this char offset instead of a table"},
                 "max_chars": {"type": "integer", "description": "prose window size (default 40000)"},
             }, ["path"]),
             handle_read_docx,
             allowed=False),
    ToolSpec("read_pdf_text",
             "Extract selectable text from the main or a supplementary PDF (methods, legends, results). "
             "Returns a paging window (default 40000 chars from offset 0); if more text remains the result "
             "reports the next offset -- call again with that offset to continue. Do NOT read a spilled "
             "result file (host file tools are unavailable).",
             _schema({
                 "path": {"type": "string"},
                 "offset": {"type": "integer", "description": "start char (default 0); use the next-offset the result reports"},
                 "max_chars": {"type": "integer", "description": "window size (default 40000)"},
             }, ["path"]),
             handle_read_pdf_text),
    ToolSpec("survey_sources",
             "INVENTORY every supplement in ONE call so you can LOCATE where data lives before opening files. "
             "Pass the paper directory (recursed) or a single file. Returns, per .xlsx/.pdf/.docx: xlsx -> each "
             "sheet's name, rows x cols, and HEADER ROW (so a peptide/immunogenicity table is recognizable by "
             "its columns even when the file/sheet name is unhelpful); pdf -> page count + first text "
             "(flags image-only figure PDFs); docx -> table count + caption + first table header. "
             "CALL THIS FIRST on a paper with many supplements (don't guess which file holds the manifest). "
             "Byte-capped; dropped files are listed (re-call on a subdir/file to see them).",
             _schema({
                 "path": {"type": "string", "description": "paper directory (recursed) or a single supplement file"},
                 "max_chars": {"type": "integer", "description": "byte cap for the digest (default 14000)"},
             }, ["path"]),
             handle_survey_sources),
    ToolSpec("read_figure",
             "Read numbers off a figure when they exist in NO table/source-data/text. TWO STEPS: "
             "(1) call with path+page (no region) to SEE the page image and locate the panel; "
             "(2) call again with region=[x0,y0,x1,y1] (fractions of the page, top-left origin) "
             "to get a zoomed, legible crop, then read the values. RECORD conservatively: put the "
             "read estimate in Measurement.raw, set value=null (or a number only for a clean simple "
             "chart), tier='reported', confidence<=2; attach Provenance(kind='figure', "
             "needs_review=true); set the row's quoted_text to a VERBATIM fragment of the figure "
             "CAPTION (get it via read_pdf_text) and section_ref to the panel (e.g. 'Figure 1G'). "
             "NEVER fabricate: if a value is unreadable, leave value=null.",
             _schema({
                 "path": {"type": "string"},
                 "page": {"type": "integer"},
                 "what": {"type": "string"},
                 "region": {"type": "array", "items": {"type": "number"},
                            "description": "optional [x0,y0,x1,y1] fractions of the page"},
             }, ["path", "page", "what"]),
             handle_read_figure),
    ToolSpec("validate",
             "(DEPRECATED — small records only; prefer init_record/add_entities/finalize) Check a candidate ExtractedPaper JSON against the schema. Returns 'VALID' or the errors to fix.",
             _schema({"candidate_json": {"type": "string"}}, ["candidate_json"]),
             handle_validate),
    ToolSpec("save_extraction",
             "(DEPRECATED — small records only; prefer init_record/add_entities/finalize) Validate AND write the final JSON in one shot. Fails if it exceeds the output limit on large papers.",
             _schema({"candidate_json": {"type": "string"}, "out_path": {"type": "string"}},
                     ["candidate_json", "out_path"]),
             handle_save_extraction),
    ToolSpec("init_record",
             "Start a new record from paper-level fields (pmid, journal, year, title, cohort_size, "
             "indication_summary, + optional doi/pmcid/nct_id/n_enrolled). Entity lists are added later "
             "with add_entities. Creates an on-disk partial that survives across turns.",
             _schema({"out_path": {"type": "string"},
                      "paper_meta_json": {"type": "string", "description": "JSON object of paper-level fields"}},
                     ["out_path", "paper_meta_json"]),
             handle_init_record),
    ToolSpec("add_entities",
             "Append a batch (<=~50) of entities to one section of the in-progress record. section is one "
             "of: patients, immunizing_peptides, epitopes, pools, evidence, survival_outcomes. Each item is "
             "validated against its schema; if any item is invalid the WHOLE batch is rejected and the "
             "partial is unchanged. Call repeatedly to add all rows (do NOT omit non-immunogenic peptides). "
             "For evidence, add one row per REPORTED response, not one per peptide. For bulk peptide tables "
             "prefer add_table and omit patient_paper_id.",
             _schema({"out_path": {"type": "string"},
                      "section": {"type": "string"},
                      "items_json": {"type": "string", "description": "JSON list of entity objects"}},
                     ["out_path", "section", "items_json"]),
             handle_add_entities),
    ToolSpec("set_safety_summary",
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
             _schema({"out_path": {"type": "string"},
                      "safety_json": {"type": "string", "description": "JSON object of SafetySummary fields"}},
                     ["out_path", "safety_json"]),
             handle_set_safety_summary),
    ToolSpec("clear_entities",
             "Reset one section to empty (to correct a mistake), then re-add it with add_entities.",
             _schema({"out_path": {"type": "string"}, "section": {"type": "string"}},
                     ["out_path", "section"]),
             handle_clear_entities),
    ToolSpec("partial_status",
             "Report current per-section counts + the paper metadata. Use to check progress, or to resume "
             "after a context summary (the on-disk partial is the source of truth, not your memory).",
             _schema({"out_path": {"type": "string"}}, ["out_path"]),
             handle_partial_status),
    ToolSpec("finalize",
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
             finalize_input_schema(),
             handle_finalize),
    ToolSpec("add_table",
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
             _schema({
                 "out_path": {"type": "string"},
                 "section": {"type": "string"},
                 "path": {"type": "string", "description": "path to the .xlsx file"},
                 "mapping_json": {"type": "string", "description": "the column->field mapping (see description)"},
                 "sheet": {"type": "string", "description": "single sheet name; omit for the first sheet"},
                 "sheets": {"type": "array", "items": {"type": "string"},
                            "description": "MULTI-SHEET: list of sheet names; the same mapping is applied to each (each row gets a reserved __sheet__ column). Overrides `sheet`."},
                 "header_row": {"type": "integer",
                                "description": "0-based header row; omit to auto-skip leading title rows"},
             }, ["out_path", "section", "path", "mapping_json"]),
             handle_add_table),
    ToolSpec("build_pools",
             "DETERMINISTICALLY build one per-patient ExtractedPeptidePool from a per-patient peptide-ASSIGNMENT "
             "sheet (one row per patient x peptide, e.g. 33064988 'Vaccine peptides': 'Patient Alias' | "
             "'Peptide Sequence'). It groups the sheet by patient_col and, for each patient, sets "
             "member_peptide_ids to the ALREADY-LOADED immunizing_peptides matched by sequence. Call this ONCE "
             "(after loading the peptides with add_table) INSTEAD of hand-building pools with add_entities -- the "
             "group-by is identical every run, so it removes the per-patient pool variance. patient_col/peptide_col "
             "may be a header NAME or 0-based index.",
             _schema({
                 "out_path": {"type": "string"},
                 "path": {"type": "string", "description": "path to the .xlsx file with the assignment table"},
                 "patient_col": {"type": ["string", "integer"], "description": "patient column (header name or 0-based index)"},
                 "peptide_col": {"type": ["string", "integer"], "description": "peptide-sequence column (header name or 0-based index)"},
                 "sheet": {"type": "string", "description": "sheet name; omit for the first sheet"},
                 "section_ref": {"type": "string", "description": "provenance section_ref for the pools (e.g. 'Supp Table 8; Figure 4A')"},
                 "quoted_text_template": {"type": "string", "description": "optional; may use {patient} and {n} (member count)"},
                 "header_row": {"type": "integer", "description": "0-based header row; omit to auto-skip title rows"},
             }, ["out_path", "path", "patient_col", "peptide_col"]),
             handle_build_pools),
    ToolSpec("build_pool_evidence",
             "DETERMINISTICALLY emit ONE pool-target immunogenic evidence row per MONITORED patient. Use for a "
             "paper whose per-patient pool immunogenicity is stated uniformly in TEXT ('de novo responses in all "
             "patients') but quantified only in per-patient FIGURES with no backing table (e.g. 33064988 Fig 4A / "
             "Supp Fig 3) -- the agent kept building these by hand and the count swung 0<->N run-to-run. Call this "
             "ONCE, AFTER build_pools. The MONITORED set = the patients the paper INDIVIDUALLY shows (NOT every "
             "vaccinated patient): pass patients=[...] OR sheets_path+sheet_pattern to derive them from the "
             "per-patient tab names (e.g. sheet_pattern='IAP-(.+)'). Each row references the patient's existing "
             "pool, magnitude=null + needs_review=true (figure-magnitude backfilled out-of-band), provenance "
             "kind='figure'. Faithfulness: it will NOT invent a row for a monitored patient that has no pool.",
             _schema({
                 "out_path": {"type": "string"},
                 "patients": {"type": "array", "items": {"type": "string"},
                              "description": "explicit monitored-patient ids; OR use sheets_path+sheet_pattern"},
                 "sheets_path": {"type": "string", "description": "xlsx whose per-patient tab names enumerate the monitored set"},
                 "sheet_pattern": {"type": "string", "description": "regex with one group over sheet names, e.g. 'IAP-(.+)'"},
                 "assay": {"type": "string", "description": "assay (default 'elispot')"},
                 "outcome": {"type": "string", "description": "outcome (default 'immunogenic')"},
                 "section_ref": {"type": "string", "description": "e.g. 'Figure 4A; Supplemental Figure 3'"},
                 "provenance_locator": {"type": "string", "description": "figure locator (default same as section_ref)"},
                 "quoted_text_template": {"type": "string", "description": "optional; may use {patient}"},
             }, ["out_path"]),
             handle_build_pool_evidence),
    ToolSpec("build_crossreactivity_evidence",
             "DETERMINISTICALLY load a mutant-vs-WT cross-reactivity TABLE (33064988 Supplemental Table 5, in "
             "mmc1.pdf: 'Peptide ID | Mutant seq | WT seq | Cross reactive to WT') into one MinimalEpitope + one "
             "epitope-target immunogenic evidence row per listed peptide. The table is READABLE text, but the "
             "agent kept hand-building these via add_entities and dropping them run-to-run. Call this ONCE, AFTER "
             "the immunizing peptides are loaded (it links each epitope to its parent peptide by sequence "
             "containment). mutation_specific is set from the 'cross reactive to WT' column (No => mutant-specific); "
             "mhc_class is inferred from length (<=11 -> I else II) and the epitope is flagged needs_review. A row "
             "whose parent peptide isn't loaded is skipped (no orphan invented).",
             _schema({
                 "out_path": {"type": "string"},
                 "pdf_path": {"type": "string", "description": "PDF holding the cross-reactivity table (e.g. mmc1.pdf)"},
                 "section_ref": {"type": "string", "description": "default 'Supplemental Table 5; Figure 4B'"},
                 "assay": {"type": "string", "description": "assay (default 'elispot')"},
                 "provenance_locator": {"type": "string", "description": "default 'Supplemental Table 5'"},
             }, ["out_path", "pdf_path"]),
             handle_build_crossreactivity_evidence),
]


def get_tool(name: str) -> ToolSpec:
    for t in TOOLS:
        if t.name == name:
            return t
    raise KeyError(name)


def allowed_tools() -> list[ToolSpec]:
    return [t for t in TOOLS if t.allowed]


def allowed_mcp_names() -> list[str]:
    return [mcp_tool_name(t.name) for t in allowed_tools()]


async def dispatch(name: str, args: dict, *, vision: bool = True) -> ToolResult:
    """Run a registry tool. Unknown names are refused (confinement)."""
    try:
        spec = get_tool(name)
    except KeyError:
        return _err(f"ERROR: tool {name!r} is not available. Use ONLY the antvac tools.")
    if not spec.allowed:
        return _err(f"ERROR: tool {name!r} is registered but not on the allowlist.")
    if spec.name == "read_figure" and not vision:
        return _err(
            "ERROR read_figure: this profile has vision=false (text-only model). "
            "Do not invent bar heights. Leave figure-only values as value=null and quarantine."
        )
    return await spec.handler(args or {})
