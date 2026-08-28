"""Deterministic loader for immunogenicity / reactivity-matrix sheets.

WHY: Hu 33479501's evidence lane swung because each run authored a fresh `per_cell`
mapping. Prompt text, a finalize-time column dump, and a survey_sources dump all
arrived too late or were overridden. The census already knows n_expected, assay_cols,
data_start, patient/imp/assay-peptide columns. This module consumes that census and
emits one evidence row per positive assay CELL.

GRAIN (docs/PLAN_2026-08-24_evidence_swing.md §5, signed 2026-08-25):
  - assay-peptide column present (Hu 4b) -> target_kind=epitope on the assay peptide,
    parent = immunizing peptide. Mint a class-II epitope only when none exists;
    mhc_class_inferred=True + needs_review=True. Do not mint class-I 8-11mers
    (those need a predicted_affinity slot); hang on a unique existing epitope or
    fall back to the immunizing peptide.
  - no assay-peptide column (Hu 4a, Supp11b) -> immunizing peptide, or a unique
    existing nested epitope.

Does not invent patients. A row whose patient token does not resolve to an
already-loaded ExtractedPatient is skipped and reported.
"""
from __future__ import annotations

import re
import zipfile

from . import tcr_source_census as _census

_PEPTIDE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{8,55}$")
_PATIENT_NUM = re.compile(r"^\d+$")
_ID_PATIENT = re.compile(r"^(\d+)\s*-")
_PATIENT_TOKEN = re.compile(r"^(?:Pt|P|Patient)?[\s._-]*(\d+)$", re.I)
_QT_ASSAY_PEP = re.compile(r"assay peptide ([ACDEFGHIKLMNPQRSTVWY]{8,55})", re.I)
_QT_IMP_PEP = re.compile(r"immunizing peptide ([ACDEFGHIKLMNPQRSTVWY]{8,55})", re.I)

# Binary ELISpot matrix cells are 1/0 calls, not spot counts. Never set value.
def table_call_magnitude(call: int) -> dict:
    return {"raw": "1" if call == 1 else "0", "unit": "unknown", "tier": "reported"}


def _call_bit(raw) -> int | None:
    n = _census._norm_call(raw)
    if n in _census._POSITIVE_CALL:
        return 1
    if n in _census._NEGATIVE_CALL:
        return 0
    return None


def canon_patient(p) -> str:
    s = str(p or "").strip()
    m = _PATIENT_TOKEN.match(s)
    return f"Pt{m.group(1)}" if m else s


def _cell(row: dict, col) -> str:
    if col is None:
        return ""
    v = row.get(col)
    if v is None:
        return ""
    if isinstance(v, float):
        return str(int(v)) if v.is_integer() else str(v)
    return str(v).strip()


def _sheet_rows(zf: zipfile.ZipFile, sheet_name: str) -> list[dict]:
    members = dict(_census._sheet_index(zf))
    member = members.get(sheet_name)
    if not member:
        return []
    rows = _census._probe_rows_full(zf, member)
    needed = {v[1] for r in rows for v in r.values() if isinstance(v, tuple)}
    strings = _census._resolve_shared_strings(zf, needed) if needed else {}
    out = []
    for r in rows:
        out.append({c: (strings.get(v[1], "") if isinstance(v, tuple) else v)
                    for c, v in r.items()})
    return out


def classify_assay_column(label: str) -> dict:
    """Map a (possibly merged) header label to the evidence discriminator triple."""
    lab = (label or "").lower()
    if any(t in lab for t in ("pre-vax", "pre vax", "prevax", "baseline", "week 0", "wk 0")):
        phase = "pre_vaccine"
    elif any(t in lab for t in ("year 3", "yr3", "yr 3", "y3", "3-4.5", "long-term",
                                 "long term", "memory", "month 44", "44 months",
                                 "47 months", "year 4")):
        phase = "memory"
    elif any(t in lab for t in ("week 16", "w16", "wk 16", "wk16", "post-vax",
                                 "post vax", "post_vaccine")):
        phase = "post_vaccine"
    else:
        phase = "unspecified"

    if "ex vivo" in lab or "ex-vivo" in lab or "exvivo" in lab:
        stim = "ex_vivo"
    elif any(t in lab for t in ("pre-stim", "pre stim", "prestim",
                                 "after pre-stimulation", "expanded", "in vitro")):
        stim = "in_vitro_expanded"
    else:
        stim = "unspecified"

    if "minigene" in lab:
        fmt = "minigene"
    elif any(t in lab for t in ("autologous", "tumor", "tumour")):
        fmt = "autologous_tumor"
    elif any(t in lab for t in ("puls", "peptide", "apc")):
        fmt = "peptide_pulsed"
    else:
        fmt = "unspecified"
    return {"timepoint_phase": phase, "assay_stimulation": stim,
            "assay_antigen_format": fmt, "timepoint_label": (label or "")[:64]}


def subset_from_sheet(sheet_name: str) -> str | None:
    hay = (sheet_name or "").lower()
    if "cd8" in hay:
        return "cd8"
    if "cd4" in hay:
        return "cd4"
    return None


def _match_patient(raw: str, rec_ids: set[str]) -> str | None:
    canon = canon_patient(raw)
    if canon in rec_ids:
        return canon
    m = re.search(r"(\d+)$", canon)
    if not m:
        return None
    hits = [p for p in rec_ids if re.search(rf"(?:P|Pt)?{m.group(1)}$", p or "")]
    return hits[0] if len(hits) == 1 else None


def _imps_for(rec: dict, patient: str | None) -> list[dict]:
    rows = rec.get("immunizing_peptides") or []
    if not patient:
        return list(rows)
    pref = [p for p in rows if p.get("patient_paper_id") == patient]
    return pref or list(rows)


def _match_seq(seq: str, rows: list[dict], patient: str | None = None) -> dict | None:
    if not seq:
        return None
    u = seq.upper()
    hits = [r for r in rows if (r.get("sequence") or "").upper() == u]
    if patient:
        pref = [r for r in hits if r.get("patient_paper_id") == patient]
        if len(pref) == 1:
            return pref[0]
        if len(pref) > 1:
            return None
    return hits[0] if len(hits) == 1 else None


def _unique_parent(asp: str, rec: dict, patient: str | None) -> dict | None:
    if not asp:
        return None
    cands = []
    for p in _imps_for(rec, patient):
        seq = p.get("sequence") or ""
        if asp in seq and asp != seq:
            cands.append(p)
    return cands[0] if len(cands) == 1 else None


def _header_labels(rows: list[dict], data_start: int) -> dict:
    return _census._combined_header_labels(rows, data_start, lambda v: v if v is not None else "")


def iter_scored_cells(spec, zf: zipfile.ZipFile, counts: dict | None = None):
    """Yield (row, col, patient_raw, imp_seq, asp_seq, id_token, call) for call in {1, 0}.

    n.d. / blank / anything outside the positive and negative call sets is not a row.
    `counts` (if given) is updated with n_positives / n_zeros / n_nd_skipped for every
    assay cell, including cells skipped later for missing patient.
    """
    if counts is None:
        counts = {}
    counts.setdefault("n_positives", 0)
    counts.setdefault("n_zeros", 0)
    counts.setdefault("n_nd_skipped", 0)
    rows = _sheet_rows(zf, spec.sheet)
    data = rows[spec.data_start:]
    sheet_peptides = {p for p in (_cell(r, spec.imp_col) for r in data) if p}
    carried = None
    for r in data:
        explicit = _cell(r, spec.patient_col)
        if _PATIENT_NUM.match(explicit) or _PATIENT_TOKEN.match(explicit):
            carried = explicit
        for col in spec.assay_cols:
            bit = _call_bit(_cell(r, col))
            if bit is None:
                counts["n_nd_skipped"] += 1
                continue
            if bit == 1:
                counts["n_positives"] += 1
            else:
                counts["n_zeros"] += 1
            pt = explicit if (_PATIENT_NUM.match(explicit) or _PATIENT_TOKEN.match(explicit)) else ""
            if not pt:
                m = _ID_PATIENT.match(_cell(r, spec.id_col))
                pt = m.group(1) if m else (carried or "")
            if not pt:
                continue
            imp = _cell(r, spec.imp_col)
            asp = _cell(r, spec.assay_pep_col)
            if not imp and asp:
                hits = {p for p in sheet_peptides if asp in p and asp != p}
                if len(hits) == 1:
                    imp = next(iter(hits))
            yield r, col, pt, imp, asp, _cell(r, spec.id_col), bit


def iter_positive_cells(spec, zf: zipfile.ZipFile):
    """Yield (row_dict, assay_col, patient_raw, imp_seq, asp_seq, id_token) per positive cell."""
    for row, col, pt, imp, asp, id_tok, call in iter_scored_cells(spec, zf):
        if call == 1:
            yield row, col, pt, imp, asp, id_tok


def build_items(rec: dict, paper_dir: str, sheet: str | None = None) -> tuple[list, list, dict]:
    """(evidence_items, epitope_items_to_mint, report).

    `report` counts loaded / skipped per sheet. Does not mutate `rec`.
    """
    census = _census.census_tcr_sheets(paper_dir, with_structure=True)
    specs = list(census.reactivity)
    if sheet:
        specs = [s for s in specs if s.sheet == sheet]
        if not specs:
            return [], [], {"error": f"no reactivity matrix named {sheet!r} under {paper_dir}"}
    rec_patients = {p.get("paper_local_id") for p in (rec.get("patients") or []) if p.get("paper_local_id")}
    existing_epitope_ids = {e.get("paper_local_id") for e in (rec.get("epitopes") or []) if e.get("paper_local_id")}
    existing_keys = {_ev_key(e) for e in (rec.get("evidence") or [])}

    evidence, minted, minted_ids = [], [], set(existing_epitope_ids)
    report = {"sheets": [], "n_cells": 0, "n_evidence": 0, "n_minted_epitopes": 0,
              "skipped_no_patient": 0, "skipped_no_target": 0, "skipped_dup": 0,
              "n_positives": 0, "n_zeros": 0, "n_nd_skipped": 0,
              "n_immunogenic": 0, "n_not_immunogenic": 0}

    by_file: dict = {}
    for spec in specs:
        by_file.setdefault(spec.path, []).append(spec)

    for path, group in by_file.items():
        with zipfile.ZipFile(path) as zf:
            for spec in group:
                rows = _sheet_rows(zf, spec.sheet)
                labels = _header_labels(rows, spec.data_start)
                subset = subset_from_sheet(spec.sheet)
                n_sheet = 0
                cell_counts = {"n_positives": 0, "n_zeros": 0, "n_nd_skipped": 0}
                for _row, col, pt_raw, imp_seq, asp_seq, id_tok, call in iter_scored_cells(
                        spec, zf, cell_counts):
                    report["n_cells"] += 1
                    patient = _match_patient(pt_raw, rec_patients)
                    if not patient:
                        report["skipped_no_patient"] += 1
                        continue
                    disc = classify_assay_column(labels.get(col, ""))
                    cue = subset.upper() if subset else "T cell"
                    quoted = (
                        f"{cue} reactivity by IFN-g ELISPOT: {patient}"
                        + (f", immunizing peptide {imp_seq}" if imp_seq else "")
                        + (f", assay peptide {asp_seq}" if asp_seq else "")
                        + f", scored {call} under {labels.get(col, f'column {col}')}"
                    )[:2000]
                    item = {
                        "patient_paper_id": patient,
                        "assay": "elispot",
                        "outcome": "immunogenic" if call == 1 else "not_immunogenic",
                        "magnitude": table_call_magnitude(call),
                        "vaccine_induced": disc["timepoint_phase"] != "pre_vaccine",
                        "section_ref": spec.sheet,
                        "quoted_text": quoted,
                        "assay_detail": (labels.get(col) or f"assay column {col}")[:200],
                        "timepoint_phase": disc["timepoint_phase"],
                        "timepoint_label": disc["timepoint_label"],
                        "assay_stimulation": disc["assay_stimulation"],
                        "assay_antigen_format": disc["assay_antigen_format"],
                        "provenance": [{"kind": "table", "locator": spec.sheet}],
                    }
                    if subset:
                        item["t_cell_subset"] = subset
                    target = _resolve_target(rec, patient, imp_seq, asp_seq, spec,
                                             minted, minted_ids, quoted, spec.sheet, id_tok)
                    if not target:
                        report["skipped_no_target"] += 1
                        continue
                    item.update(target)
                    key = _ev_key(item)
                    if key in existing_keys:
                        report["skipped_dup"] += 1
                        continue
                    existing_keys.add(key)
                    evidence.append(item)
                    n_sheet += 1
                    if call == 1:
                        report["n_immunogenic"] += 1
                    else:
                        report["n_not_immunogenic"] += 1
                report["n_positives"] += cell_counts["n_positives"]
                report["n_zeros"] += cell_counts["n_zeros"]
                report["n_nd_skipped"] += cell_counts["n_nd_skipped"]
                report["sheets"].append({"sheet": spec.sheet, "file": spec.file,
                                         "path": spec.path,
                                         "n_expected": spec.n_expected, "n_loaded": n_sheet,
                                         "n_positives": cell_counts["n_positives"],
                                         "n_zeros": cell_counts["n_zeros"],
                                         "n_nd_skipped": cell_counts["n_nd_skipped"]})
                report["n_evidence"] += n_sheet
    report["n_minted_epitopes"] = len(minted)
    return evidence, minted, report


def _resolve_target(rec, patient, imp_seq, asp_seq, spec, minted, minted_ids,
                    quoted, section_ref, id_tok) -> dict | None:
    asp = (asp_seq or "").upper()
    imp = (imp_seq or "").upper()
    if asp and not _PEPTIDE.match(asp):
        asp = ""
    if imp and not _PEPTIDE.match(imp):
        imp = ""

    parent = _match_seq(imp, rec.get("immunizing_peptides") or [], patient) if imp else None
    if parent is None and asp:
        parent = _unique_parent(asp, rec, patient)
    # Hu Supp11b's peptide column is a nested predicted 15-mer, not the administered
    # 25-mer. Fall back to unique containment of that sequence in a loaded IMP.
    if parent is None and imp:
        parent = _unique_parent(imp, rec, patient)

    if spec.assay_pep_col is not None:
        seq = asp
    elif parent is None:
        seq = imp          # predicted assay peptide sitting in the "immunizing" column
    else:
        seq = asp or ""

    if seq:
        epi = _match_seq(seq, rec.get("epitopes") or [], None)
        if epi is None:
            epi = next((e for e in minted if (e.get("sequence") or "").upper() == seq), None)
        if epi is not None:
            return {"target_kind": "epitope", "epitope_paper_id": epi["paper_local_id"],
                    "immunizing_peptide_paper_id": None}
        if 8 <= len(seq) <= 11:
            nested = [e for e in (rec.get("epitopes") or [])
                      if (e.get("sequence") or "").upper() == seq]
            if len(nested) == 1:
                return {"target_kind": "epitope", "epitope_paper_id": nested[0]["paper_local_id"],
                        "immunizing_peptide_paper_id": None}
            # Class-I 8-11mer assayed as itself (Hu 4a mutated peptide). The
            # schema requires a predicted_affinity slot on class-I epitopes;
            # the sheet may not print nM on this tab, so the slot is lossless.
            eid = _mint_id(id_tok, seq, minted_ids)
            minted.append({
                "paper_local_id": eid,
                "sequence": seq,
                "mhc_class": "I",
                "mhc_class_inferred": True,
                "is_neoantigen": True,
                "needs_review": True,
                "parent_peptide_ids": [parent["paper_local_id"]] if parent else [],
                "predicted_affinity": {
                    "unit": "unknown",
                    "raw": "affinity not printed on this reactivity sheet",
                    "tier": "predicted",
                },
                "quoted_text": quoted[:2000],
                "section_ref": section_ref,
            })
            minted_ids.add(eid)
            rec.setdefault("epitopes", []).append(minted[-1])
            return {"target_kind": "epitope", "epitope_paper_id": eid,
                    "immunizing_peptide_paper_id": None}
        if len(seq) >= 12:
            # Mint a class-II assay peptide. Empty parents is legal: the paper reported
            # the peptide without resolving a long immunizing parent (Hu Supp11b).
            eid = _mint_id(id_tok, seq, minted_ids)
            minted.append({
                "paper_local_id": eid,
                "sequence": seq,
                "mhc_class": "II",
                "mhc_class_inferred": True,
                "is_neoantigen": True,
                "needs_review": True,
                "parent_peptide_ids": [parent["paper_local_id"]] if parent else [],
                "quoted_text": quoted[:2000],
                "section_ref": section_ref,
            })
            minted_ids.add(eid)
            rec.setdefault("epitopes", []).append(minted[-1])
            return {"target_kind": "epitope", "epitope_paper_id": eid,
                    "immunizing_peptide_paper_id": None}

    if parent is not None:
        return {"target_kind": "immunizing_peptide",
                "immunizing_peptide_paper_id": parent["paper_local_id"],
                "epitope_paper_id": None}
    return None


def _mint_id(id_tok: str, seq: str, taken: set) -> str:
    raw = (id_tok or "").strip()
    if raw and re.match(r"^[\w.-]+$", raw):
        cand = raw if raw.upper().startswith(("ASP", "EPT", "EPI")) else f"ASP_{raw}"
        if cand not in taken:
            return cand
    base = f"EPI_{seq}"
    if base not in taken:
        return base
    i = 2
    while f"{base}_{i}" in taken:
        i += 1
    return f"{base}_{i}"


def _ev_key(e: dict) -> tuple:
    return (e.get("patient_paper_id"), e.get("target_kind"), e.get("immunizing_peptide_paper_id"),
            e.get("epitope_paper_id"), e.get("pool_paper_id"), e.get("candidate_paper_id"),
            e.get("assay"), e.get("outcome"), e.get("t_cell_subset"), e.get("mhc_class"),
            e.get("assay_stimulation"), e.get("timepoint_phase"),
            e.get("assay_antigen_format"), e.get("assay_detail"))


def _seq_index(rec: dict) -> dict[str, str]:
    out = {}
    for lane in ("epitopes", "immunizing_peptides"):
        for row in rec.get(lane) or []:
            pid = row.get("paper_local_id")
            seq = (row.get("sequence") or "").upper()
            if pid and seq:
                out[pid] = seq
    return out


def _evidence_seqs(ev: dict, seqs: dict[str, str]) -> tuple[set[str], set[str]]:
    """(assay/epitope sequences, immunizing-peptide sequences). Prefer the assay peptide."""
    epi: set[str] = set()
    imp: set[str] = set()
    eid = ev.get("epitope_paper_id")
    iid = ev.get("immunizing_peptide_paper_id")
    if eid and eid in seqs:
        epi.add(seqs[eid])
    if iid and iid in seqs:
        imp.add(seqs[iid])
    qt = ev.get("quoted_text") or ""
    m = _QT_ASSAY_PEP.search(qt)
    if m:
        epi.add(m.group(1).upper())
    m = _QT_IMP_PEP.search(qt)
    if m:
        imp.add(m.group(1).upper())
    return epi, imp


_SCORED_UNDER = re.compile(r"scored [01] under (.+)$")


def _label_matches(ev: dict, label: str) -> bool:
    lab = (label or "").strip()
    if not lab:
        return False
    detail = (ev.get("assay_detail") or "").strip()
    if detail and (detail == lab[:200] or detail.startswith(lab)):
        return True
    qt = ev.get("quoted_text") or ""
    m = _SCORED_UNDER.search(qt)
    under = (m.group(1).strip() if m else "")
    if under and (under == lab or lab.startswith(under) or under.startswith(lab)):
        return True
    return False


def is_reactivity_section(section_ref: str | None, sheet_names: set[str]) -> bool:
    s = section_ref or ""
    if s in sheet_names:
        return True
    low = s.lower()
    return any(tag in low for tag in ("4a", "4b", "11b", "dataset 4", "reactivity"))


def scored_cell_census(rec: dict, paper_dir: str) -> tuple[list[dict], dict]:
    cells = []
    census = _census.census_tcr_sheets(paper_dir, with_structure=True)
    rec_patients = {p.get("paper_local_id") for p in (rec.get("patients") or []) if p.get("paper_local_id")}
    totals = {"n_positives": 0, "n_zeros": 0, "n_nd_skipped": 0, "sheets": []}
    by_file: dict = {}
    for spec in census.reactivity:
        by_file.setdefault(spec.path, []).append(spec)
    for path, group in by_file.items():
        with zipfile.ZipFile(path) as zf:
            for spec in group:
                rows = _sheet_rows(zf, spec.sheet)
                labels = _header_labels(rows, spec.data_start)
                cell_counts = {"n_positives": 0, "n_zeros": 0, "n_nd_skipped": 0}
                for _row, col, pt_raw, imp_seq, asp_seq, id_tok, call in iter_scored_cells(
                        spec, zf, cell_counts):
                    patient = _match_patient(pt_raw, rec_patients)
                    cells.append({
                        "sheet": spec.sheet,
                        "patient": patient,
                        "patient_raw": pt_raw,
                        "imp": (imp_seq or "").upper(),
                        "asp": (asp_seq or "").upper(),
                        "label": labels.get(col) or f"assay column {col}",
                        "call": call,
                        "id_tok": id_tok,
                    })
                for k in ("n_positives", "n_zeros", "n_nd_skipped"):
                    totals[k] += cell_counts[k]
                totals["sheets"].append({"sheet": spec.sheet, **cell_counts,
                                         "n_expected": spec.n_expected})
    return cells, totals


def patch_table_magnitudes(rec: dict, paper_dir: str) -> dict:
    """Write magnitude from the matching scored cell onto existing reactivity evidence.

    Does not invent raw='1' because the row is immunogenic. Unmatched rows stay unchanged.
    """
    cells, census = scored_cell_census(rec, paper_dir)
    sheet_names = {s["sheet"] for s in census["sheets"]}
    seqs = _seq_index(rec)
    patched, unmatched, ambiguous = [], [], []
    used = [False] * len(cells)
    for i, ev in enumerate(rec.get("evidence") or []):
        if not is_reactivity_section(ev.get("section_ref"), sheet_names):
            continue
        epi_seqs, imp_seqs = _evidence_seqs(ev, seqs)
        hits = []
        for j, cell in enumerate(cells):
            if used[j] or cell["patient"] != ev.get("patient_paper_id"):
                continue
            if cell["sheet"] != ev.get("section_ref"):
                continue
            if not _label_matches(ev, cell["label"]):
                continue
            asp, imp = cell["asp"], cell["imp"]
            if epi_seqs:
                if asp not in epi_seqs and not (not asp and imp in epi_seqs):
                    continue
            elif imp_seqs:
                if imp not in imp_seqs and asp not in imp_seqs:
                    continue
            hits.append(j)
        if not hits:
            unmatched.append({"index": i, "patient": ev.get("patient_paper_id"),
                              "section_ref": ev.get("section_ref"),
                              "epitope_paper_id": ev.get("epitope_paper_id"),
                              "assay_detail": ev.get("assay_detail")})
            continue
        detail = (ev.get("assay_detail") or "").strip()
        exact = [j for j in hits if detail and detail == cells[j]["label"][:200]]
        if len(exact) == 1:
            hits = exact
        calls = {cells[j]["call"] for j in hits}
        if len(calls) != 1:
            ambiguous.append({"index": i, "calls": sorted(calls),
                              "section_ref": ev.get("section_ref")})
            continue
        call = next(iter(calls))
        ev["magnitude"] = table_call_magnitude(call)
        used[hits[0]] = True
        patched.append(i)
    return {
        **census,
        "n_patched": len(patched),
        "n_unmatched": len(unmatched),
        "n_ambiguous": len(ambiguous),
        "unmatched": unmatched,
        "ambiguous": ambiguous,
        "sheet_names": sorted(sheet_names),
    }


def add_zero_evidence(rec: dict, paper_dir: str) -> dict:
    """Append missing not_immunogenic rows from scored 0-cells. Does not delete positives."""
    items, minted, report = build_items(rec, paper_dir)
    zeros = [it for it in items if it.get("outcome") == "not_immunogenic"]
    rec.setdefault("evidence", []).extend(zeros)
    # build_items/_resolve_target already appended minted epitopes onto rec
    report["n_zeros_added"] = len(zeros)
    report["n_positives_new_ignored"] = sum(1 for it in items if it.get("outcome") == "immunogenic")
    report["n_minted_for_zeros"] = len(minted)
    return report
