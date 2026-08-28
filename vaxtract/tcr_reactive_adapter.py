"""PROTOTYPE — deterministic loader for paired single-cell "neoantigen-reactive T cell" TCR tables.

Motivation (stability lab, 2026-07-06): the TCR lane is the last non-deterministic one on Keskin — the
live agent extracts 4 clonotypes one run and 0 the next (the metric's 8/0/8), and it captures only the
main-text "dominant" clones, not the full sorted-reactive set. Those clones live in clean paired
single-cell supplementary tables (one row per CELL, with the paper's own `Clone`/`Clonotype` grouping),
so reading them is a PURE FUNCTION of the sheet — byte-identical every run, and complete.

Shape handled (Keskin 30568305 MOESM7 = Pt8, MOESM10 = Pt7; two-row header, data from row 3):
    Cell designation | [Expression classification] | TRAV | TRAJ | Alpha CDR3 DNA seq |
    Alpha CDR3 amino acid seq | TRBV | TRBJ | TRBC | Beta CDR3 DNA seq | Beta CDR3 amino acid seq |
    Clone | Clonotype
One CLONOTYPE = one unique (CDR3a, CDR3b) amino-acid pair (== the paper's own Clonotype column). Cells
sharing the pair are one clonally-expanded clonotype; `cell_count` is how many cells carried it.

SCOPE is a curation choice, exposed as `min_cells` (default 2 = clonally expanded, dropping single-cell
singletons): on Keskin the tables hold 439 unique clonotypes, 40 expanded (>=2 cells), 16 at >=3 cells,
and a handful of dominant clones the paper headlines (the 49-cell Pt7 CD8 clone = "predominant CD8 clone
specific to EPT12A"). NOT wired into the extractor — this is the prototype for sign-off on that scope.
"""
from __future__ import annotations

import pathlib
import re

import openpyxl

from . import schema

_JUNCTION = re.compile(r"^C[ACDEFGHIKLMNPQRSTVWY]{2,30}[FW]$")   # IMGT junction: C…F/W


def _col(header, *names):
    for i, c in enumerate(header):
        cl = str(c).strip().lower() if c is not None else ""
        if any(n in cl for n in names):
            return i
    return None


def _aa(v):
    if not isinstance(v, str):
        return None
    s = v.strip().upper()
    return s if re.match(r"^[ACDEFGHIKLMNPQRSTVWY]{4,}$", s) else None


def _cell(row, i):
    return row[i] if i is not None and i < len(row) else None


def _chain(locus, cdr3, v, j, c=None):
    if not cdr3:
        return None
    is_j = bool(_JUNCTION.match(cdr3))
    return schema.TcrChain(
        locus=locus,
        junction_aa=cdr3 if is_j else None,
        cdr3_aa=None if is_j else cdr3,
        cdr3_scheme="imgt_junction_aa" if is_j else "cdr3_aa_no_flank",
        v_call=str(v).strip() if v else None,
        j_call=str(j).strip() if j else None,
    )


def load_reactive_clonotypes(xlsx_path, sheet, *, patient_paper_id, antigen, section_ref,
                             min_cells=2, header_row=1):
    """PURE FUNCTION: one reactive single-cell sheet -> list[schema.TcrClonotype], one per unique
    (CDR3a, CDR3b) pair with cell_count >= min_cells. Deterministic (clonotypes sorted by the pair,
    stable ids). `antigen` names the sorted specificity (e.g. 'mutant SHANK2'); it is recorded in
    quoted_text. needs_review stays false unless schema-forced (inferred_association
    or grouped observations the paper did not assert)."""
    wb = openpyxl.load_workbook(pathlib.Path(xlsx_path), read_only=True, data_only=True)
    ws = wb[sheet]
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    wb.close()
    hdr = rows[header_row]
    ca, cb = _col(hdr, "alpha cdr3 amino"), _col(hdr, "beta cdr3 amino")
    trav, traj = _col(hdr, "trav"), _col(hdr, "traj")
    trbv, trbj = _col(hdr, "trbv"), _col(hdr, "trbj")

    groups: dict = {}
    for r in rows[header_row + 1:]:
        a, b = _aa(_cell(r, ca)), _aa(_cell(r, cb))
        if not a and not b:
            continue
        key = (a, b)
        g = groups.setdefault(key, {"n": 0, "trav": _cell(r, trav), "traj": _cell(r, traj),
                                    "trbv": _cell(r, trbv), "trbj": _cell(r, trbj)})
        g["n"] += 1

    out = []
    for idx, (key, g) in enumerate(sorted(groups.items(), key=lambda kv: (kv[0][0] or "", kv[0][1] or "")), 1):
        if g["n"] < min_cells:
            continue
        a, b = key
        chains = [c for c in (_chain("TRA", a, g["trav"], g["traj"]),
                              _chain("TRB", b, g["trbv"], g["trbj"])) if c]
        pairing = "paired_ab" if a and b else ("tra_only" if a else "trb_only")
        out.append(schema.TcrClonotype(
            paper_local_id=f"TCR_{patient_paper_id}_{re.sub(r'[^A-Za-z0-9]', '', antigen)}_{idx}",
            patient_paper_id=patient_paper_id,
            receptor_modality="single_cell_tcr_seq",
            chain_pairing=pairing,
            chains=chains,
            # the schema requires specificity_evidence='not_established' until specific_for_kind is linked
            # to a target entity — that linkage (by gene/antigen) happens when this is wired into finalize
            # with the record's epitopes/candidates in hand; the sorted antigen is named in quoted_text.
            specificity_evidence="not_established",
            needs_review=False,  # antigen-sorted table row; H3 is inferred_association, not this
            quoted_text=(f"{patient_paper_id} {antigen}-reactive single-cell clonotype "
                         f"(CDR3a={a}, CDR3b={b}; {g['n']} cell{'s' if g['n'] != 1 else ''})"),
            section_ref=section_ref,
        ))
    return out


# The Keskin reactive-clone sheets (neoantigen-reactive only; MOESM10 (a) 'Tumor-associated' is excluded).
KESKIN_SHEETS = [
    ("41586_2018_792_MOESM7_ESM.xlsx", "(a)Pt8 mutSHANK2-react Tcells", "Pt8", "mutant SHANK2",
     "Supplementary Table 7a (Pt8 mutSHANK2-reactive T cells)"),
    ("41586_2018_792_MOESM7_ESM.xlsx", "(b)Pt8 mutSVEP1-react Tcells", "Pt8", "mutant SVEP1",
     "Supplementary Table 7b (Pt8 mutSVEP1-reactive T cells)"),
    ("41586_2018_792_MOESM10_ESM.xlsx", "(b)Pt7 Neoantigen-reactive CD4+", "Pt7", "pool C neoantigen (CD4)",
     "Supplementary Table 10b (Pt7 pool-C-reactive CD4+ T cells)"),
    ("41586_2018_792_MOESM10_ESM.xlsx", "(c)Pt7 Neoantigen-reactive CD8", "Pt7", "pool C neoantigen (CD8)",
     "Supplementary Table 10c (Pt7 pool-C-reactive CD8+ T cells)"),
]


def load_keskin(supps_dir, min_cells=2):
    """Load all four Keskin neoantigen-reactive sheets -> list[schema.TcrClonotype] at `min_cells`."""
    supps = pathlib.Path(supps_dir)
    out = []
    for fn, sheet, pt, antigen, ref in KESKIN_SHEETS:
        out.extend(load_reactive_clonotypes(supps / fn, sheet, patient_paper_id=pt, antigen=antigen,
                                            section_ref=ref, min_cells=min_cells))
    return out
