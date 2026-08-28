#!/usr/bin/env python3
"""Curator HTML renderer for a validated ExtractedPaper JSON.

Ships in the *core* vaxtract package (pydantic only — no agent extra). Pure
function of the JSON: never re-reads the paper. Refuses invalid records.
Empty sections omit. Optional in the extraction workflow.

Usage:
    vaxtract-report EXTRACTED.json [OUT.html]
    vaxtract report EXTRACTED.json [OUT.html]
    python -m vaxtract.report EXTRACTED.json [OUT.html]
"""
from __future__ import annotations

import collections
import html
import json
import sys

from .schema import SCHEMA_VERSION, ExtractedPaper

DASH = "—"
TCR_CLONE_CAP = 40
# Lockstep with vaxtract.agent_core.SOFT_OVERRIDES. Anything else in
# finalize_overrides_used is HARD (scale lane → needs_review).
SOFT_OVERRIDES = frozenset({
    "allow_unknown_funnel_size",
    "allow_regimen_divergence",
})


def esc(s):
    return html.escape(str(s)) if s is not None else ""


def mag_str(m):
    """Lossless magnitude label. Never invent SFC/1e6 from a raw token."""
    if not m:
        return ""
    if m.get("value") is not None:
        unit = m.get("unit") or ""
        label = {
            "sfc_per_1e6": "SFC/10\u2076",
            "percent_of_parent": "%",
            "stimulation_index": "SI",
            "unknown": "",
        }.get(unit, unit.replace("_", " "))
        return f"{m['value']:g} {label}".strip()
    raw = m.get("raw")
    unit = (m.get("unit") or "").lower()
    if raw in ("1", "0") and unit in ("unknown", ""):
        return f"table call {raw}"
    if m.get("grade"):
        return str(m["grade"]).replace("_", " ")
    return str(raw) if raw else ""


def count_needs_review(obj):
    n = 0
    if isinstance(obj, dict):
        if obj.get("needs_review") is True:
            n += 1
        for v in obj.values():
            n += count_needs_review(v)
    elif isinstance(obj, list):
        for v in obj:
            n += count_needs_review(v)
    return n


def responders(ex):
    r = set()
    for p in ex["patients"]:
        if (p.get("n_peptides_immunogenic") or 0) > 0:
            r.add(p["paper_local_id"])
    for e in ex["evidence"]:
        if e.get("outcome") in ("immunogenic", "positive"):
            r.add(e["patient_paper_id"])
    return r


def chain_cdr3(cl):
    bits = []
    for ch in cl.get("chains") or []:
        aa = ch.get("junction_aa") or ch.get("cdr3_aa") or ch.get("raw") or ""
        if aa:
            bits.append(aa)
    return ", ".join(bits) if bits else DASH


def tcr_target(row):
    return (
        row.get("epitope_paper_id")
        or row.get("immunizing_peptide_paper_id")
        or row.get("candidate_paper_id")
        or DASH
    )


def load_record(src):
    paper = ExtractedPaper(**json.loads(open(src).read()))
    return paper.model_dump(), SCHEMA_VERSION


def build_html(ex, sv):
    imp_by_id = {i["paper_local_id"]: i for i in ex["immunizing_peptides"]}
    epi_by_id = {e["paper_local_id"]: e for e in ex["epitopes"]}
    pool_by_id = {p["paper_local_id"]: p for p in ex["pools"]}
    cand_by_id = {c["paper_local_id"]: c for c in (ex.get("candidates") or [])}
    epi_by_parent = collections.defaultdict(dict)
    for e in ex["epitopes"]:
        for pid in e.get("parent_peptide_ids") or []:
            epi_by_parent[pid][e["mhc_class"]] = e

    n1 = sum(1 for e in ex["epitopes"] if e["mhc_class"] == "I")
    n2 = sum(1 for e in ex["epitopes"] if e["mhc_class"] == "II")
    n_imp_ev = sum(1 for e in ex["evidence"] if e.get("target_kind") == "immunizing_peptide")
    nmag = sum(1 for e in ex["evidence"] if e.get("magnitude"))
    n_tcr = len(ex.get("tcr_clonotypes") or [])
    n_nr = count_needs_review(ex)
    nresp = len(responders(ex))
    ov = list(ex.get("finalize_overrides_used") or [])
    hard_ov = [o for o in ov if o not in SOFT_OVERRIDES]
    empty_lanes = (not ex["immunizing_peptides"]) and (not ex["evidence"]) and bool(ex["patients"])

    def target_descr(e):
        kind = e.get("target_kind")
        if kind == "immunizing_peptide":
            i = imp_by_id.get(e.get("immunizing_peptide_paper_id"), {})
            key = i.get("paper_local_id")
        elif kind == "epitope":
            i = epi_by_id.get(e.get("epitope_paper_id"), {})
            key = (i.get("parent_peptide_ids") or [None])[0]
        elif kind == "candidate":
            c = cand_by_id.get(e.get("candidate_paper_id"), {})
            return (c.get("gene_symbol") or DASH, c.get("mutation") or DASH, "", "")
        elif kind == "pool":
            p = pool_by_id.get(e.get("pool_paper_id"), {})
            return (f"pool ({len(p.get('member_peptide_ids') or [])} peptides)", "", "", "")
        else:
            return (DASH, DASH, "", "")
        e1 = epi_by_parent.get(key, {}).get("I", {})
        e2 = epi_by_parent.get(key, {}).get("II", {})
        return (
            i.get("gene_symbol") or DASH,
            i.get("mutation") or DASH,
            e1.get("sequence") or i.get("sequence") or "",
            e2.get("sequence") or "",
        )

    def sec_health():
        dirty = bool(ov or empty_lanes or n_nr or ex.get("companion_paper_ref")
                     or ex.get("peptide_manifest_present") is False)
        if not dirty:
            return ""
        rows = ""
        def row(k, v):
            nonlocal rows
            rows += f"<tr><td class='mono'>{esc(k)}</td><td>{v}</td></tr>"
        if ov:
            chips = []
            for o in ov:
                tag = "HARD" if o not in SOFT_OVERRIDES else "SOFT"
                chips.append(f"<span class='b {'non' if tag=='HARD' else 'soft'}'>{esc(tag)} {esc(o)}</span>")
            row("finalize overrides", " ".join(chips))
        if ex.get("data_resolution"):
            row("data resolution", esc(ex["data_resolution"]))
        if ex.get("peptide_manifest_present") is not None:
            row("peptide manifest present", "yes" if ex["peptide_manifest_present"] else "no")
        if ex.get("tcr_seq_status"):
            row("TCR seq status", esc(ex["tcr_seq_status"]))
        row("needs_review", str(n_nr))
        if ex.get("companion_paper_ref"):
            row("companion paper", esc(ex["companion_paper_ref"]))
        if empty_lanes:
            man = ex.get("peptide_manifest_present")
            note = ("Peptide and evidence lanes are empty. "
                    + ("A sequence manifest was present in the source; sequences are not in this record."
                       if man is True else
                       "No sequence-bearing peptide entities were written."))
            row("empty lanes", esc(note))
        return ("<section><h2>Record health</h2>"
                "<table><thead><tr><th>field</th><th>value</th></tr></thead>"
                f"<tbody>{rows}</tbody></table></section>")

    def sec_safety():
        s = ex.get("safety_summary")
        if not s:
            return ""
        g3 = {True: "yes", False: "no"}.get(s.get("any_grade3plus_related"), DASH)
        irae = {True: "yes", False: "no"}.get(s.get("irae_present"), DASH)
        raw = f"<p class='dim'>{esc(s.get('raw'))}</p>" if s.get("raw") else ""
        return ("<section><h2>Safety</h2>"
                "<table><thead><tr><th>max related grade</th><th>grade ≥3 related</th>"
                "<th class='num'>n with related AE</th><th>irAE</th></tr></thead><tbody>"
                f"<tr><td class='num'>{s.get('max_related_grade', DASH)}</td>"
                f"<td>{g3}</td><td class='num'>{s.get('n_patients_with_related_ae', DASH)}</td>"
                f"<td>{irae}</td></tr></tbody></table>{raw}</section>")

    def sec_survival():
        so = ex.get("survival_outcomes") or []
        if not so:
            return ""
        EP = {"rfs": "RFS", "os": "OS", "landmark_rfs": "landmark RFS", "dfs": "DFS",
              "pfs": "PFS", "efs": "EFS", "ttr": "TTR", "other": "other"}
        rows = ""
        for s in so:
            med = "<b>not reached</b>" if s.get("not_reached") else (
                f"{s['median_value']:g} {s.get('time_unit','months')[:2]}"
                if s.get("median_value") is not None else DASH)
            comp = (f"HR {s['hazard_ratio']:g} ({esc(s.get('hr_ci',''))}); P={s.get('p_value')}"
                    if s.get("hazard_ratio") is not None else DASH)
            rows += (f"<tr><td class='mono strong'>{EP.get(s['endpoint'], s['endpoint'])}</td>"
                     f"<td>{esc(s.get('arm_label') or DASH)}</td>"
                     f"<td class='num'>{s.get('n_patients', DASH)}</td>"
                     f"<td class='num'>{med}</td><td class='mono dim'>{comp}</td>"
                     f"<td class='dim'>{esc(s.get('stratifier') or DASH)}</td></tr>")
        return ("<section><h2>Survival / time-to-event outcomes</h2>"
                "<table><thead><tr><th>endpoint</th><th>arm</th><th class='num'>n</th>"
                "<th class='num'>median</th><th>comparison</th><th>stratifier</th></tr></thead>"
                f"<tbody>{rows}</tbody></table></section>")

    def sec_preclinical():
        rows = ""
        for p in ex["patients"]:
            for pe in (p.get("preclinical_efficacy") or []):
                rows += (f"<tr><td class='pt'>{esc(p['paper_local_id'])}</td>"
                         f"<td class='mono'>{esc(pe['readout'])}</td>"
                         f"<td class='mono strong'>{esc(pe['result'])}</td>"
                         f"<td>{esc(pe.get('combination', DASH))}"
                         f"{(' / '+esc(pe['combination_detail'])) if pe.get('combination_detail') else ''}</td>"
                         f"<td class='dim'>{esc(pe.get('setting', DASH))}</td>"
                         f"<td class='mono dim'>{esc(pe.get('statistic') or DASH)}</td></tr>")
        if not rows:
            return ""
        return ("<section><h2>Preclinical antitumor efficacy</h2>"
                "<table><thead><tr><th>cohort</th><th>readout</th><th>result</th>"
                "<th>arm</th><th>setting</th><th>statistic</th></tr></thead>"
                f"<tbody>{rows}</tbody></table></section>")

    def _benefit_row(who, b):
        tp = esc(b.get("timepoint_label") or (b.get("timepoint_phase") or "").replace("_", " ") or DASH)
        ar = {True: "yes", False: "no"}.get(b.get("associated_with_response"), DASH)
        return (f"<tr><td class='pt'>{who}</td>"
                f"<td class='mono strong'>{esc(b['readout'].replace('_',' '))}</td>"
                f"<td class='mono'>{esc(b['direction'])}</td>"
                f"<td class='dim'>{tp}</td><td class='dim'>{ar}</td>"
                f"<td class='dim'>{esc(b.get('note') or DASH)}</td></tr>")

    def sec_benefit():
        rows = ""
        for p in ex["patients"]:
            for b in (p.get("clinical_benefit_signals") or []):
                rows += _benefit_row(esc(p["paper_local_id"]), b)
        for b in (ex.get("clinical_benefit_signals") or []):
            rows += _benefit_row("<i>cohort</i>", b)
        if not rows:
            return ""
        return ("<section><h2>Clinical benefit signals</h2>"
                "<table><thead><tr><th>cohort</th><th>readout</th><th>direction</th>"
                "<th>timepoint</th><th>assoc. response</th><th>note</th></tr></thead>"
                f"<tbody>{rows}</tbody></table></section>")

    def sec_mutations():
        mut = ex.get("neoantigen_mutations") or []
        if not mut:
            return ""
        rows = ""
        for m in mut:
            vaf = " → ".join(
                f"{(v.get('timepoint_label') or v.get('timepoint_phase') or '')}: {v['value']:g}"
                for v in (m.get("vaf") or []) if v.get("value") is not None)
            hla = ", ".join(m.get("hla_restrictions") or [])
            rows += (f"<tr><td class='pt'>{esc(m.get('patient_paper_id', DASH))}</td>"
                     f"<td class='mono strong'>{esc(m.get('gene_symbol') or DASH)}</td>"
                     f"<td class='mono dim'>{esc(m.get('genomic_change') or DASH)}</td>"
                     f"<td class='mono'>{esc(m.get('status') or DASH)}</td>"
                     f"<td class='dim'>{esc(m.get('clonality') or DASH)}</td>"
                     f"<td class='mono dim'>{esc(vaf or DASH)}</td>"
                     f"<td class='mono dim'>{esc(hla or DASH)}</td></tr>")
        return ("<section><h2>Neoantigen mutations</h2>"
                "<table><thead><tr><th>patient</th><th>gene</th><th>change</th>"
                "<th>status</th><th>clonality</th><th>VAF</th><th>HLA</th></tr></thead>"
                f"<tbody>{rows}</tbody></table></section>")

    def sec_funnel():
        cands = ex.get("candidates") or []
        reported = [
            ("predicted", ex.get("n_predicted_reported"), len(cands)),
            ("selected", ex.get("n_selected_reported"),
             sum(1 for c in cands if c.get("candidate_status") in ("selected", "administered"))),
            ("immunogenic (paper)", ex.get("n_immunogenic_reported"),
             sum(1 for e in ex["evidence"] if e.get("outcome") in ("immunogenic", "positive"))),
            ("tested negative (paper)", ex.get("n_tested_negative_reported"),
             sum(1 for e in ex["evidence"] if e.get("outcome") in ("not_immunogenic", "negative"))),
        ]
        has_rep = any(r[1] is not None for r in reported)
        statuses = collections.Counter(c.get("candidate_status") or "?" for c in cands)
        if not has_rep and not cands:
            return ""
        rows = ""
        for name, paper_n, rec_n in reported:
            if paper_n is None and rec_n == 0:
                continue
            pn = DASH if paper_n is None else paper_n
            rows += (f"<tr><td>{esc(name)}</td><td class='num'>{pn}</td>"
                     f"<td class='num'>{rec_n}</td></tr>")
        status_line = ""
        if statuses:
            bits = ", ".join(f"{k} {v}" for k, v in sorted(statuses.items()))
            status_line = f"<p class='dim'>Candidates in record: {esc(bits)} (status counts, not a dump).</p>"
        return ("<section><h2>Candidate funnel</h2>"
                "<table><thead><tr><th>count</th><th class='num'>paper stated</th>"
                "<th class='num'>in this record</th></tr></thead>"
                f"<tbody>{rows}</tbody></table>{status_line}</section>")

    def sec_cohort():
        R = responders(ex)
        rows = ""
        mixed = len({p.get("species", "human") for p in ex["patients"]}) > 1
        has_setting = any(p.get("trial_setting") for p in ex["patients"])
        for p in ex["patients"]:
            pid = p["paper_local_id"]
            r = pid in R
            sp = f"<td class='dim'>{esc(p.get('species','human'))}</td>" if mixed else ""
            st = (f"<td class='dim'>{esc((p.get('trial_setting') or DASH).replace('_',' '))}</td>"
                  if has_setting else "")
            badge = ("<span class='b resp'>responder</span>" if r
                     else "<span class='b non'>non-responder</span>")
            nr = " <span class='unv'>needs review</span>" if p.get("needs_review") else ""
            plat = (p.get("vaccine_platform") or "").replace("_", " ")
            rows += (f"<tr class='{'r' if r else ''}'><td class='pt'>{esc(pid)}{nr}</td>{sp}{st}"
                     f"<td>{badge}</td><td class='dim'>{esc(plat or DASH)}</td>"
                     f"<td class='num'>{p.get('n_peptides_synthesized', DASH)}</td>"
                     f"<td class='num'>{p.get('n_peptides_administered', DASH)}</td>"
                     f"<td class='num strong'>{p.get('n_peptides_immunogenic', DASH)}</td>"
                     f"<td class='num'>{len(p.get('hla_alleles') or [])}</td></tr>")
        sph = "<th>species</th>" if mixed else ""
        sth = "<th>setting</th>" if has_setting else ""
        table = ("<table><thead><tr><th>patient</th>" + sph + sth +
                 "<th>response</th><th>platform</th>"
                 "<th class='num'>synthesized</th><th class='num'>administered</th>"
                 "<th class='num'>immunogenic</th><th class='num'>HLA</th></tr></thead>"
                 f"<tbody>{rows}</tbody></table>")
        npat = len(ex["patients"])
        op = " open" if npat <= 12 else ""
        nr = len(R)
        return (f"<section><h2>Cohort</h2>"
                f"<details class='pcard'{op}><summary>{npat} patients "
                f"<span class='pn'>{nr} responder{'' if nr == 1 else 's'}</span></summary>"
                f"<div class='pbody'>{table}</div></details></section>")

    def sec_evidence():
        ev = ex.get("evidence") or []
        if not ev:
            return ""
        assays = collections.Counter(e.get("assay") or "?" for e in ev)
        outcomes = collections.Counter(e.get("outcome") or "?" for e in ev)
        mix = (", ".join(f"{k} {v}" for k, v in assays.most_common())
               + " · " + ", ".join(f"{k} {v}" for k, v in outcomes.most_common()))
        by_pt = collections.defaultdict(list)
        for e in ev:
            by_pt[e.get("patient_paper_id") or "?"].append(e)
        open_default = len(by_pt) <= 6
        op = " open" if open_default else ""
        blocks = ""
        for pid in sorted(by_pt, key=lambda x: (len(x), x)):
            rows = ""
            for e in by_pt[pid]:
                g, mut, ep1, ep2 = target_descr(e)
                ms = mag_str(e.get("magnitude"))
                mb = f"<span class='mg'>{esc(ms)}</span>" if ms else DASH
                oc = esc((e.get("outcome") or DASH).replace("_", " "))
                rows += (f"<tr><td class='mono'>{esc(g)}</td><td class='mono'>{esc(mut)}</td>"
                         f"<td class='mono hl'>{esc(ep1) or DASH}</td>"
                         f"<td class='mono'>{esc(ep2) or DASH}</td>"
                         f"<td class='dim'>{esc(e.get('assay') or DASH)}</td>"
                         f"<td class='mono'>{oc}</td><td class='mono'>{mb}</td></tr>")
            npos = sum(1 for e in by_pt[pid] if e.get("outcome") in ("immunogenic", "positive"))
            blocks += (f"<details class='pcard'{op}><summary>{esc(pid)} "
                       f"<span class='pn'>{len(by_pt[pid])} rows · {npos} immunogenic</span></summary>"
                       "<div class='pbody'><table class='na'><thead><tr><th>gene</th><th>mutation</th>"
                       "<th>MHC-I</th><th>MHC-II</th><th>assay</th><th>outcome</th>"
                       "<th>magnitude</th></tr></thead>"
                       f"<tbody>{rows}</tbody></table></div></details>")
        controls = ("<div class='rc'><span class='dim'>" + esc(mix) + "</span>"
                    "<span class='rcb'><button type='button' onclick=\"_pcards(true)\">Expand all</button>"
                    "<button type='button' onclick=\"_pcards(false)\">Collapse all</button></span></div>")
        return ("<section><h2>Evidence</h2>"
                "<p>All evidence rows in the record (immunogenic and named negatives). "
                "Magnitude is lossless: a number with unit, a table call 1/0, a grade, or raw.</p>"
                + controls + blocks + "</section>")

    def sec_tcr():
        methods = ex.get("tcr_seq_methods") or []
        deps = ex.get("data_depositions") or []
        flags = ex.get("neoantigen_tcr_flags") or []
        clones = ex.get("tcr_clonotypes") or []
        if not (methods or deps or flags or clones or ex.get("tcr_seq_status")):
            return ""
        parts = [f"<p class='dim'>tcr_seq_status={esc(ex.get('tcr_seq_status') or 'unset')}</p>"]
        if methods:
            rows = ""
            for m in methods:
                rows += (f"<tr><td class='mono'>{esc(m.get('method_local_id'))}</td>"
                         f"<td>{esc(m.get('modality'))}</td><td>{esc(m.get('platform'))}</td>"
                         f"<td>{esc(m.get('clonotype_definition'))}</td>"
                         f"<td class='dim'>{esc(', '.join(m.get('patient_scope') or []) or DASH)}</td></tr>")
            parts.append("<h3>Methods</h3><table><thead><tr><th>id</th><th>modality</th>"
                         "<th>platform</th><th>clonotype definition</th><th>patients</th>"
                         "</tr></thead>" f"<tbody>{rows}</tbody></table>")
        if deps:
            rows = ""
            for d in deps:
                rows += (f"<tr><td>{esc(d.get('repository'))}</td>"
                         f"<td class='mono'>{esc(d.get('accession'))}</td>"
                         f"<td>{esc(d.get('data_type'))}</td><td>{esc(d.get('access'))}</td></tr>")
            parts.append("<h3>Depositions</h3><table><thead><tr><th>repository</th>"
                         "<th>accession</th><th>data type</th><th>access</th></tr></thead>"
                         f"<tbody>{rows}</tbody></table>")
        if flags:
            rows = ""
            for f in flags:
                rows += (f"<tr><td class='pt'>{esc(f.get('patient_paper_id') or DASH)}</td>"
                         f"<td>{esc(f.get('target_kind'))}</td>"
                         f"<td class='mono'>{esc(tcr_target(f))}</td>"
                         f"<td class='mono'>{esc(f.get('tcr_identified'))}</td></tr>")
            parts.append("<h3>Per-target flags</h3><table><thead><tr><th>patient</th>"
                         "<th>target kind</th><th>target</th><th>TCR identified</th>"
                         "</tr></thead>" f"<tbody>{rows}</tbody></table>")
        if clones:
            shown = clones[:TCR_CLONE_CAP]
            more = len(clones) - len(shown)
            rows = ""
            for c in shown:
                nr = " <span class='unv'>needs review</span>" if c.get("needs_review") else ""
                rows += (f"<tr><td class='mono'>{esc(c.get('paper_local_id'))}{nr}</td>"
                         f"<td class='pt'>{esc(c.get('patient_paper_id') or DASH)}</td>"
                         f"<td class='mono'>{esc(chain_cdr3(c))}</td>"
                         f"<td class='mono'>{esc(tcr_target(c))}</td>"
                         f"<td class='dim'>{esc(c.get('specificity_evidence') or DASH)}</td></tr>")
            extra = (f"<p class='dim'>{more} more clonotypes in the JSON (not listed).</p>"
                     if more else "")
            op = " open" if len(shown) <= 12 else ""
            parts.append(
                f"<h3>Clonotypes</h3><details class='pcard'{op}>"
                f"<summary>{len(clones)} clonotypes "
                f"<span class='pn'>showing {len(shown)}</span></summary>"
                "<div class='pbody'><table class='na'><thead><tr><th>id</th><th>patient</th>"
                "<th>CDR3</th><th>target</th><th>specificity</th></tr></thead>"
                f"<tbody>{rows}</tbody></table>{extra}</div></details>")
        return "<section><h2>TCR</h2>" + "".join(parts) + "</section>"

    def sec_curator_notes():
        notes = ex.get("curator_notes") or []
        if not notes:
            return ""
        KLAB = {"challenge": "CHALLENGE", "decision": "DECISION",
                "caveat": "CAVEAT", "highlight": "HIGHLIGHT"}
        cards = ""
        for n in notes:
            unv = "<span class='unv'>unverified</span>" if n.get("needs_review") else ""
            ref_html = ""
            if n.get("refs"):
                chips = "".join(f"<span class='ref'>{esc(r)}</span>" for r in n["refs"])
                ref_html = f"<div class='refs'>{chips}</div>"
            cards += (f"<div class='note {esc(n['kind'])}'>"
                      f"<div class='ntag'>{KLAB.get(n['kind'], n['kind'].upper())}{unv}</div>"
                      f"<div class='nc'><p>{esc(n['text'])}</p>{ref_html}</div></div>")
        return ("<section><h2>Curator notes</h2>"
                "<p class='dim'>Rendered verbatim from the record (not auto-generated).</p>"
                + cards + "</section>")

    class2_note = (
        "<p class='dim' style='margin:12px 2px 0'>Epitope counts are predicted minimal "
        "epitopes (all class-I here). Class-II / CD4 responses are recorded as evidence "
        "on the long immunizing peptides, not as minimal epitopes.</p>"
        if (n2 == 0 and n_imp_ev) else ""
    )
    meta = " · ".join(filter(None, [
        f"PMID {ex['pmid']}" if ex.get("pmid") else "",
        f"DOI {ex['doi']}" if ex.get("doi") else "",
        f"NCT {ex['nct_id']}" if ex.get("nct_id") else "",
    ]))
    glance = [
        ("patients", len(ex["patients"])),
        ("peptides", len(ex["immunizing_peptides"])),
        ("epitopes", f"{n1}+{n2}" if (n1 or n2) else 0),
        ("evidence", len(ex["evidence"])),
        ("TCR clonotypes", n_tcr),
        ("needs review", n_nr),
    ]
    gcards = "".join(
        f"<div class='g'><div class='n'>{v}</div><div class='l'>{k}</div></div>"
        for k, v in glance)

    STYLE = """
:root{--paper:#FBF7F0;--ink:#26221E;--ink2:#5A534A;--rule:#E4DCCD;--card:#fff;--teal:#1F6F6B;--teal-w:#E2F0EE;--blue:#2C5A86;--blue-w:#E5EEF6;--green:#2F6B4F;--kras:#FBEFD6}
*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font-family:"Newsreader",Georgia,serif;font-size:18px;line-height:1.6}
.wrap{max-width:980px;margin:0 auto;padding:54px 30px 90px}
h1,h2,h3,h4{font-family:"Fraunces",Georgia,serif;font-weight:600;line-height:1.12}
code,.mono{font-family:"IBM Plex Mono",monospace}
.masthead{border-bottom:3px solid var(--ink);padding-bottom:22px;margin-bottom:28px}
.kicker{font-family:"IBM Plex Mono",monospace;font-size:12px;letter-spacing:.16em;text-transform:uppercase;color:var(--teal);font-weight:600}
h1{font-size:40px;font-weight:900;margin:.28em 0 .2em;letter-spacing:-.01em}.sub{font-size:20px;color:var(--ink2);font-style:italic;max-width:62ch}
.meta{font-family:"IBM Plex Mono",monospace;font-size:12px;color:var(--ink2);margin-top:16px}
.sv{display:inline-block;background:var(--green);color:#fff;border-radius:5px;padding:2px 9px;font-weight:600}
section{margin:42px 0}h2{font-size:13px;font-family:"IBM Plex Mono",monospace;letter-spacing:.16em;text-transform:uppercase;color:var(--ink2);border-bottom:1px solid var(--rule);padding-bottom:9px;margin-bottom:18px;font-weight:600}
h3{font-family:"IBM Plex Mono",monospace;font-size:12px;letter-spacing:.08em;text-transform:uppercase;color:var(--ink2);margin:18px 0 8px;font-weight:600}
p{margin:.5em 0 1em}
.glance{display:grid;grid-template-columns:repeat(auto-fit,minmax(130px,1fr));gap:14px;margin:22px 0}
.g{background:var(--card);border:1px solid var(--rule);border-radius:11px;padding:15px}.g .n{font-family:"Fraunces";font-weight:900;font-size:30px;color:var(--teal)}.g .l{font-size:13px;color:var(--ink2);margin-top:6px}
table{width:100%;border-collapse:collapse;margin:8px 0;font-size:15px}
th{font-family:"IBM Plex Mono",monospace;font-size:11px;letter-spacing:.04em;text-transform:uppercase;color:var(--ink2);text-align:left;border-bottom:1.5px solid var(--ink);padding:7px 9px;font-weight:600}
td{padding:6px 9px;border-bottom:1px solid var(--rule);vertical-align:top}
td.pt{font-family:"Fraunces";font-weight:600}td.num{text-align:right;font-family:"IBM Plex Mono",monospace}td.strong{font-weight:600;color:var(--teal)}tr.r{background:#FCFAF5}
.b{font-family:"IBM Plex Mono",monospace;font-size:11px;padding:2px 8px;border-radius:20px;font-weight:600;display:inline-block;margin:1px 3px 1px 0}.b.resp{background:var(--teal);color:#fff}.b.non{background:#EDE6D8;color:var(--ink2)}.b.soft{background:var(--teal-w);color:var(--teal)}
.rc{display:flex;justify-content:space-between;align-items:center;margin:6px 0 10px;gap:12px;flex-wrap:wrap}
.rcb button{font-family:"IBM Plex Mono",monospace;font-size:11px;color:var(--teal);background:var(--teal-w);border:1px solid var(--rule);border-radius:6px;padding:4px 10px;margin-left:7px;cursor:pointer}
.rcb button:hover{background:#D2E7E4}
details.pcard{background:var(--card);border:1px solid var(--rule);border-radius:11px;margin:10px 0;overflow:hidden}
details.pcard>summary{list-style:none;cursor:pointer;display:flex;justify-content:space-between;align-items:center;padding:13px 17px;font-family:"Fraunces",Georgia,serif;font-weight:600;font-size:18px}
details.pcard>summary::-webkit-details-marker{display:none}
details.pcard>summary::after{content:"\\25B8";color:var(--ink2);font-size:13px;margin-left:12px}
details.pcard[open]>summary::after{content:"\\25BE"}
details.pcard[open]>summary{border-bottom:1px solid var(--rule)}
details.pcard>summary:hover{background:#FCFAF5}
.pcard .pn{font-family:"IBM Plex Mono",monospace;font-size:12px;color:var(--teal);font-weight:600}
.pcard .pbody{padding:4px 17px 14px}
table.na td.mono{font-size:13px}.hl{background:#E7F1EE;border-radius:3px}.dim{color:var(--ink2);font-size:12.5px}
.mg{background:var(--teal-w);color:var(--teal);padding:1px 7px;border-radius:5px;font-size:12px;font-family:"IBM Plex Mono",monospace}
.note{display:flex;gap:13px;background:var(--card);border:1px solid var(--rule);border-radius:11px;padding:13px 16px;margin:10px 0}
.note .ntag{flex:0 0 100px;font-family:"IBM Plex Mono",monospace;font-size:10px;font-weight:600;letter-spacing:.04em;padding-top:2px;display:flex;flex-direction:column;gap:4px}
.note.highlight{border-left:4px solid var(--green)}.note.highlight .ntag{color:var(--green)}
.note.decision{border-left:4px solid var(--blue)}.note.decision .ntag{color:var(--blue)}
.note.caveat,.note.challenge{border-left:4px solid #B5701A}.note.caveat .ntag,.note.challenge .ntag{color:#B5701A}
.note .nc p{margin:0 0 6px;font-size:15.5px;color:#403a32;line-height:1.55}
.unv{font-family:"IBM Plex Mono",monospace;font-size:9px;background:#F6E0D8;color:#9A3B1A;padding:1px 5px;border-radius:4px;font-weight:600}
.refs{display:flex;gap:5px;flex-wrap:wrap}.ref{font-family:"IBM Plex Mono",monospace;font-size:10px;background:var(--blue-w);color:var(--blue);padding:1px 6px;border-radius:4px}
footer{margin-top:50px;border-top:3px solid var(--ink);padding-top:16px;font-family:"IBM Plex Mono",monospace;font-size:12px;color:var(--ink2);line-height:1.7}footer b{color:var(--ink)}
@media(max-width:720px){h1{font-size:30px}body{font-size:16.5px}}
"""
    HEAD = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>antVacDB · {esc(ex.get("pmid",""))}</title>'
            '<link rel="preconnect" href="https://fonts.googleapis.com">'
            '<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,400;9..144,600;9..144,900&family=Newsreader:ital,opsz@0,6..72;1,6..72&family=IBM+Plex+Mono:wght@400;500;600&display=swap" rel="stylesheet">'
            f"<style>{STYLE}</style></head><body><div class='wrap'>")
    MAST = ("<header class='masthead'><div class='kicker'>antVacDB · extraction review</div>"
            f"<h1>{esc(ex.get('title',''))}</h1>"
            f"<div class='sub'>{esc(ex.get('journal',''))} {ex.get('year','')} — {esc(ex.get('indication_summary',''))}</div>"
            f"<div class='meta'>{meta} · schema <span class='sv'>v{sv}</span></div></header>")
    GLANCE = f"<section><h2>At a glance</h2><div class='glance'>{gcards}</div>{class2_note}</section>"
    ov_txt = (", HARD overrides: " + ", ".join(hard_ov)) if hard_ov else ""
    FOOT = (f"<footer><b>Extraction.</b> {len(ex['patients'])} patients · "
            f"{len(ex['immunizing_peptides'])} immunizing peptides · {n1} MHC-I + {n2} MHC-II epitopes · "
            f"{len(ex['pools'])} pools · {len(ex['evidence'])} evidence rows ({nmag} with magnitude) · "
            f"{n_tcr} TCR clonotypes · {n_nr} needs_review{ov_txt}.<br>"
            f"Generated from the validated record by make_report.py; schema "
            f"<b>v{sv}</b>. Pure function of the JSON — no source re-reading. "
            f"Optional curator view; not a load or gold sign-off.</footer>")
    SCRIPT = ("<script>function _pcards(o){document.querySelectorAll('details.pcard')"
              ".forEach(function(d){d.open=o});}</script>")
    return (HEAD + MAST + GLANCE + sec_health() + sec_safety() + sec_survival()
            + sec_preclinical() + sec_benefit() + sec_mutations() + sec_funnel()
            + sec_cohort() + sec_evidence() + sec_tcr() + sec_curator_notes()
            + FOOT + SCRIPT + "</div></body></html>")


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        print("usage: vaxtract-report EXTRACTED.json [OUT.html]", file=sys.stderr)
        return 2
    src = argv[0]
    out = argv[1] if len(argv) > 1 else src.replace(".json", "_review.html")
    try:
        ex, sv = load_record(src)
    except Exception as e:
        print(f"REFUSING to render: record does not validate ({e})")
        return 1
    body = build_html(ex, sv)
    open(out, "w").write(body)
    n_tcr = len(ex.get("tcr_clonotypes") or [])
    print(f"wrote {out}  ({len(body)} bytes)  schema v{sv}  | "
          f"peptides={len(ex['immunizing_peptides'])} "
          f"evidence={len(ex['evidence'])} tcr={n_tcr}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
