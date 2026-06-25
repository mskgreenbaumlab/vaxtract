#!/usr/bin/env python3
"""antVacDB — human-readable report generator.

Renders the editorial HTML review from a VALIDATED ExtractedPaper JSON. The report
is a pure function of the JSON: it never re-reads the paper or supplementary files,
so the human view and the machine view can never disagree. Every section is
conditional on what the record actually contains, so one generator serves every
paper (human/mouse, peptide/RNA/DNA, with or without survival/efficacy/magnitude).

Usage:
    python3 make_report.py EXTRACTED.json [OUT.html]

Pipeline position:  extraction -> validated JSON -> (this) -> report.html
"""
from __future__ import annotations
import sys, json, html, re, os, collections

# Em dash, as a name so it never appears as a backslash escape *inside* an
# f-string expression (a SyntaxError on Python < 3.12).
DASH = "—"

# --- 0. load + VALIDATE (the report only ever runs on a record the schema accepts)
SRC = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else SRC.replace(".json", "_review.html")
HERE = os.path.dirname(os.path.abspath(SRC)) or "."
sys.path.insert(0, HERE); sys.path.insert(0, os.getcwd())
try:
    import schema
    from schema import ExtractedPaper
    ex = ExtractedPaper(**json.loads(open(SRC).read())).model_dump()
    SV = schema.SCHEMA_VERSION
except Exception as e:
    print(f"REFUSING to render: record does not validate ({e})"); sys.exit(1)

def esc(s): return html.escape(str(s)) if s is not None else ""

# --- indices -----------------------------------------------------------------
imp_by_id = {i["paper_local_id"]: i for i in ex["immunizing_peptides"]}
epi_by_id = {e["paper_local_id"]: e for e in ex["epitopes"]}
pool_by_id = {p["paper_local_id"]: p for p in ex["pools"]}
# epitopes grouped by parent IMP + class
epi_by_parent = collections.defaultdict(dict)
for e in ex["epitopes"]:
    for pid in e.get("parent_peptide_ids", []):
        epi_by_parent[pid][e["mhc_class"]] = e

def responders():
    r = set()
    for p in ex["patients"]:
        if (p.get("n_peptides_immunogenic") or 0) > 0: r.add(p["paper_local_id"])
    for e in ex["evidence"]:
        if e.get("outcome") in ("immunogenic", "positive"): r.add(e["patient_paper_id"])
    return r

def mag_str(m):
    if not m: return ""
    if m.get("value") is not None:
        u = {"sfc_per_1e6": "SFC/10\u2076", "percent_of_parent": "%", "stimulation_index": "SI"}.get(m["unit"], m["unit"])
        return f"{m['value']:g} {u}"
    if m.get("grade"): return m["grade"].replace("_", " ")
    mm = re.search(r"\[([\d,\s.]+)\]", m.get("raw", "") or "")
    return (mm.group(1).replace(" ", "") + " SFC/10\u2076") if mm else "(set)"

# --- section builders (each returns "" when its data is absent) ---------------
def sec_survival():
    so = ex.get("survival_outcomes") or []
    if not so: return ""
    EP = {"rfs": "RFS", "os": "OS", "landmark_rfs": "landmark RFS", "dfs": "DFS",
          "pfs": "PFS", "efs": "EFS", "ttr": "TTR", "other": "other"}
    rows = ""
    for s in so:
        med = "<b>not reached</b>" if s.get("not_reached") else (
              f"{s['median_value']:g} {s.get('time_unit','months')[:2]}" if s.get("median_value") is not None else "\u2014")
        comp = (f"HR {s['hazard_ratio']:g} ({esc(s.get('hr_ci',''))}); P={s.get('p_value')}"
                if s.get("hazard_ratio") is not None else "\u2014")
        rows += (f"<tr><td class='mono strong'>{EP.get(s['endpoint'], s['endpoint'])}</td>"
                 f"<td>{esc(s.get('arm_label') or DASH)}</td><td class='num'>{s.get('n_patients',DASH)}</td>"
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
                     f"<td class='mono'>{esc(pe['readout'])}</td><td class='mono strong'>{esc(pe['result'])}</td>"
                     f"<td>{esc(pe.get('combination',DASH))}{(' / '+esc(pe['combination_detail'])) if pe.get('combination_detail') else ''}</td>"
                     f"<td class='dim'>{esc(pe.get('setting',DASH))}</td>"
                     f"<td class='mono dim'>{esc(pe.get('statistic') or DASH)}</td></tr>")
    if not rows: return ""
    return ("<section><h2>Preclinical antitumor efficacy</h2>"
            "<table><thead><tr><th>cohort</th><th>readout</th><th>result</th><th>arm</th>"
            "<th>setting</th><th>statistic</th></tr></thead>"
            f"<tbody>{rows}</tbody></table></section>")

def sec_benefit():
    # v2.10 P22: response->benefit bridge signals (per-patient + paper-level cohort-aggregate).
    rows = ""
    for p in ex["patients"]:
        for b in (p.get("clinical_benefit_signals") or []):
            rows += _benefit_row(esc(p["paper_local_id"]), b)
    for b in (ex.get("clinical_benefit_signals") or []):
        rows += _benefit_row("<i>cohort</i>", b)
    if not rows: return ""
    return ("<section><h2>Clinical benefit signals (response→benefit bridge)</h2>"
            "<table><thead><tr><th>cohort</th><th>readout</th><th>direction</th>"
            "<th>timepoint</th><th>assoc. response</th><th>note</th></tr></thead>"
            f"<tbody>{rows}</tbody></table></section>")

def _benefit_row(who, b):
    tp = esc(b.get("timepoint_label") or (b.get("timepoint_phase") or "").replace("_", " ") or DASH)
    ar = {True: "yes", False: "no"}.get(b.get("associated_with_response"), DASH)
    return (f"<tr><td class='pt'>{who}</td>"
            f"<td class='mono strong'>{esc(b['readout'].replace('_',' '))}</td>"
            f"<td class='mono'>{esc(b['direction'])}</td>"
            f"<td class='dim'>{tp}</td><td class='dim'>{ar}</td>"
            f"<td class='dim'>{esc(b.get('note') or DASH)}</td></tr>")

def sec_mutations():
    # v2.11 P20: gene-level neoantigen mutations (clonality / VAF dynamics / HLA).
    mut = ex.get("neoantigen_mutations") or []
    if not mut: return ""
    rows = ""
    for m in mut:
        vaf = " → ".join(f"{(v.get('timepoint_label') or v.get('timepoint_phase') or '')}: "
                         f"{v['value']:g}" for v in (m.get("vaf") or []) if v.get("value") is not None)
        hla = ", ".join(m.get("hla_restrictions") or [])
        rows += (f"<tr><td class='pt'>{esc(m.get('patient_paper_id',DASH))}</td>"
                 f"<td class='mono strong'>{esc(m.get('gene_symbol') or DASH)}</td>"
                 f"<td class='mono dim'>{esc(m.get('genomic_change') or DASH)}</td>"
                 f"<td class='mono'>{esc(m.get('status') or DASH)}</td>"
                 f"<td class='dim'>{esc(m.get('clonality') or DASH)}</td>"
                 f"<td class='mono dim'>{esc(vaf or DASH)}</td>"
                 f"<td class='mono dim'>{esc(hla or DASH)}</td></tr>")
    return ("<section><h2>Neoantigen mutations (clonality / antigen dynamics)</h2>"
            "<table><thead><tr><th>patient</th><th>gene</th><th>change</th><th>status</th>"
            "<th>clonality</th><th>VAF</th><th>HLA</th></tr></thead>"
            f"<tbody>{rows}</tbody></table></section>")

def sec_cohort():
    R = responders(); rows = ""; mixed = len({p.get("species","human") for p in ex["patients"]}) > 1
    has_setting = any(p.get("trial_setting") for p in ex["patients"])
    for p in ex["patients"]:
        pid = p["paper_local_id"]; r = pid in R
        sp = f"<td class='dim'>{esc(p.get('species','human'))}</td>" if mixed else ""
        st = f"<td class='dim'>{esc((p.get('trial_setting') or DASH).replace('_',' '))}</td>" if has_setting else ""
        badge = (f"<span class='b resp'>responder</span>" if r else "<span class='b non'>non-responder</span>")
        rows += (f"<tr class='{'r' if r else ''}'><td class='pt'>{esc(pid)}</td>{sp}{st}<td>{badge}</td>"
                 f"<td class='num'>{p.get('n_peptides_synthesized',DASH)}</td>"
                 f"<td class='num strong'>{p.get('n_peptides_immunogenic',DASH)}</td>"
                 f"<td class='num'>{len(p.get('hla_alleles') or [])}</td></tr>")
    sph = "<th>species</th>" if mixed else ""
    sth = "<th>setting</th>" if has_setting else ""
    table = ("<table><thead><tr><th>patient</th>" + sph + sth +
             "<th>response</th><th class='num'>peptides</th><th class='num'>immunogenic</th>"
             "<th class='num'>HLA</th></tr></thead>" f"<tbody>{rows}</tbody></table>")
    # Collapsible like the per-patient cards: the cohort table is one row per patient,
    # so it collapses by default once the cohort is large and stays open when small.
    npat = len(ex["patients"]); op = " open" if npat <= 12 else ""
    nr = len(R)
    return (f"<section><h2>Cohort</h2>"
            f"<details class='pcard'{op}><summary>{npat} patients "
            f"<span class='pn'>{nr} responder{'' if nr == 1 else 's'}</span></summary>"
            f"<div class='pbody'>{table}</div></details></section>")

def target_descr(e):
    """Resolve an evidence row's target to (gene, mutation, mhc-I epi, mhc-II epi)."""
    if e["target_kind"] == "immunizing_peptide":
        i = imp_by_id.get(e["immunizing_peptide_paper_id"], {}); key = i.get("paper_local_id")
    elif e["target_kind"] == "epitope":
        i = epi_by_id.get(e["epitope_paper_id"], {}); key = (i.get("parent_peptide_ids") or [None])[0]
    else:
        p = pool_by_id.get(e["pool_paper_id"], {}); 
        return (f"pool ({len(p.get('member_peptide_ids',[]))} peptides)", "", "", "")
    e1 = epi_by_parent.get(key, {}).get("I", {}); e2 = epi_by_parent.get(key, {}).get("II", {})
    return (i.get("gene_symbol", "\u2014"), i.get("mutation", "\u2014"),
            e1.get("sequence", ""), e2.get("sequence", ""))

def sec_responders():
    R = responders()
    by_pt = collections.defaultdict(list)
    for e in ex["evidence"]:
        if e.get("outcome") in ("immunogenic", "positive"):
            by_pt[e["patient_paper_id"]].append(e)
    if not by_pt: return ""
    # Per-patient cards are COLLAPSIBLE (native <details>). They collapse by default
    # once the cohort is large so the section stays scannable as patient count grows;
    # for a small cohort they start open (no clicking needed). The expand/collapse-all
    # controls flip every card at once. open=open_default.
    open_default = len(by_pt) <= 6
    op = " open" if open_default else ""
    blocks = ""
    for pid in sorted(by_pt, key=lambda x: (len(x), x)):
        rows = ""
        for e in by_pt[pid]:
            g, mut, ep1, ep2 = target_descr(e); ms = mag_str(e.get("magnitude"))
            mb = f"<span class='mg'>{ms}</span>" if ms else ""
            kr = " class='kras'" if g == "KRAS" else ""
            rows += (f"<tr{kr}><td class='mono'>{esc(g)}</td><td class='mono'>{esc(mut)}</td>"
                     f"<td class='mono hl'>{esc(ep1) or DASH}</td><td class='mono'>{esc(ep2) or DASH}</td>"
                     f"<td class='mono'>{mb}</td></tr>")
        blocks += (f"<details class='pcard'{op}><summary>{esc(pid)} "
                   f"<span class='pn'>{len(by_pt[pid])} immunogenic</span></summary>"
                   "<div class='pbody'><table class='na'><thead><tr><th>gene</th><th>mutation</th>"
                   "<th>MHC-I epitope</th><th>MHC-II epitope</th><th>magnitude</th></tr></thead>"
                   f"<tbody>{rows}</tbody></table></div></details>")
    controls = ("<div class='rc'><span class='dim'>" + str(len(by_pt)) + " responder"
                + ("s" if len(by_pt) != 1 else "") + "</span>"
                "<span class='rcb'><button type='button' onclick=\"_pcards(true)\">Expand all</button>"
                "<button type='button' onclick=\"_pcards(false)\">Collapse all</button></span></div>")
    return ("<section><h2>Immunogenic responses, by patient</h2>"
            "<p>Resolved from the validated evidence rows. Magnitude shown where the record carries one "
            "(numeric, ordinal grade, or a lossless set). Click a patient to expand.</p>"
            + controls + blocks + "</section>")

def sec_curator_notes():
    """Render curator_notes VERBATIM. No model call, no generation — the text is
    printed exactly as stored; this section is the reason the field exists."""
    notes = ex.get("curator_notes") or []
    if not notes: return ""
    KLAB = {"challenge": "CHALLENGE", "decision": "DECISION", "caveat": "CAVEAT", "highlight": "HIGHLIGHT"}
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
    return ("<section><h2>Curator notes &amp; decisions</h2>"
            "<p class='dim'>Curatorial commentary, rendered verbatim from the record (not "
            "auto-generated). Items marked <span class='unv'>unverified</span> await human "
            "sign-off.</p>" + cards + "</section>")


# --- assemble ----------------------------------------------------------------
n1 = sum(1 for e in ex["epitopes"] if e["mhc_class"] == "I")
n2 = sum(1 for e in ex["epitopes"] if e["mhc_class"] == "II")
# When the minimal-epitope layer is all class-I but the paper reports CD4/class-II responses
# (recorded as evidence on the long immunizing peptides, per the schema's P4 design), say so —
# else "0 MHC-II epitopes" reads as "no class-II biology", which would be misleading (Keskin).
n_imp_evidence = sum(1 for e in ex["evidence"] if e.get("target_kind") == "immunizing_peptide")
class2_note = (
    "<p style='color:var(--ink2);font-size:14.5px;font-style:italic;margin:12px 2px 0'>"
    "Epitope counts are predicted <b>minimal</b> epitopes (all class-I here). Class-II / CD4 "
    "T-cell responses are recorded as evidence on the long immunizing peptides — see the "
    "evidence rows — not as minimal epitopes.</p>"
    if (n2 == 0 and n_imp_evidence) else ""
)
nmag = sum(1 for e in ex["evidence"] if e.get("magnitude"))
nresp = len(responders())
meta = " · ".join(filter(None, [
    f"PMID {ex['pmid']}" if ex.get("pmid") else "", f"DOI {ex['doi']}" if ex.get("doi") else "",
    f"NCT {ex['nct_id']}" if ex.get("nct_id") else ""]))
glance = [("patients", len(ex["patients"])), ("responders", nresp),
          ("epitopes", f"{n1}+{n2}" if (n1 or n2) else 0), ("magnitudes", nmag)]
gcards = "".join(f"<div class='g'><div class='n'>{v}</div><div class='l'>{k}</div></div>" for k, v in glance)

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
p{margin:.5em 0 1em}.lead{font-size:20px}
.glance{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:22px 0}
.g{background:var(--card);border:1px solid var(--rule);border-radius:11px;padding:15px}.g .n{font-family:"Fraunces";font-weight:900;font-size:30px;color:var(--teal)}.g .l{font-size:13px;color:var(--ink2);margin-top:6px}
table{width:100%;border-collapse:collapse;margin:8px 0;font-size:15px}
th{font-family:"IBM Plex Mono",monospace;font-size:11px;letter-spacing:.04em;text-transform:uppercase;color:var(--ink2);text-align:left;border-bottom:1.5px solid var(--ink);padding:7px 9px;font-weight:600}
td{padding:6px 9px;border-bottom:1px solid var(--rule);vertical-align:top}
td.pt{font-family:"Fraunces";font-weight:600}td.num{text-align:right;font-family:"IBM Plex Mono",monospace}td.strong{font-weight:600;color:var(--teal)}tr.r{background:#FCFAF5}
.b{font-family:"IBM Plex Mono",monospace;font-size:11px;padding:2px 8px;border-radius:20px;font-weight:600}.b.resp{background:var(--teal);color:#fff}.b.non{background:#EDE6D8;color:var(--ink2)}
.rc{display:flex;justify-content:space-between;align-items:center;margin:6px 0 10px}
.rcb button{font-family:"IBM Plex Mono",monospace;font-size:11px;color:var(--teal);background:var(--teal-w);border:1px solid var(--rule);border-radius:6px;padding:4px 10px;margin-left:7px;cursor:pointer}
.rcb button:hover{background:#D2E7E4}
details.pcard{background:var(--card);border:1px solid var(--rule);border-radius:11px;margin:10px 0;overflow:hidden}
details.pcard>summary{list-style:none;cursor:pointer;display:flex;justify-content:space-between;align-items:center;padding:13px 17px;font-family:"Fraunces",Georgia,serif;font-weight:600;font-size:18px}
details.pcard>summary::-webkit-details-marker{display:none}
details.pcard>summary::after{content:"\\25B8";color:var(--ink2);font-size:13px;margin-left:12px;transition:transform .15s}
details.pcard[open]>summary::after{content:"\\25BE"}
details.pcard[open]>summary{border-bottom:1px solid var(--rule)}
details.pcard>summary:hover{background:#FCFAF5}
.pcard .pn{font-family:"IBM Plex Mono",monospace;font-size:12px;color:var(--teal);font-weight:600}
.pcard .pbody{padding:4px 17px 14px}
table.na td.mono{font-size:13px}.hl{background:#E7F1EE;border-radius:3px}.dim{color:var(--ink2);font-size:12.5px}
.mg{background:var(--teal-w);color:var(--teal);padding:1px 7px;border-radius:5px;font-size:12px;font-family:"IBM Plex Mono",monospace}
tr.kras td{background:var(--kras)}
.note{display:flex;gap:13px;background:var(--card);border:1px solid var(--rule);border-radius:11px;padding:13px 16px;margin:10px 0}
.note .ntag{flex:0 0 100px;font-family:"IBM Plex Mono",monospace;font-size:10px;font-weight:600;letter-spacing:.04em;padding-top:2px;display:flex;flex-direction:column;gap:4px}
.note.highlight{border-left:4px solid var(--green)}.note.highlight .ntag{color:var(--green)}
.note.decision{border-left:4px solid var(--blue)}.note.decision .ntag{color:var(--blue)}
.note.caveat,.note.challenge{border-left:4px solid #B5701A}.note.caveat .ntag,.note.challenge .ntag{color:#B5701A}
.note .nc p{margin:0 0 6px;font-size:15.5px;color:#403a32;line-height:1.55}
.unv{font-family:"IBM Plex Mono",monospace;font-size:9px;background:#F6E0D8;color:#9A3B1A;padding:1px 5px;border-radius:4px;font-weight:600}
.refs{display:flex;gap:5px;flex-wrap:wrap}.ref{font-family:"IBM Plex Mono",monospace;font-size:10px;background:var(--blue-w);color:var(--blue);padding:1px 6px;border-radius:4px}
footer{margin-top:50px;border-top:3px solid var(--ink);padding-top:16px;font-family:"IBM Plex Mono",monospace;font-size:12px;color:var(--ink2);line-height:1.7}footer b{color:var(--ink)}
@media(max-width:720px){.glance{grid-template-columns:repeat(2,1fr)}h1{font-size:30px}body{font-size:16.5px}}
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
        f"<div class='meta'>{meta} · schema <span class='sv'>v{SV}</span></div></header>")

GLANCE = f"<section><h2>At a glance</h2><div class='glance'>{gcards}</div>{class2_note}</section>"

FOOT = (f"<footer><b>Extraction.</b> {len(ex['patients'])} patients · "
        f"{len(ex['immunizing_peptides'])} immunizing peptides · {n1} MHC-I + {n2} MHC-II epitopes · "
        f"{len(ex['pools'])} pools · {len(ex['evidence'])} evidence rows ({nmag} with magnitude) · "
        f"{len(ex.get('survival_outcomes') or [])} survival outcomes.<br>"
        f"Generated from the validated record by make_report.py; conforms to antVacDB schema "
        f"<b>v{SV}</b>. The report is a pure function of the JSON — no source re-reading.</footer>")

# tiny, dependency-free toggle for the collapsible per-patient cards (expand/collapse all)
SCRIPT = ("<script>function _pcards(o){document.querySelectorAll('details.pcard')"
          ".forEach(function(d){d.open=o});}</script>")

body = (HEAD + MAST + GLANCE + sec_survival() + sec_preclinical() + sec_benefit() + sec_mutations()
        + sec_cohort() + sec_responders() + sec_curator_notes() + FOOT + SCRIPT + "</div></body></html>")
open(OUT, "w").write(body)
print(f"wrote {OUT}  ({len(body)} bytes)  schema v{SV}  | sections: "
      f"survival={'Y' if (ex.get('survival_outcomes')) else '-'} "
      f"preclinical={'Y' if any(p.get('preclinical_efficacy') for p in ex['patients']) else '-'} "
      f"magnitudes={nmag}")
