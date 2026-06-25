#!/usr/bin/env python3
"""schema_overview.py — render a human-readable HTML overview of the extraction schema.

For collaborators (not machines): introspects the packet's pydantic models
(cancervac_packet/schema.py) + controlled vocabularies (vocab.py) and emits a
single self-contained HTML page describing WHAT gets extracted from each paper —
the entities, their fields/types, and the allowed category values.

Usage:
    python tools/schema_overview.py [OUT.html]      # default: schema_overview.html
"""
from __future__ import annotations

import html
import inspect
import pathlib
import sys
import typing

PKG = pathlib.Path(__file__).resolve().parents[1] / "cancervac_packet"
sys.path.insert(0, str(PKG))

import schema  # noqa: E402
import vocab  # noqa: E402
from pydantic import BaseModel  # noqa: E402

OUT = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(__file__).resolve().parents[1] / "schema_overview.html"

# Models in document order: top-level container first, then the entities it holds,
# then the value objects. Skip the private base classes (_Frozen/_Extracted/_PeptideCore).
ORDER = [
    "ExtractedPaper", "ExtractedPatient", "ImmunizingPeptide", "MinimalEpitope",
    "ExtractedEvidence", "ExtractedPeptidePool", "SurvivalOutcome",
    "PreclinicalEfficacy", "ConcomitantImmunosuppression", "CuratorNote",
    "Measurement", "ResponseMagnitude", "Provenance",
]
MODELS = {n: getattr(schema, n) for n in ORDER if hasattr(schema, n)}
MODEL_NAMES = set(MODELS)


def esc(x) -> str:
    return html.escape(str(x))


def anchor(name: str) -> str:
    return "m-" + name.lower()


def _constraint_hint(metadata) -> str:
    """Pull the human-useful bits (pattern, length bounds) out of pydantic
    StringConstraints/Field metadata; drop the StringConstraints(...=None) noise."""
    pat = mn = mx = None
    for m in metadata:
        pat = getattr(m, "pattern", None) or pat
        if getattr(m, "min_length", None) is not None:
            mn = m.min_length
        if getattr(m, "max_length", None) is not None:
            mx = m.max_length
    bits = []
    if pat:
        bits.append(f"format <code>{esc(pat)}</code>")
    if mn is not None and mx is not None:
        bits.append(f"len {mn}–{mx}")
    elif mx is not None:
        bits.append(f"≤{mx} chars")
    elif mn is not None:
        bits.append(f"≥{mn} chars")
    return f" <span class='fmt'>{' · '.join(bits)}</span>" if bits else ""


def type_str(ann) -> str:
    """Render a type annotation as readable HTML (Literals inline, models linked)."""
    if ann is type(None):
        return "null"
    # Annotated[str, StringConstraints(...)] -> render the base type + a small hint,
    # not the full StringConstraints(...=None) repr.
    if hasattr(ann, "__metadata__"):
        return type_str(ann.__origin__) + _constraint_hint(ann.__metadata__)

    origin = typing.get_origin(ann)
    args = typing.get_args(ann)
    if origin is typing.Literal:
        vals = " · ".join(f"<code>{esc(a)}</code>" for a in args)
        return f"<span class='lit'>one of: {vals}</span>"
    if origin in (typing.Union, getattr(__import__("types"), "UnionType", None)):
        non_none = [a for a in args if a is not type(None)]
        rendered = " or ".join(type_str(a) for a in non_none)
        return rendered + (" <span class='opt'>(optional)</span>" if type(None) in args else "")
    if origin in (list, typing.List):
        return "list of " + type_str(args[0])
    if origin in (dict, typing.Dict):
        return "map"
    if inspect.isclass(ann) and issubclass(ann, BaseModel) and ann.__name__ in MODEL_NAMES:
        n = ann.__name__
        return f"<a href='#{anchor(n)}'>{esc(n)}</a>"
    if inspect.isclass(ann):
        return f"<span class='ty'>{esc(ann.__name__)}</span>"
    return f"<span class='ty'>{esc(ann)}</span>"


def field_rows(model) -> str:
    rows = []
    for fname, fi in model.model_fields.items():
        req = "<span class='req'>required</span>" if fi.is_required() else "<span class='optn'>optional</span>"
        default = ""
        if not fi.is_required() and fi.default is not None and repr(fi.default) != "PydanticUndefined":
            default = f"<span class='def'>default {esc(fi.default)}</span>"
        desc = esc(fi.description) if fi.description else ""
        rows.append(
            f"<tr><td class='fn'><code>{esc(fname)}</code></td>"
            f"<td class='ft'>{type_str(fi.annotation)}</td>"
            f"<td class='fr'>{req} {default}</td>"
            f"<td class='fd'>{desc}</td></tr>"
        )
    return "\n".join(rows)


def model_card(name: str, model) -> str:
    doc = inspect.getdoc(model) or ""
    return f"""
    <section class="card" id="{anchor(name)}">
      <h3>{esc(name)} <span class="nf">{len(model.model_fields)} fields</span></h3>
      <p class="doc">{esc(doc)}</p>
      <table class="fields">
        <thead><tr><th>field</th><th>type</th><th>presence</th><th>notes</th></tr></thead>
        <tbody>{field_rows(model)}</tbody>
      </table>
    </section>"""


def vocab_section() -> str:
    tuples = [(k, v) for k, v in sorted(vars(vocab).items())
              if isinstance(v, tuple) and k.isupper()]
    blocks = []
    for name, vals in tuples:
        chips = " ".join(f"<code>{esc(x)}</code>" for x in vals)
        blocks.append(f"<div class='vocab'><div class='vn'>{esc(name)}</div><div class='vc'>{chips}</div></div>")
    return "\n".join(blocks)


def toc() -> str:
    links = " · ".join(f"<a href='#{anchor(n)}'>{esc(n)}</a>" for n in MODELS)
    return f"<p class='toc'>{links} · <a href='#vocab'>Controlled vocabularies</a></p>"


HTML = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>antVacDB extraction schema (v{schema.SCHEMA_VERSION})</title>
<style>
  :root {{ --ink:#1a1a1a; --dim:#666; --line:#e2e2e2; --accent:#0b5; --bg:#fafafa; }}
  * {{ box-sizing:border-box; }}
  body {{ font-family:Arial,Helvetica,sans-serif; color:var(--ink); background:var(--bg);
         margin:0; line-height:1.5; }}
  .wrap {{ max-width:1040px; margin:0 auto; padding:32px 24px 80px; }}
  h1 {{ font-size:30px; margin:0 0 4px; }}
  .sub {{ color:var(--dim); margin:0 0 18px; }}
  .intro {{ background:#fff; border:1px solid var(--line); border-radius:10px; padding:16px 20px; margin:18px 0; }}
  .intro h2 {{ font-size:16px; margin:0 0 8px; }}
  .intro ul {{ margin:6px 0 0; padding-left:20px; }}
  .toc {{ font-size:13px; color:var(--dim); line-height:2; }}
  .toc a {{ color:#06c; text-decoration:none; }}
  h2.sec {{ font-size:20px; border-bottom:2px solid var(--ink); padding-bottom:4px; margin:34px 0 10px; }}
  .card {{ background:#fff; border:1px solid var(--line); border-radius:10px; padding:16px 20px; margin:14px 0; }}
  .card h3 {{ font-size:18px; margin:0 0 6px; }}
  .nf {{ font-size:12px; font-weight:normal; color:var(--dim); }}
  .doc {{ color:#333; margin:0 0 12px; white-space:pre-wrap; }}
  table.fields {{ border-collapse:collapse; width:100%; font-size:13px; }}
  table.fields th {{ text-align:left; color:var(--dim); font-weight:600; border-bottom:1px solid var(--line);
                     padding:4px 8px; }}
  table.fields td {{ border-bottom:1px solid #f0f0f0; padding:5px 8px; vertical-align:top; }}
  td.fn code {{ font-weight:700; }}
  code {{ font-family:"SF Mono",Menlo,Consolas,monospace; background:#f3f3f3; padding:1px 5px;
          border-radius:4px; font-size:12px; }}
  a {{ color:#06c; }}
  .lit code {{ background:#eef7f0; }}
  .opt, .optn {{ color:var(--dim); font-style:italic; }}
  .req {{ color:var(--accent); font-weight:600; }}
  .def {{ color:var(--dim); }}
  .ty {{ color:#555; }}
  .fmt {{ color:var(--dim); font-size:11px; }}
  .fmt code {{ background:#f7f7f7; font-size:10px; word-break:break-all; }}
  .vocab {{ display:flex; gap:14px; padding:8px 0; border-bottom:1px solid #f0f0f0; align-items:baseline; }}
  .vn {{ flex:0 0 230px; font-weight:700; font-size:13px; }}
  .vc {{ flex:1; }} .vc code {{ margin:2px; display:inline-block; }}
</style></head><body><div class="wrap">
  <h1>antVacDB — what we extract from each paper</h1>
  <p class="sub">Pydantic schema <b>v{schema.SCHEMA_VERSION}</b> · {len(MODELS)} entity types · {len([k for k,v in vars(vocab).items() if isinstance(v,tuple) and k.isupper()])} controlled vocabularies. One paper → one <code>ExtractedPaper</code> record.</p>

  <div class="intro">
    <h2>How to read this</h2>
    <ul>
      <li>Every paper becomes one <b>ExtractedPaper</b>, which holds lists of the entities below
          (patients, peptides, epitopes, evidence rows, survival outcomes…).</li>
      <li><b>Two peptide levels:</b> an <b>ImmunizingPeptide</b> is the long peptide actually given;
          a <b>MinimalEpitope</b> is the short MHC-binding stretch within it (one long peptide can yield several).</li>
      <li><b>Nothing is invented.</b> A value not stated in the paper is left empty. Numbers read off a
          figure are kept lossless (the raw reading is preserved, flagged low-confidence, never a clean number).</li>
      <li><b>Controlled vocabularies</b> (bottom) are the only allowed values for the categorical fields —
          they keep extractions comparable across papers.</li>
    </ul>
    {toc()}
  </div>

  <h2 class="sec">Entities</h2>
  {''.join(model_card(n, m) for n, m in MODELS.items())}

  <h2 class="sec" id="vocab">Controlled vocabularies</h2>
  <p class="sub">The allowed values for each categorical axis (the extractor cannot use anything else).</p>
  {vocab_section()}
</div></body></html>"""

OUT.write_text(HTML)
print(f"wrote {OUT}  ({len(HTML)} bytes)  schema v{schema.SCHEMA_VERSION}  · {len(MODELS)} models")
