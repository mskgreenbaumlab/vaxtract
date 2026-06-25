"""schema_digest.py — a compact, SDK-free field reference generated from the live schema.

Embedded in the agent's system prompt INSTEAD of the 70 KB schema.py source: it lists
every entity model's fields, types (Literals inline), required/optional, and the easy
constraints — but no docstrings or validator code. The add_entities/finalize tools
enforce the full contract, so the prompt only needs the shape. Generated at runtime so
it can never drift from schema.py.
"""
from __future__ import annotations

import inspect
import types as _types
import typing

# Document order: container, the entities it holds, then the value objects.
ENTITY_MODELS = [
    "ExtractedPaper", "ExtractedPatient", "ImmunizingPeptide", "MinimalEpitope",
    "ExtractedEvidence", "ExtractedPeptidePool", "SurvivalOutcome", "PreclinicalEfficacy",
    "ConcomitantImmunosuppression", "CuratorNote", "Measurement", "ResponseMagnitude", "Provenance",
    "NeoantigenCandidate", "PrioritizationScore", "VaccineDelivery", "CohortLatency",
]


def _annotated_metadata(ann):
    """Collect constraint metadata objects (e.g. StringConstraints) carried on the
    annotation itself, walking through Optional/list/Annotated wrappers. Needed because
    pydantic surfaces Annotated[str, StringConstraints(...)] constraints on the
    annotation rather than on FieldInfo.metadata."""
    found = []
    meta = getattr(ann, "__metadata__", None)
    if meta:
        found.extend(meta)
    for a in typing.get_args(ann):
        found.extend(_annotated_metadata(a))
    return found


def _type_str(ann) -> str:
    origin = typing.get_origin(ann)
    args = typing.get_args(ann)
    if ann is type(None):
        return "null"
    if origin is typing.Annotated or getattr(ann, "__metadata__", None) is not None:
        return _type_str(args[0])
    if origin is typing.Literal:
        return "one of[" + "|".join(str(a) for a in args) + "]"
    if origin in (typing.Union, getattr(_types, "UnionType", None)):
        non_none = [a for a in args if a is not type(None)]
        s = " or ".join(_type_str(a) for a in non_none)
        return s + ("?" if type(None) in args else "")
    if origin in (list, typing.List):
        return "list[" + _type_str(args[0]) + "]"
    if origin in (dict, typing.Dict):
        return "map"
    if inspect.isclass(ann):
        return ann.__name__
    return str(ann)


def _constraints(field_info) -> str:
    pairs = []
    seen = set()
    metas = list(field_info.metadata or []) + _annotated_metadata(field_info.annotation)
    for m in metas:
        for attr, label in (("pattern", "pattern="), ("ge", ">="), ("le", "<="),
                            ("gt", ">"), ("lt", "<"), ("max_length", "maxlen="),
                            ("min_length", "minlen=")):
            v = getattr(m, attr, None)
            if v is not None:
                tok = f"{label}{v}"
                if tok not in seen:
                    seen.add(tok)
                    pairs.append(tok)
    return (" {" + ", ".join(pairs) + "}") if pairs else ""


def build_schema_digest(schema_module, vocab_module=None) -> str:
    from pydantic import BaseModel

    out = [
        f"SCHEMA v{schema_module.SCHEMA_VERSION} — compact field reference.",
        "Notation: X? = optional; one of[a|b] = allowed values; list[X] = list of X.",
        "The add_entities/finalize tools enforce the FULL rules (regex, ranges, cross-entity).",
        "",
    ]
    for name in ENTITY_MODELS:
        model = getattr(schema_module, name, None)
        if not (inspect.isclass(model) and issubclass(model, BaseModel)):
            continue
        doc = inspect.getdoc(model) or ""
        purpose = doc.splitlines()[0] if doc else ""
        out.append(f"## {name} — {purpose}")
        for fname, fi in model.model_fields.items():
            req = "req" if fi.is_required() else "opt"
            out.append(f"  {fname}: {_type_str(fi.annotation)} [{req}]{_constraints(fi)}")
        out.append("")
    return "\n".join(out)
