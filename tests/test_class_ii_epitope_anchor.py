"""v2.12 P23: class-II MinimalEpitope minting gated on a quoted restriction anchor; mouse class-II regex."""
import re
import pytest
from cancervac_packet.schema import MinimalEpitope, _MHC_PATTERN


def _epi(**kw):
    base = dict(quoted_text="x", section_ref="Table S5", paper_local_id="E1",
                sequence="SIINFEKLAA", is_neoantigen=True, mhc_class="II")
    base.update(kw)
    return MinimalEpitope(**base)


def test_class_ii_named_allele_ok():
    e = _epi(hla_allele="HLA-DRB1*04:01")
    assert e.mhc_class == "II"

def test_class_ii_prose_cue_ok():
    e = _epi(quoted_text="class II-restricted CD4 epitope")
    assert e.mhc_class == "II"

def test_class_ii_no_anchor_rejected():
    with pytest.raises(ValueError, match="restriction anchor"):
        _epi(quoted_text="predicted minimal epitope")

def test_class_ii_inferred_flag_exempts_anchor():
    # explicit audited heuristic call (e.g. a length-convention table adapter): no anchor needed,
    # but it MUST be flagged for review.
    e = _epi(quoted_text="MHC class inferred from length", mhc_class_inferred=True, needs_review=True)
    assert e.mhc_class == "II" and e.mhc_class_inferred is True

def test_inferred_flag_requires_needs_review():
    # the escape hatch is only honest if it routes to a human
    with pytest.raises(ValueError, match="flagged for review"):
        _epi(quoted_text="MHC class inferred from length", mhc_class_inferred=True, needs_review=False)

def test_class_i_unaffected_by_anchor_rule():
    # class-I only needs its affinity slot; build a minimal valid class-I epitope
    from cancervac_packet.schema import Measurement
    e = MinimalEpitope(quoted_text="x", section_ref="T", paper_local_id="E2",
                       sequence="SIINFEKL", is_neoantigen=True, mhc_class="I",
                       predicted_affinity=Measurement(unit="unknown", raw="n/a"))
    assert e.mhc_class == "I"

@pytest.mark.parametrize("allele,ok", [
    ("I-Ab", True), ("H-2IAb", True), ("I-Ag7", True), ("H-2 IAb", True),
    ("HLA-DPA1*01:03/DPB1*02:01", True), ("HLA-DRB1*04:01", True),
])
def test_mouse_and_paired_class_ii_regex(allele, ok):
    assert bool(re.match(_MHC_PATTERN, allele)) is ok
