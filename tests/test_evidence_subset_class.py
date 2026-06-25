"""v2.12 P23: CD4/CD8 t_cell_subset + evidence mhc_class — enums, optionality, verbatim-token faithfulness."""
import pytest
from cancervac_packet.schema import (
    ExtractedEvidence, TCellSubset, EvidenceMhcClass, SCHEMA_VERSION,
)
from typing import get_args


def _ev(**kw):
    base = dict(quoted_text="x", section_ref="Results", patient_paper_id="P1",
                target_kind="epitope", epitope_paper_id="E1", assay="elispot", outcome="immunogenic")
    base.update(kw)
    return ExtractedEvidence(**base)


def test_schema_version_bumped():
    assert SCHEMA_VERSION == "2.15.0"

def test_enum_members():
    assert set(get_args(TCellSubset)) == {"cd4", "cd8", "bulk_or_unknown"}
    assert set(get_args(EvidenceMhcClass)) == {"class_i", "class_ii", "not_determined"}

def test_fields_default_none():
    e = _ev()
    assert e.t_cell_subset is None and e.mhc_class is None

def test_unknown_values_need_no_cue():
    e = _ev(t_cell_subset="bulk_or_unknown", mhc_class="not_determined")
    assert e.t_cell_subset == "bulk_or_unknown"

def test_cd8_with_cue_ok():
    e = _ev(quoted_text="CD8+ cytotoxic T-cell response by IFN-g ELISpot", t_cell_subset="cd8")
    assert e.t_cell_subset == "cd8"

def test_cd8_without_cue_rejected():
    with pytest.raises(ValueError, match="cd8"):
        _ev(quoted_text="response detected at week 16", t_cell_subset="cd8")

def test_cd4_with_cue_ok():
    e = _ev(quoted_text="CD4+ helper response to the long peptide", t_cell_subset="cd4")
    assert e.t_cell_subset == "cd4"

def test_cd4_without_cue_rejected():
    with pytest.raises(ValueError, match="cd4"):
        _ev(quoted_text="positive ELISpot", t_cell_subset="cd4")

def test_class_ii_with_cue_ok():
    e = _ev(quoted_text="HLA-DRB1*04:01-restricted response", mhc_class="class_ii")
    assert e.mhc_class == "class_ii"

def test_class_ii_without_cue_rejected():
    with pytest.raises(ValueError, match="class_ii"):
        _ev(quoted_text="reactive in ELISpot", mhc_class="class_ii")

def test_class_i_via_named_allele_field():
    # cue may live on hla_allele, not quoted_text
    e = _ev(quoted_text="strong response", hla_allele="HLA-A*02:01", mhc_class="class_i")
    assert e.mhc_class == "class_i"

def test_orthogonal_both_set():
    e = _ev(quoted_text="CD8+ cytotoxic, HLA-A*02:01-restricted", t_cell_subset="cd8", mhc_class="class_i")
    assert (e.t_cell_subset, e.mhc_class) == ("cd8", "class_i")


# ---- v2.13.1 (Fable review): murine H-2 cues accepted; sibling-provenance bleed no longer counts ----

def test_murine_h2_cd8_cue_ok():
    # a faithful mouse CD8 row cued only by an H-2 class-I allele must NOT be rejected
    e = _ev(quoted_text="H-2Kb-restricted response by IFN-g ELISpot", t_cell_subset="cd8")
    assert e.t_cell_subset == "cd8"

def test_murine_h2_cd4_cue_ok():
    e = _ev(quoted_text="I-Ab-restricted CD4 help", t_cell_subset="cd4")
    assert e.t_cell_subset == "cd4"

def test_provenance_only_cue_no_longer_counts():
    # a CD8 cue living ONLY in a sibling provenance quote (possibly about a different target) no longer
    # licenses the label — the cue must be on the row's own quoted_text/assay_detail/hla_allele
    with pytest.raises(ValueError, match="cd8"):
        _ev(quoted_text="response detected at week 12", t_cell_subset="cd8",
            provenance=[{"kind": "table", "locator": "Table 1", "quoted_text": "CD8+ cytotoxic"}])
