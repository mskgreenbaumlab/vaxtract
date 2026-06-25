"""`derive_mutation_specific` — the source-agnostic mutant-vs-WT specificity primitive.

Generalizes the scientifically meaningful axis (mutant-preferential reactivity) out of the
33064988-specific cross-reactivity adapter, enforcing the measured-two-arm rule in ONE tested
place. A prediction or a single arm can NEVER establish specificity -> None (refuse, don't fabricate).
Reviewed with Fable 5 (2026-06-09).
"""
import agent_core
import pytest

dms = agent_core.derive_mutation_specific


@pytest.mark.parametrize("mut,wt,measured,expected", [
    # measured two-arm comparisons -> a real claim
    (True,  False, True,  True),    # mutant reactive, WT not -> mutant-preferential
    (True,  True,  True,  False),   # both reactive -> cross-reactive (a real biological result)
    # refusals (None): never fabricate a comparative claim
    (True,  False, False, None),    # not a measured assay (in-silico) -> never evidence
    (True,  None,  True,  None),    # single-arm: WT not tested
    (None,  False, True,  None),    # mutant arm not observed
    (False, False, True,  None),    # no positive response -> specificity N/A (lives on immunogenicity axis)
    (False, True,  True,  None),    # mutant negative -> still N/A here
])
def test_contract(mut, wt, measured, expected):
    assert dms(mut, wt, measured=measured) is expected


def test_false_and_none_not_collapsed():
    # cross-reactive (False) is a measured finding; un-derivable (None) is absence of one. Distinct.
    assert dms(True, True, measured=True) is False
    assert dms(True, None, measured=True) is None


def test_33064988_adapter_mapping_unchanged():
    # the adapter feeds mutant_reactive=True (epitope is immunogenic) + wt_reactive=cross_reactive_flag.
    # No -> not cross-reactive -> mutant-specific True;  Yes -> cross-reactive -> False. (regression lock)
    assert dms(True, False, measured=True) is True    # 'cross-reactive to WT: No'
    assert dms(True, True, measured=True) is False     # 'cross-reactive to WT: Yes'
