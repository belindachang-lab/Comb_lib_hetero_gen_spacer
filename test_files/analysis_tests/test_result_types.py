import pytest

from fixtures_and_helpers import *

ALL_DMS = [
    STANDARD_DM,
    RAND_DM
]


@pytest.mark.parametrize('dimer_man', ALL_DMS)
def test_get_max_homo_unique(dimer_man) -> None:
    """Test that <dimer_man> successfully isolates the most stable of its
    dimers."""
    dimer_sets = (dimer_man.get_for_max_homo(),
                  dimer_man.get_rev_max_homo(),
                  dimer_man.get_hetero_max())
    for hset in dimer_sets:
        for i, dimer in enumerate(hset):
            for other_dimer in hset[i + 1:]:
                assert not dimer.same_primers(other_dimer)


@pytest.mark.parametrize('dimer_man', ALL_DMS)
def test_get_max_homo_correct_num(dimer_man) -> None:
    """Test that <dimer_manager> successfully isolates the most stable of its
    dimers."""
    assert len(dimer_man.get_for_max_homo()) <= 4 and \
                  len(dimer_man.get_rev_max_homo()) <= 4
    assert len(dimer_man.get_hetero_max()) <= 16


@pytest.mark.parametrize('dimer_man', ALL_DMS)
def test_get_max_homo_most_comp(dimer_man) -> None:
    """Test that <dimer_manager> successfully isolates the most stable of its
    dimers."""


