from typing import List
import copy
from statistics import mean
from src.meta_tools.analysis_tools import Heterodimer, Homodimer


class DimerManager:
    """
    A class designed to provide insight into the attributes of some sets of homo
    and heterodimers.

    === Private Attributes ===
    _for_homos:
        Homodimers formed by some forward primer in some primerset.
    _rev_homos:
        Homodimers formed by some reverse primer in some primerset.
    _heteros:
        Heterodimers formed by some primers in some primerset.
    """
    _for_homos: List[Homodimer]
    _rev_homos: List[Homodimer]
    _heteros: List[Heterodimer]

    def __init__(self, for_homos: List[Homodimer], rev_homos: List[Homodimer],
                 heteros: List[Heterodimer]) -> None:
        """Creates a dimer manager with the given values."""
        self._for_homos = for_homos
        self._rev_homos = rev_homos
        self._heteros = heteros

    def get_for_dimers_cpy(self) -> List[Homodimer]:
        """Returns a list of all the currently stored forward homodimers."""
        return copy.deepcopy(self._for_homos)

    def get_rev_dimers_cpy(self) -> List[Homodimer]:
        """Returns a list of all the currently stored reverse homodimers."""
        return copy.deepcopy(self._rev_homos)

    def get_hetero_dimers_cpy(self) -> List[Heterodimer]:
        """Returns a list of all the currently stored heterodimers."""
        return copy.deepcopy(self._heteros)

    def get_mean_for_comp(self) -> float:
        """Returns the mean complimentary amoung the currently stored forward
        homodimers."""
        for_homo_scores = [homo.get_num_comp() for homo in self._for_homos]
        return mean(for_homo_scores)
    def get_mean_rev_comp(self) -> float:
        """Returns the mean complimentary amoung the currently stored reverse
        homodimers."""

    def get_mean_hetero_comp(self) -> float:
        """Returns the mean complimentary amoung the currently stored
        heterodimers."""

    def __str__(self) -> str:
        """Returns a string representation of this result."""
        out_str = ''
        out_str += "Forward Homodimers:\n"
        for homodimer in self._for_homos:
            out_str += '\t' + str(homodimer) + '\n'

        out_str += "Reverse Homodimers:\n"
        for homodimer in self._rev_homos:
            out_str += '\t' + str(homodimer) + '\n'

        out_str += "Heterodimers:\n"
        for heterodimer in self._heteros:
            out_str += '\t' + str(heterodimer) + '\n'

        return out_str

    def get_for_max_homo(self) -> List[Homodimer]:
        """Returns the most stable forward homodimers for each conformation."""
        to_remove = []
        max_homo = copy.deepcopy(self._for_homos)
        # For each unique combinations of homodimers, see if they specify the
        # same template primer.
        for i in range(len(max_homo)):
            if i in to_remove:
                # We've marked this as having less complimentarity.
                continue
            for j in range(i + 1, len(max_homo)):
                h1 = max_homo[i]
                h2 = max_homo[j]
                same_primer = h1.get_primer_ind() == h2.get_primer_ind()
                # These primers are the same, so remove the least stable.
                if same_primer:
                    if h1.get_num_comp() < h2.get_num_comp():
                        to_remove.append(i)
                    else:
                        to_remove.append(j)
        # Ensure uniqueness.
        to_remove = list(set(to_remove))
        for ind in to_remove[::-1]:
            max_homo.pop(ind)

        return max_homo

    def get_rev_max_homo(self) -> List[Homodimer]:
        """Returns the most stable reverse homodimers for each conformation."""
        to_remove = []
        max_homo = copy.deepcopy(self._rev_homos)
        # For each unique combinations of homodimers, see if they specify the
        # same template primer.
        for i in range(len(max_homo)):
            if i in to_remove:
                # We've marked this as having less complimentarity.
                continue
            for j in range(i + 1, len(max_homo)):
                h1 = max_homo[i]
                h2 = max_homo[j]
                same_primer = h1.get_primer_ind() == h2.get_primer_ind()
                # These primers are the same, so remove the least stable.
                if same_primer:
                    if h1.get_num_comp() < h2.get_num_comp():
                        to_remove.append(i)
                    else:
                        to_remove.append(j)
        # Ensure uniqueness.
        to_remove = list(set(to_remove))
        for ind in to_remove[::-1]:
            max_homo.pop(ind)

        return max_homo

    def get_hetero_max(self) -> List[Heterodimer]:
        """Returns the most stable hetero for each conformation."""
        to_remove = []
        max_heteros = copy.deepcopy(self._heteros)
        # Repeat for heterodimers.
        for i in range(len(max_heteros)):
            if i in to_remove:
                # We've marked this as having less complimentarity.
                continue
            for j in range(i + 1, len(max_heteros)):
                h1 = max_heteros[i]
                h2 = max_heteros[j]

                same_primer = h1.get_forward_ind() == h2.get_forward_ind() and \
                              h1.get_reverse_ind() == h2.get_reverse_ind()
                # These primers are the same, so remove the least stable.
                if same_primer:
                    if h1.get_num_comp() < h2.get_num_comp():
                        to_remove.append(i)
                    else:
                        to_remove.append(j)
        # Ensure uniqueness.
        to_remove = list(set(to_remove))
        for ind in to_remove[::-1]:
            max_heteros.pop(ind)
        return max_heteros



    def prune_mismatch_heteros(self) -> None:
        """Removes heterodimers that don't have matching forward and reverse
         indices."""
        pass


class TFResult(DimerManager):
    pass



class InternalResult(DimerManager):
    pass
