from typing import List

class OligoDimer:
    """A linear dimer structure formed between two oligonucleotides"""
    pass
    _reverse_offset: int
    _num_comp: int
    _con_comp: List[int]
    _max_con_comp: int
    def __init__(self, reverse_offset: int, num_comp: int,
                 continuous_comps: List[int]) -> None:
        """Initialises an empty Homodimer object."""
        self._reverse_offset = reverse_offset
        self._num_comp = num_comp
        self._con_comp = continuous_comps[:]

    # Getters and setters

    def set_reverse_offset(self, offset: int) -> None:
        self._reverse_offset = offset
        return

    def set_num_comp(self, num_comp: int) -> None:
        self._num_comp = num_comp
        return

    #


class Homodimer(object):
    """An object used to store the homodimer results of a TF analysis.

    === Private Attributes ===:
        _primer_ind:
            The index of the input primer that produced the homodimer.

        _is_forward:
            The direction of the primer

        _reverse_offset:
            By reversing this primer and shifting <reverse_offset>
            bases to the right, this particular homodimer is produced.

        _num_comp:
            The number of complementary bases in this homodimer.

        _con_comp:
            A set of integers representing sets of adjacent complimentary bases.


        Example Usage:
        >R3
        5-gtgactggagttcagacgtgtgctcttccgatctatcagagccctgatgtgtgaagatccaccaa->
             |  | | ||     |    |   ||||| || |||||   |    |     || | |  |
           <-acaccgaaccacctagaagtgtgtagtcccgagactatctagccttctcgtgtgcagacttgagg-5

        _primer_ind = 3
        _is_forward = False
        _reverse_offset = 3
        _num_comp = 26
        _con_comp = [4, 4]
    """

    _primer_ind: int
    _is_forward: bool
    _reverse_offset: int
    _num_comp: int
    _con_comp: List[int]

    def __init__(self, primer_ind: int, is_forward: bool, reverse_offset: int,
                 num_comp: int, continuous_comps: List[int]) -> None:
        """Initialises an empty Homodimer object."""
        self._primer_ind = primer_ind
        self._is_forward = is_forward
        self._reverse_offset = reverse_offset
        self._num_comp = num_comp
        self._con_comp = continuous_comps

    def __str__(self):
        """Returns a string representation of this heterodimer"""
        return "Homodimer {dir}{num: d}:Comp: {comp: d}, Offset: {off: d}". \
            format(
            dir=['R', 'F'][self._is_forward], num=self._primer_ind,
            comp=self._num_comp, off=self._reverse_offset)

    def same_primers(self, other) -> bool:
        """Returns whether <self> and <other> refer to the same pairing of
        primers."""
        return self._primer_ind == other._primer_ind and \
               self._is_forward == other._is_forward

        # Getters and setters.

    def set_primer_ind(self, ind: int) -> None:
        self._primer_ind = ind
        return

    def set_is_forward(self, is_forward: bool) -> None:
        self._is_forward = is_forward
        return

    def set_reverse_offset(self, offset: int) -> None:
        self._reverse_offset = offset
        return

    def set_num_comp(self, num_comp: int) -> None:
        self._num_comp = num_comp
        return

    def get_primer_ind(self) -> int:
        return self._primer_ind

    def get_is_forward(self) -> bool:
        return self._is_forward


    def get_copy(self):
        """Returns a deep copy of this object"""
        hd = Homodimer(self._primer_ind, self._is_forward, self._reverse_offset,
                       self._num_comp, self._con_comp)
        return hd


class Heterodimer:
    """An object used to store the heterodimer results of a TF analysis.

    === Private Attributes ===:
        _forward_ind: The index of the input forward primer.
        _reverse_ind : The index of the input reverse primer.

        _reverse_offset: By reversing the reverse primer (so that it is 3'-5'_
            and shifting <reverse_offset> bases to the right, this particular
            heterodimer is produced.

        _num_comp: The number of complementary bases in this homodimer.
    """
    _forward_ind: int
    _reverse_ind: int
    _reverse_offset: int
    _num_comp: int

    def __init__(self, forward_ind: int, reverse_ind: int, reverse_offset: int,
                 num_comp: int) -> None:
        """Initialises an empty Heterodimer object."""
        self._forward_ind = forward_ind
        self._reverse_ind = reverse_ind
        self._reverse_offset = reverse_offset
        self._num_comp = num_comp

    def __str__(self):
        """Returns a string representation of this heterodimer"""
        return "Heterodimer F{F}-R{R}. Comp: {comp: d}, Offset: {off: d}".format(
            F=self._forward_ind, R=self._reverse_ind, comp=self._num_comp,
            off=self._reverse_offset)

    def same_primers(self, other) -> bool:
        """Returns whether <self> and <other> refer to the same pairing of
        primers."""
        return self._forward_ind == other._forward_ind and \
               self._reverse_ind == other._reverse_ind

    # Getters and setters.
    def set_forward_ind(self, ind: int) -> None:
        self._forward_ind = ind
        return

    def set_reverse_ind(self, ind: int) -> None:
        self._reverse_ind = ind
        return

    def set_reverse_offset(self, offset: int) -> None:
        self._reverse_offset = offset
        return

    def set_num_comp(self, num_comp: int) -> None:
        self._num_comp = num_comp
        return

    def get_forward_ind(self) -> int:
        return self._forward_ind

    def get_reverse_ind(self) -> int:
        return self._reverse_ind

    def get_reverse_offset(self) -> int:
        return self._reverse_offset

    def get_num_comp(self) -> int:
        return self._num_comp

    def get_copy(self):
        """Returns a deep copy of this object"""
        hd = Heterodimer(self._forward_ind, self._reverse_ind,
                         self._reverse_offset, self._num_comp)
        return hd


RIGOUR = 'rig'
SETNUM = 'set'
SETSIZE = 'sze'
REPLIC = 'rep'
FT = '.txt'
DELIM = '@'


class Parameter:
    """ A class designed to specify a parameter of a set of primers to be
    generated, or their computed properties.

    === Private Attributes ===
    param:
            The type of this parameter. One of a limited set of parameters
            types.

    val:
            The value of this parameter."""

    _param: str
    _val: int

    def __init__(self, param: str, val: int) -> None:
        """Initialises the Parameter"""
        self._param = param
        self._val = val

    def set_param(self, param: str) -> None:
        self._param = param

    def set_val(self, val: int) -> None:
        self._val = val

    def get_param(self) -> str:
        return self._param

    def get_val(self) -> int:
        return self._val


class Result:
    """The results from a single thermofischer primer analysis.

    === Private Attributes ===
    _parameters: The parameters with which this run was conducted.

    _homodimers: The homodimers found in this run.

    _heterodimers: The heterodimers found in this run."""

    _parameters: List[Parameter]
    _homodimers: List[Homodimer]
    _heterodimers: List[Heterodimer]

    def __init__(self, parameters: List[Parameter], homodimers: List[Homodimer],
                 heterodimers: List[Heterodimer]) -> None:
        """Initialises a Result with the given values."""
        self._parameters = parameters
        self._homodimers = homodimers
        self._heterodimers = heterodimers

    def __str__(self) -> str:
        """Returns a string representation of this result."""
        out_str = ''
        out_str += "Homodimers:\n"
        for homodimer in self._homodimers:
            out_str += '\t' + str(homodimer) + '\n'

        out_str += "Heterodimers:\n"
        for heterodimer in self._heterodimers:
            out_str += '\t' + str(heterodimer) + '\n'

        return out_str

    def optimise_dimers(self) -> None:
        """Removes all but the most stable homodimers for each conformation."""
        to_remove = []
        # For each unique combinations of homodimers, see if they specify the
        # same template primer.
        for i in range(len(self._homodimers)):
            for j in range(i + 1, len(self._homodimers)):
                h1 = self._homodimers[i]
                h2 = self._homodimers[j]
                same_primer = h1.get_primer_ind() == h2.get_primer_ind() and \
                              h1.get_is_forward() == h2.get_is_forward()
                # These primers are the same, so remove the least stable.
                if same_primer:
                    if h1.get_num_comp() < h2.get_num_comp():
                        to_remove.append(i)
                    else:
                        to_remove.append(j)
        # Ensure uniqueness.
        to_remove = list(set(to_remove))
        for ind in to_remove[::-1]:
            self._homodimers.pop(ind)

        to_remove = []
        # Repeat for heterodimers.
        for i in range(len(self._heterodimers)):
            for j in range(i + 1, len(self._heterodimers)):
                h1 = self._heterodimers[i]
                h2 = self._heterodimers[j]

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
            self._heterodimers.pop(ind)

    def prune_mismatch_heteros(self) -> None:
        """Removes heterodimers that don't have matching forward and reverse
         indices."""
        to_remove = []
        for i in enumerate(self._heterodimers):
            hd = i[1]
            ind = i[0]
            if hd.get_forward_ind() != hd.get_reverse_ind():
                to_remove.append(ind)

        to_remove.sort(reverse=True)
        for ind in to_remove:
            self._heterodimers.pop(ind)

    def prune(self) -> None:
        """Removes superfluous dimers from this set."""
        self.prune_mismatch_heteros()
        self.optimise_dimers()

    def get_homodimers(self) -> List[Homodimer]:
        return self._homodimers

    def get_heterodimers(self) -> List[Heterodimer]:
        return self._heterodimers

    def get_parameters(self) -> List[Parameter]:
        return self._parameters

    def get_homo_scores(self) -> List[int]:
        """Returns the scores for the homodimers found in this run."""
        homo_scores = []
        for hd in self._homodimers:
            homo_scores.append(hd.get_num_comp())
        return homo_scores

    def get_hetero_scores(self) -> List[int]:
        """Returns the scores for the heterodimers found in this run."""
        hetero_scores = []
        for hd in self._heterodimers:
            hetero_scores.append(hd.get_num_comp())
        return hetero_scores

    def get_param_val(self, param_name: str) -> int:
        for param in self._parameters:
            if param.get_param() == param_name:
                return param.get_val()
        raise ValueError("Param does not exist, or is not a parameter of this "
                         "result.")
