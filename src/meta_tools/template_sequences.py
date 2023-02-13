from Bio.Seq import Seq
from typing import List, Tuple, Union


def mix_to_seq(seqs_or_strs: List[Union[Seq or str]]) -> None:
    """Converts any str value in <seq_to_strs> into a Seq."""
    for i, seq in enumerate(seqs_or_strs):
        if isinstance(seq, str):
            seqs_or_strs[i] = Seq(seq)


class GenSet:
    """Stores primers to be used in testing of primer generation algorithms.

    All attributes are public.
    Stores some attributes used for the construction of primer sets.

    """
    name: str
    _for_bindings: List[Seq]
    _rev_bindings: List[Seq]
    _for_adapters: List[Seq]
    _rev_adapters: List[Seq]

    def __init__(self, name: str, for_bindings: List[Union[Seq or str]],
                 rev_bindings: List[Union[Seq or str]],
                 for_adapters: List[Union[Seq or str]],
                 rev_adapters: List[Union[Seq or str]]) -> None:
        """Initialises and empty GenSet"""
        self._for_bindings = for_bindings[:]
        self._rev_bindings = rev_bindings[:]
        self._for_adapters = for_adapters[:]
        self._rev_adapters = rev_adapters[:]
        mix_to_seq(self._for_bindings)
        mix_to_seq(self._rev_bindings)
        mix_to_seq(self._for_adapters)
        mix_to_seq(self._rev_adapters)

        if not len(self._for_bindings) == len(self._rev_bindings) and \
                len(self._for_adapters) == len(self._rev_adapters):
            raise ValueError("Given lists of sequences must have the"
                             "same length.")
        self.name = name
        return

    def num_bind(self) -> int:
        return len(self._rev_bindings)

    def num_adapt(self) -> int:
        return len(self._rev_adapters)

    def unpack(self) -> Tuple[List[Seq], List[Seq], List[Seq], List[Seq]]:
        return self._for_bindings, self._rev_bindings, \
               self._for_adapters, self._rev_adapters

    def has_binding_num(self, ind: int) -> bool:
        return min(len(self._for_bindings), len(self._rev_bindings)) < ind




fb = [Seq('ACACTCTTTCCCTACACGACGCTCTTCCGATCT')]
rb = [Seq('GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT')]
fa = [Seq('GTCAACCCAGCAGCTTACGCTGCCCTGGGTGCCGC')]
ra = [Seq('TGTGGCTTGGTGGATCTTCACACATCAGGGCTCTGAT')]

STANDARD_SET = GenSet('Standard', fb, rb, fa, ra)

fb = [Seq('GTCAACCCAGCAGCTTACGCTGCCCTGGGTGCCTACA'),
      Seq('TCACCACAACGATGTACACCTCCATGCAC'),
      Seq('CTTGGCTGCAATCTGGAAGGATTCTTTGC'),
      Seq('TCTCTGGTCACTGGTCGTTCTGGCTATTG'),
      Seq('AGCCCATCAGCAACTTCCGCTTCGGAGAG'),
      Seq('AACTACATCCTGCTGAATCTCGCGGTGGCCGACC')]
rb = [Seq('TGTGGCTTGGTGGATCTTCACACATCAGGGCTCTGAT').reverse_complement(),
      Seq('TTTGCCAAGAGTTCCTCCATCTACAACCC').reverse_complement(),
      Seq('AGACCACCCAGAGGGCTGAGAGGGAAG').reverse_complement(),
      Seq('ATCATCATGGTCTTCGCCTTCCTGATG').reverse_complement(),
      Seq('CCGAGACCACCCAGAGGGCTGAGAGGGAA').reverse_complement(),
      Seq('AATTGTGTTCTTCTGCTACGGCCGTCTGCTGTGCGCC').reverse_complement()]

# Reuse adapters.
RANDOM_RHO_SET = GenSet('Random Rho', fb, rb, fa, ra)

fa = [Seq('ACACT')]
ra = [Seq('GTGAC')]

fb = [Seq('GTCAA'),
      Seq('TCACC'),
      Seq('CTTGG'),
      Seq('TCTCT'),
      Seq('AGCCC'),
      Seq('AACTA')]
rb = [Seq('TGTGG').reverse_complement(),
      Seq('TTTGC').reverse_complement(),
      Seq('AGACC').reverse_complement(),
      Seq('ATCAT').reverse_complement(),
      Seq('CCGAG').reverse_complement(),
      Seq('AATTG').reverse_complement()]

SMALL_BINDING_SET = GenSet('Small Binding', fb, rb, fa, ra)

fa = [Seq("TTTAGGACCGGGGAATCGGTAGGCGAACGA")]
ra = [Seq('CTAGTGCTACGGAGTAGGCATTTATATTTA')]
fb = [Seq('CGTATCCCGAATTAGCCAGTTTGGATGGGT')]
rb = [Seq('ACATCTTCCTTCACAGGGTTCCAGTCACCA')]

OPTIMAL_SET = GenSet('Optimal Set', fb, rb, fa, ra)

# SET OF PRIMERS RECEIVED FROM ALEX 01/02/2022



fa = [Seq('ACACTCTTTCCCTACACGACGCTCTTCCGATCT')]
ra = [Seq('GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT')]

fb = ["GTCATCTTCTTCTGCTACG"]
rb = ["GCCGGAATGGTCATGAAGA"]
name = "HuRho_1"
AS1 = GenSet(name, fb, rb, ra, fa)

fb = ["GAGTCCTTCGTCATCTACATG"]
rb = ["AGGGTGGTGATCATGCAG"]
name = "EnviRho_1"
AS2 = GenSet(name, fb, rb, ra, fa)

fb = ["GTCGGTAAAACTCGTGCCAGC"]
rb = ["CATAGTGGGGTATCTAATCCCAGTTTG"]
name = "Mitofish_UF"
MITOFISH = GenSet(name, fb, rb, ra, fa)

fb = ["ACACCGCCCGTCACTCT"]
rb = ["CTTCCGGTACACTTACCATG"]
name = "Teleo"
AS4 = GenSet(name, fb, rb, ra, fa)

fb = ["CTTGGTCATTTAGAGGAAGTAA"]
rb = ["GCTGCGTTCTTCATCGATGC"]
name = "ITS"
AS5 = GenSet(name, fb, rb, ra, fa)

fb = ["SGCCTGTTTACCAAAAACATCAC"]
rb = ["CTCCATAGGGTCTTCTCGTCTT"]
name = "sixteen"
AS6 = GenSet(name, fb, rb, ra, fa)

ALEX_SETS = [AS1, AS2, MITOFISH, AS4, AS5, AS6]
