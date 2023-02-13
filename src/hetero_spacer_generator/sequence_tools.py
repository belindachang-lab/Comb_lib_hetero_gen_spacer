from Bio.Seq import Seq
from typing import List, Tuple, Iterable, Callable
from src.hetero_spacer_generator.defaults import DEGEN_TO_POSSIBLE

COMP = {
    'A':'T',
    'T':'A',
    'C':'G',
    'G':'C'
}
def get_complement(seq: Seq) -> Seq:

    new_seq = ''
    for base in seq:
        new_seq += COMP[base]

    return Seq(new_seq)


def order_seqs(seq1: Seq, seq2: Seq) -> Tuple[Seq, Seq, bool]:
    """Returns a tuple of two seqs where len(tuple[0]) >= len(tuple[1]).
    Returns true iff the sequences were inverted."""
    if len(seq1) >= len(seq2):
        return seq1, seq2, False
    else:
        return seq2, seq1, True


def is_degen_base(base: str) -> bool:
    """Returns whether the base is degenerate."""
    if base in "ATCG":
        return False
    return True


def is_degen(seq: Seq) -> bool:
    """Returns whether <seq> contains degenerate bases."""
    for base in seq:
        if is_degen_base(base):
            return True
    return False


def remove_degen(seq: Seq) -> Seq:
    """Returns a modified version of <seq> that doesn't contain degenerate
    bases."""
    seq_str = ''
    for base in seq:
        if is_degen_base(base):
            seq_str += 'I'
        else:
            seq_str += base

    return Seq(seq_str)


def comp_seqs_any_overlap(seq1: Seq, seq2: Seq,
                          comp: Callable[[Seq, Seq], int]) -> None:
    """Compares seq1 and seq2 over every possible alignment of the sequences
    using comp. Returns the max val produced by comp."""
    max_val = 0
    l_seq, s_seq, flipped = order_seqs(seq1, seq2)
    n = len(s_seq)
    m = len(l_seq)
    # Where l_seq[i + s] aligns with s_seq[i] for all valid i. Iterate over all
    # s (short for shift) that allow indices valid with the above.
    for s in range(1 - n, 0):
        size = n + s
        max_val = max(comp(s_seq[-size:], l_seq[:size]), max_val)
    for s in range(0, m - n):
        size = n
        max_val = max(comp(s_seq, m[s:s + size]), max_val)
    for s in range(m - n, m - 1):
        size = m + s
        max_val = max(comp(s_seq[:size], m[s:]), max_val)


def iterate_seqs(oligo: Seq, seqs: Iterable[Seq],
                 comp_method: Callable[[Seq, Seq, int], int]) -> int:
    """Compares all completely overlapping combinations of <oligo> over each of
     <seq> using <comp_method>, returns the greatest complementarity found."""
    max_comp = 0
    oligo_len = len(oligo)
    for seq in seqs:
        # Each possible alignment along seq
        seq_len = len(seq)
        if seq_len >= oligo_len:  # Iterate over alignments of oligo over seq.
            larger_seq = seq
            smaller_seq = oligo
            num_alignments = len(seq) - oligo_len
        else:  # Iterate over alignments of seq over oligo.
            larger_seq = oligo
            smaller_seq = seq
            num_alignments = oligo_len - len(seq)

        for offset in range(0, num_alignments):
            local_max_comp = comp_method(smaller_seq, larger_seq, offset)
            if local_max_comp > max_comp:
                max_comp = local_max_comp
    return max_comp


def get_max_complementarity(oligo: Seq, seqs: Iterable[Seq]) -> int:
    """Returns the highest number of complimentary bases between <oligo> at any
    possible alignment on any of <seqs>, while <oligo> and <seqs> are still
    fully overlapped. Considers the two seq to be in opposing directions.

    >>> get_max_complementarity(Seq("ATCGA"), [Seq("ATTCAGACC")])
     4"""
    return iterate_seqs(oligo, seqs, get_site_complementarity)


def get_max_complementarity_consec(oligo: Seq, seqs: Iterable[Seq]) -> int:
    """Returns the highest number of consicutive complimentary bases between
    <oligo> at any possible alignment on any of <seqs>, while <oligo> and <seqs>
    are still fully overlapped. Considers the two seq to be in opposing
    directions.

    >>> get_max_complementarity(Seq("ATCGA"), [Seq("ATTCAGACC")])
     4"""
    return iterate_seqs(oligo, seqs, get_site_complementarity_consec)


def get_site_complementarity_consec(oligo: Seq, seq: Seq, offset: int) -> int:
    """Gets the number of complementary bases of <oligo> at the site along
    <seq> specified by  the <offset> where oligo and seq are in opposite
    directions.

    Precondition:
    offset + len(oligo) <= len(seq)
    >>> get_site_complementarity(Seq("ATCG"), Seq("TTAG CC"), 0)
     1
    >>> get_site_complementarity(Seq("ATCG"), Seq("T TAGC C"), 1)
     4
     """
    max_comp = 0
    run_comp = 0
    oligo_compliment = oligo.complement()

    # At each base, check for equality between the oligo's compliment and seq
    for i in range(len(oligo)):
        if compare_bases_degenerate(oligo_compliment[i], seq[i + offset]):
            run_comp += 1
        else:
            if max_comp < run_comp:
                max_comp = run_comp
            run_comp = 0
    if max_comp < run_comp:
        max_comp = run_comp
    return max_comp


def get_site_complementarity(oligo: Seq, seq: Seq, offset: int) -> int:
    """Gets the number of complementary bases of <oligo> at the site along
    <seq> specified by  the <offset>.

    Precondition:
    offset + len(oligo) <= len(seq)
    >>> get_site_complementarity(Seq("ATCG"), Seq("TTAG CC"), 0)
     1
    >>> get_site_complementarity(Seq("ATCG"), Seq("T TAGC C"), 1)
     4
     """
    num_complement = 0
    oligo_compliment = oligo.complement()

    # At each base, check for equality between the oligo's compliment and seq
    for i in range(len(oligo)):
        if compare_bases_degenerate(oligo_compliment[i], seq[i + offset]):
            num_complement += 1
    return num_complement


def compare_bases(b1: str, b2: str) -> bool:
    """Returns whether the two bases are complementary."""
    return b1 == b2


def compare_bases_can_ignore(b1: str, b2: str) -> bool:
    """Returns whether the two bases are complementary."""
    if b1 == 'I' or b2 == "I":
        return False
    return b1 == b2


def compare_bases_degenerate(b1: str, b2: str) -> bool:
    """Returns whether the two bases are complementary."""
    if b1 == b2:
        return True

    b1 = DEGEN_TO_POSSIBLE[b1]
    b2 = DEGEN_TO_POSSIBLE[b2]

    if b1 == 'N' or b2 == 'N':
        return True

    for p_base_1 in b1:
        for p_base_2 in b2:
            if p_base_2 == p_base_1:
                return True

    return False


class SeqAnalyzer:
    """Responsible for analysing features of sequences."""

    # Whether the degenerate bases have been specified as allowed or disallowed.
    _degeneracy_specified: bool
    # Whether to expect degeneracy in input sequences. False if not specified.
    _expect_degeneracy: bool
    # Method used to compare individual bases.
    _base_comp_method: Callable[[str, str], bool]
    # Whether this analysed should count the number of times its core comparison
    # method has been called.
    _counting: bool
    # Number of times this function's core comparison methods have been called.
    _count: int

    def __init__(self, degen: bool = None) -> None:
        """Initialises this SeqAnalyser. If degeneracy is to be expected,
        runtime will drastically increase."""
        if degen is not None:
            self.expect_degeneracy(degen)
        else:
            self._degeneracy_specified = False
            self._expect_degeneracy = False
        self._counting = False
        self._count = 0

    def start_counting_comparisons(self) -> None:
        """Begins counting the number of sequence analyses performed."""
        self._count = 0
        self._counting = True

    def stop_counting_comparisons(self) -> int:
        """Returns the number of sequence analyses performed since
        start_counting_comparisons was last called."""
        self._counting = True
        return self._count

    def expect_degeneracy(self, expect_degeneracy: bool) -> None:
        """Specifies whether degenerate bases should be expected. It is
        recommended to set this to false as it helps runtime."""
        self._expect_degeneracy = expect_degeneracy
        self._degeneracy_specified = True
        self.set_base_comp_method()

    def set_base_comp_method(self) -> None:
        """Sets the method used to compare bases."""
        if self._degeneracy_specified:
            if self._expect_degeneracy:
                self._base_comp_method = compare_bases_degenerate
            else:
                self._base_comp_method = compare_bases

    def degen_check(self, seqs: List[Seq]) -> None:
        """Makes sure that the correct base comparison method is set."""
        if self._degeneracy_specified:
            pass
        else:
            degen = False
            for seq in seqs:
                if is_degen(seq):
                    degen = True
                    break
            if degen:
                self._base_comp_method = compare_bases_degenerate
            else:
                self._base_comp_method = compare_bases

    def comp_seqs_any_overlap(self, seq1: Seq, seq2: Seq,
                              comp: Callable[[str, str], int],
                              f_skip: int = 0, r_skip: int = 0) -> int:
        """Compares <seq1> and <seq2> over every possible alignment of the
        sequences using comp. Returns the max val produced by comp.
        <f_skip> is the number of alignments to skip from

        Note:
            For proper binding values, ensure that all input sequences are
            5'-3', but the function will work if this is
            switched."""

        self.degen_check([seq1, seq2])

        if self._counting:
            self._count += 1

        max_val = 0
        l_seq, s_seq, flipped = order_seqs(seq1, seq2)
        n = len(s_seq)
        m = len(l_seq)
        s_seq = get_complement(s_seq)
        if flipped:
            f_skip, r_skip = r_skip, f_skip
        # Where l_seq[i + s] aligns with s_seq[i] for all valid i. Iterate over all
        # shift that allow indices valid with the above.

        # Convert to strings
        l_seq = l_seq.__str__()
        # Convert second seq to 3'-5'.
        s_seq = s_seq.__str__()[::-1]

        # Check combinations with full overlap first.
        for shift in range(0, m - n):
            size = n
            max_val = max(comp(s_seq, l_seq[shift:shift + size]), max_val)


        # Check combinations with partial overlap. There's no point checking
        # overlaps with size less than max_val, so skip them

        # 5' of s_seq overlaps 5' of l_seq
        for shift in range(0, 1 - n + f_skip, -1):
            size = n + shift
            if max_val > size:
                continue
            max_val = max(comp(s_seq[-size:], l_seq[:size]), max_val)

        # 3' of s_seq overlaps 3' of l_seq
        for shift in range(m - n, m - 1 - r_skip):
            size = m - shift
            if max_val > size:
                continue
            max_val = max(comp(s_seq[:size], l_seq[shift:]), max_val)

        return max_val

    def get_consec_complementarity(self, seq1: str, seq2: str) -> int:
        """Returns the number of consecutive complementary bases between <seq1>
        and <seq2>.

        Assumes than len(seq1) == len(seq2)."""
        max_consec = 0
        running_consec = 0
        length = len(seq1)

        for i in range(length):
            if self._base_comp_method(seq1[i], seq2[i]):
                running_consec += 1
            else:
                if max_consec < running_consec:
                    max_consec = running_consec
                running_consec = 0

        if max_consec < running_consec:
            return running_consec
        return max_consec

    def get_non_consec_complementarity(self, seq1: str, seq2: str) -> int:
        """Returns the number of complementary between <seq1> and <seq2>.

        Assumes than len(seq1) == len(seq2)."""

        num_comp = 0
        for i in range(len(seq1)):
            if self._base_comp_method(seq1[i], seq2[i]):
                num_comp += 1
        return num_comp


