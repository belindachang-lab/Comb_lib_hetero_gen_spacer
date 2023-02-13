import timeit
from pathlib import Path
from typing import List, Tuple, Union
import math
import numpy as np
from copy import deepcopy
from random import randint, choices
from statistics import stdev
import time

# Apple trees and goat cheese.
BASE_TO_INT = {
    'A': 0,
    'T': 1,
    'G': 2,
    'C': 3,
    '-': 4
}


def compute_column_score(bases: Union[List[int], np.ndarray],
                         total: int) -> int:
    """Given the base frequency of each type of base in an alignment column,
    ordered: A, T, C, G, unassigned (matching indexing with BASE_TO_INT),
    returns some score quantifying the level (or potential level) of diversity
    in that column. Total is the total number of bases in that column.

    Precondition:
        <bases> must match the formatting given, being of length 5."""

    b_inds = [0, 1, 2, 3]
    b_vals = np.array(bases[1:], dtype=np.int8)

    # Sort from least to most occurrences.
    list.sort(b_inds, key=lambda a: b_vals[a])

    unassigned = bases[0]
    while unassigned > 0:
        # Assign base to the lowest abundance base.
        b_vals[b_inds[0]] += 1
        unassigned -= 1

        # Adjust position of incremented base to make sure indices are properly
        # sorted.
        for i in range(3):
            if b_vals[b_inds[i]] > b_vals[b_inds[i + 1]]:
                b_inds[i], b_inds[i + 1] = b_inds[i + 1], b_inds[i]
            else:
                break
    # Calculate percent representation.
    percent_reps = [max(val / total * 100, 1) for val in b_vals]

    return int(np.prod(percent_reps))


class NumpyBindingAlign:
    """A class designed to enable rapid comparison of diversity at within a
    short alignment.

    _og_align:
        The original alignment of DNA sequences.

    _numpy_align:
        An 2D array representation of the alignment, with bases being converted
        to numbers.

    _spacer_sizes:
        The sizes of the heterogeneity spacers used for each of the sequences in
        the alignment. This should be an alias of some other list.

    _num_seqs:
        The number of sequences in the alignment.

    _num_heterogeneity:
        Size of the heterogeneity region."""

    _og_align: List[str]

    _numpy_align: np.ndarray

    _spacer_sizes: List[int]

    _num_seqs: int

    _num_hetero: int

    def __init__(self, align: List[str], spacer_sizes: List[int],
                 num_hetero: int) -> None:
        """Converts the alignment to a numpy array and stores it internally.

        Preconditions:
            Number of bases of heterogeneity must be less than length of the
            shortest sequence in alignment."""
        self._og_align = deepcopy(align)
        self._num_hetero = num_hetero
        self._num_seqs = len(align)

        lens = [len(seq) for seq in align]
        min_len = min(lens)
        max_len = max(lens)

        if min_len < num_hetero:
            raise ValueError("Number of bases of heterogeneity must be less "
                             "than the length of the shortest sequence in "
                             "alignment.")

        # Create integer representation of the alignment.
        int_rep_arr = []
        for i, seq in enumerate(align):
            int_rep_arr.append([])
            j = max_len
            for j, char in enumerate(seq):
                num_rep = BASE_TO_INT[char]
                int_rep_arr[i].append(num_rep)
            for k in range(j, max_len - 1):
                int_rep_arr[i].append(BASE_TO_INT['-'])

        # Convert to numpy array.
        self._numpy_align = np.array(int_rep_arr, dtype=np.int8)

        # Create spacer sizes. ALIASING!
        if len(spacer_sizes) != self._num_seqs:
            raise ValueError('Spacer list should be same length as alignment.')

        self._spacer_sizes = spacer_sizes

    def find_min_div(self) -> int:
        """Returns the lowest diversity value found in the  first
        <self._num_hetero> bases."""
        # Set to theoretical max.
        min_score = 25 ** 25
        # Find diversity score of each column.
        for col in range(self._num_hetero):
            # Number of each type of base in the current column,
            #                     A, T, C, G, unassigned
            # In accordance with BASE_TO_INT
            num_bases = np.array([0, 0, 0, 0, 0], dtype=np.int32)

            # Count number of each nucleotide in each sequence within the
            # column.
            for b in range(self._num_seqs):
                # Base lies in heterogeneity spacer.
                if col < self._spacer_sizes[b]:
                    num_bases[0] += 1
                    continue

                # Accounting for the shift introduced by heterogeneity spacers.
                shifted_col = col - self._spacer_sizes[b]
                # Increment base count.
                num_bases[self._numpy_align[b, shifted_col] + 1] += 1

            # Compute score.
            score = compute_column_score(num_bases, self._num_seqs)
            if score < min_score:
                min_score = score

        return min_score


def gen_random_seq_str(length: int, GC: float = 0.5) -> str:
    """Generates a random Seq of <length>"""
    AT = 1 - GC
    bases = ['A', 'T', 'C', 'G']

    return ''.join(choices(bases, weights=[AT / 2, AT / 2, GC / 2, GC / 2],
                           k=length))

def get_runtime(num_seqs: int, num_hetero: int,) -> float:
    """Gets the runtime for some number of find_max_div operations. Returns the
    average runtime of each operation."""

    # Inputs for a standard use case.
    len_seqs = 20
    reps = 1000
    seqs = [gen_random_seq_str(len_seqs) for _ in range(num_seqs)]
    spacer_sizes = [0 for _ in range(num_seqs)]

    aln = NumpyBindingAlign(seqs, spacer_sizes, num_hetero)

    def run() -> None:
        for i in range(len(spacer_sizes)):
            spacer_sizes[i] = randint(0, num_hetero)
        score = aln.find_min_div()

    t = timeit.timeit(lambda: run(), number=reps)
    sec_per_op = t / reps
    op_per_sec = 1 / sec_per_op
    s = ['Runtime:', str(t), '\n',
         'Sec per op:', str(sec_per_op), '\n',
         'Ops per sec:', str(op_per_sec)
         ]

    print(' '.join(s))
    return sec_per_op

def get_time_string(secs: int or float) -> str:
    """Convert the seconds value into a more readable time string."""
    h = secs // (60 ** 2)
    secs -= h * 60 ** 2
    m = secs // 60
    secs -= m * 60
    secs = int(secs)
    s = []
    if h > 0:
        s += [str(h), 'hours,']
    if m > 0:
        s += [str(m), 'minutes, and', str(secs), 'seconds']
    elif secs > 0:
        s += [str(secs), 'seconds']
    else:
        s += ['less than a second']

    return ' '.join(s)


def get_num_comp(num_seqs: int, len_hetero: int, time_per_comp: float) -> float:
    """Prints the time required to fully explore all possible heterogeneity
    spacers as described"""
    t = (len_hetero + 1) ** num_seqs * time_per_comp

    s = [
        'To fully explore an alignment of ', str(num_seqs), ' sequences with a ',
        'heterogeneity region of size ', str(len_hetero), ' it would take ',
        get_time_string(t), '.'
    ]

    print(''.join(s))


def gen_csv(num_seqs: int, num_hetero: int, num_reps: int,
            output_filepath: Path) -> None:
    """Outputs several files containing the lowest diversity scores from
    randomly generated alignments matching the above specifications."""
    len_seqs = 20

    t0 = time.time()

    print('Starting generation...')

    def incr_spacer_sizes(lst: List[int]) -> None:
        """Increments the spacer size."""
        for i, val in enumerate(lst):
            if val == num_hetero:
                lst[i] = 0
            else:
                lst[i] += 1
                break

    for r in range(num_reps):

        out = output_filepath / ('R' + str(r) + '.csv')
        file = open(out, 'w')
        seqs = [gen_random_seq_str(len_seqs) for _ in range(num_seqs)]

        with open(output_filepath / ('R' + str(r) + '.txt'), 'w') as txt:
            for seq in seqs:
                txt.write(seq + '\n')

        spacer_sizes = [0 for _ in range(num_seqs)]
        aln = NumpyBindingAlign(seqs, spacer_sizes, num_hetero)

        for _ in range((num_hetero + 1) ** num_seqs):

            me = sum(spacer_sizes) / len(spacer_sizes)
            sd = stdev(spacer_sizes)
            score = aln.find_min_div()
            file.write(', '.join([str(me), str(sd), str(score)]) + '\n')
            incr_spacer_sizes(spacer_sizes)

        file.close()

        print('Replicate #' + str(r), 'complete.')

        if r == 0:
            t1 = time.time()
            rt1 = t1 - t0
            print("Estimated remaining runtime: ",
                  get_time_string(rt1 * (num_reps - 1)))

    tf = time.time()
    print('Complete. Runtime: ', get_time_string(tf - t0))







if __name__ == '__main__':
    NUM_SEQS = 5
    NUM_HETERO = 12
    NUM_REPS = 4
    PATH = Path('C:\\Users\\bfern\\OneDrive - University of Toronto\\Desktop\\Hetero Output\\Het12 Seq5')

    tpc = get_runtime(NUM_SEQS, NUM_HETERO)
    get_num_comp(NUM_SEQS, NUM_HETERO, tpc)
    gen_csv(NUM_SEQS, NUM_HETERO, NUM_REPS, PATH)

    NUM_SEQS = 6
    NUM_HETERO = 12
    NUM_REPS = 3
    PATH = Path('C:\\Users\\bfern\\OneDrive - University of Toronto\\Desktop\\Hetero Output\\Het12 Seq6')

    tpc = get_runtime(NUM_SEQS, NUM_HETERO)
    get_num_comp(NUM_SEQS, NUM_HETERO, tpc)
    gen_csv(NUM_SEQS, NUM_HETERO, NUM_REPS, PATH)

    NUM_SEQS = 7
    NUM_HETERO = 11
    NUM_REPS = 1
    PATH = Path('C:\\Users\\bfern\\OneDrive - University of Toronto\\Desktop\\Hetero Output\\Het11 Seq7')

    tpc = get_runtime(NUM_SEQS, NUM_HETERO)
    get_num_comp(NUM_SEQS, NUM_HETERO, tpc)
    gen_csv(NUM_SEQS, NUM_HETERO, NUM_REPS, PATH)



