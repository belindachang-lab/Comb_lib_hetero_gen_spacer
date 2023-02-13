from pathlib import Path
from typing import Tuple, List

from src.meta_tools.analysis_tools import Heterodimer, Homodimer


def SeqMisalignedError(BaseException):
    """Thrown when the parser becomes misaligned with the sequence."""
    pass


def parse_tf_output(filepath: Path or str) -> Tuple[List[Homodimer],
                                                    List[Heterodimer]]:
    # Open file containing output from ThermoFisher multiple primer analyser.
    output_path = Path(Path(filepath))
    with open(output_path, 'r') as output_file:
        output_lines = output_file.readlines()

    # Will contain parsed homo/heterodimers
    homodimers = []
    heterodimers = []

    i = 0

    def cur_line() -> str:
        # Returns the current line according to i with newline removed.
        return output_lines[i][0:-1]

    mode = "und"  # undefined
    while i < len(output_lines):

        # We've entered the self dimer region of the file.
        if "Self-Dimers:" in cur_line():
            mode = 'sd'  # self dimer

        if "Cross Primer Dimers:" in cur_line():
            mode = "hd"  # heterodimers

        # Parse homodimers.
        while mode == 'sd' and i < len(output_lines):

            if 'dimers for' in cur_line():  # we're on a line specifying some homodimers for a primer.

                # Collect attributes for this set of dimers.
                num_homodimers = int(cur_line()[0])
                primer_ind = int(cur_line()[-1])
                is_forward = cur_line()[-2] == 'F'
                i += 1
                # Insert known values. We'll fix those 0's later.
                template_dimer = Homodimer(primer_ind, is_forward, 0, 0)

                # Parse each homodimer below.
                for _ in range(num_homodimers):
                    # We know that there's 4 lines in each homodimer.
                    homodimer_lines = output_lines[i: i + 4]
                    i += 4

                    # Quantify offset and complementarity.
                    comp = homodimer_lines[1].count('|')
                    offset = - homodimer_lines[0].count(' ') + \
                             homodimer_lines[2].count(' ')

                    # Catch misalignment
                    if homodimer_lines[3] != '\n':
                        raise SeqMisalignedError(
                            "Sequence misaligned at line {ind: d}".format(
                                ind=i))

                    template_dimer.set_num_comp(comp)
                    template_dimer.set_reverse_offset(offset)
                    homodimers.append(template_dimer.get_copy())
            else:
                i += 1

            if "Cross Primer Dimers:" in cur_line():
                mode = "hd"  # heterodimers

        while mode == 'hd' and i < len(output_lines):

            # At the start of a block specifying a
            if 'with' in cur_line():

                p1, p2 = cur_line()[:2], cur_line()[-2:]
                # Invert the reverse offset if necessary.
                invert = p1[0] == 'R'
                # Skip line showing leading primer. p1 is leading primer.
                i += 2

                # Find inds of forward and reverse primers.
                if p1[0] == 'F':
                    f_ind = int(p1[1])
                    r_ind = int(p2[1])
                else:
                    r_ind = int(p1[1])
                    f_ind = int(p2[1])

                # Format template dimer
                template_dimer = Heterodimer(f_ind, r_ind, 0, 0)

                # Only parse  the heterodimers for this pair of primers. Stop
                # if we've reached the end of the file or we're in the next
                # set of primers. We aren't going to go beyond the end of the
                # file.
                while i + 3 < len(output_lines) and not "with" in output_lines[
                    i]:
                    # Read lines containing heterodimer data.
                    heterodimer_lines = output_lines[i:i + 4]
                    i += 4
                    if invert:
                        # Calculate size of region beyond the 3' end of the rev
                        # primer.
                        reverse_offset = len(heterodimer_lines[2]) - \
                                         len(heterodimer_lines[0])
                    else:
                        reverse_offset = heterodimer_lines[2].count(' ')

                    comp = heterodimer_lines[1].count('|')
                    # Catch misalignment
                    if heterodimer_lines[3] != '\n':
                        raise SeqMisalignedError(
                            "Sequence misaligned at line {ind: d}".format(
                                ind=i))

                    # Add dimer to list of dimers.
                    template_dimer.set_reverse_offset(reverse_offset)
                    template_dimer.set_num_comp(comp)
                    heterodimers.append(template_dimer.get_copy())
            else:
                i += 1
        i += 1
    return homodimers, heterodimers
