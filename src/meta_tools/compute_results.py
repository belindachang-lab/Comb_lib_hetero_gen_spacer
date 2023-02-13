from src.meta_tools.parse_tf_out import parse_tf_output
from src.meta_tools.analysis_tools import DELIM, FT, Parameter, Result
from typing import List, Tuple, Dict
from pathlib import Path
import os

# Parses some input file from a TF run.


# Pre-set Parameter types parameters must be alpha and of length 3.


def parse_parameters(filename: str, delim: str = DELIM,
                     file_type: str = FT) -> List[Parameter]:
    """Returns all parameters specified in a <filename>. Assumes that parameters
    are seperated by <delim>."""
    filename = filename.strip()
    received_filetype = filename[-len(file_type):]

    if file_type != received_filetype:
        raise ValueError("File must have filetype matching the one given.")

    params = []
    filename = filename[:-len(file_type)]
    raw_params = filename.split(delim)
    for raw_param in raw_params:
        param = Parameter()

        # Extract and check param type.
        param_type = raw_param[0:3]
        if len(param_type) > 3 or not param_type.isalpha():
            raise ValueError("Received invalid input parameter: " + raw_param)
        param.set_param(param_type)

        # Extract and check param value.
        param_val = raw_param[3:]
        try:
            param.set_val(int(param_val))
        except ValueError:
            raise ValueError("Received invalid input parameter: " + raw_param)
        params.append(param)

    return params


def get_all_results(dirpath: str or Path) -> List[Result]:
    """Parses a set of formatted TF results in <filepath> and stores them in a
    list of Result objects, which are then returned."""
    # Convert to path if neccesary
    if type(dirpath) == str:
        dirpath = Path(dirpath)

    filenames = os.listdir(dirpath)
    results = []
    # Compute result for each file in directory
    for filename in filenames:
        file_path = dirpath / filename
        file_params = parse_parameters(filename)
        file_homodimers, file_heterodimers = parse_tf_output(file_path)
        results.append(Result(file_params, file_homodimers, file_heterodimers))

    return results


def create_empty_output_files(parameters: Dict[str, Tuple[int, int]],
                              filepath: Path = '', filetype: str = FT,
                              make_files: bool = True) -> None or List[str]:
    """Creates empty text files in <filepath> with names given by <parameters>.
    A file name will be produced for every combination of parameters. Will not
    makes files iff not <make_files>, instead will return string of filenames.

    Precondition:
        make_files == True => <filepath> is a valid filepath.

    parameters: <str> being the name of the parameter and <Tuple[int, int]>
        being the range of values that parameter can span [int, int]. First
        index must be less than or equal to second index.


    >>> import os
    >>> filepath = Path("./Desktop") # an empty file
    >>> parameters = {RIGOUR: (1, 3), SETNUM: (1, 3)}
    >>> create_empty_output_files(parameters, filepath)
    >>>os.listdir(filepath)
    """
    # Calculate number of outfiles.
    ranges = parameters.values()
    num_outfiles = 1
    for range_inds in ranges:
        num_outfiles = num_outfiles * (range_inds[1] - range_inds[0] + 1)

    # Generate empty outfiles names.
    filenames = []
    for i in range(num_outfiles):
        filenames.append('')

    num_params = len(parameters.keys())
    delim = DELIM
    repeat_length = num_outfiles
    # Append parameters to strings.
    for param in enumerate(parameters.keys()):
        range_inds = parameters[param[1]]
        param_name = param[1]
        param_ind = param[0]
        # We don't need a delimieter on the last param.
        if param_ind == num_params - 1:
            delim = ''

        num_repeats = num_outfiles // repeat_length
        # Append to filenames in tree-like fashion.
        repeat_length = repeat_length // (range_inds[1] - range_inds[0] + 1)
        for i in range(num_repeats):
            repeat_start = i * repeat_length * \
                           (range_inds[1] - range_inds[0] + 1)
            next_repeat_ind = 0
            # Each allowable value of a parameter according to the range.
            for param_num in range(range_inds[0], range_inds[1] + 1):
                last_repeat_ind = next_repeat_ind
                next_repeat_ind += repeat_length
                # Repeat to construct tree.
                param_str = param_name + str(param_num) + delim
                for i in range(last_repeat_ind, next_repeat_ind):
                    filenames[repeat_start + i] += param_str

    # Append filetype
    for i in range(len(filenames)):
        filenames[i] += filetype

    # We needn't make the files.
    if not make_files:
        return filenames

    # Make files.
    for filename in filenames:
        path = filepath / filename
        with open(path, 'w') as new_file:
            new_file.write('\n')

    return


def make_and_get_results(params: Dict[str, Tuple[int, int]], filepath: Path) \
        -> List[Result]:
    """Promps the user to make empty result files, then waits for them to enter
    data, and returns results."""

    inpt = input("Enter 'Y' if you'd like to create empty result files in " +
                 str(filepath) + '\n' + "Enter 'S' to skip this step.\n>")
    if inpt == 'S':
        pass
    elif inpt.strip() != 'Y':
        quit()

    create_empty_output_files(params, filepath)

    while input != 'R':
        inpt = input("Enter 'R' if ready to parse data in the above path or 'Q'"
                     " to quit.")
        if inpt == 'Q':
            quit()

    return get_all_results(filepath)
