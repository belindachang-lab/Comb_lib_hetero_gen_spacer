from src.meta_tools.analysis_tools import RIGOUR, SETNUM, SETSIZE
from src.meta_tools.compute_results import *


def get_complex_params_files() -> Tuple[Dict[str, Tuple[int, int]], List[str]]:
    """Returns a dictionary of parameters and their ranges and a list of
    filenames generated with those parameters."""
    parameters = {RIGOUR: (1, 3), SETNUM: (2, 4), SETSIZE: (5, 10)}
    filenames = create_empty_output_files(parameters, make_files=False)
    return parameters, filenames



# Test basic functionality of create_empty_output_files.
def test_create_empty_output_files() -> None:
    parameters = {RIGOUR: (1, 2), SETNUM: (1, 2)}
    filenames = create_empty_output_files(parameters, make_files=False)
    assert filenames == ['rig1@set1.txt', 'rig1@set2.txt',
                         'rig2@set1.txt', 'rig2@set2.txt']

# Tests the get_complex_params_files fixture, and in doing so test_files that
# create_empty_output_files can return complicated filenames.
def test_generate_complex_files() -> None:
    file_names = get_complex_params_files()[1]
    assert len(file_names) == 3 * 3 * 6

def test_compute_results() -> None:
    params, filenames = get_complex_params_files()

    params_array = []
    for filename in filenames:
        params_array.append(parse_parameters(filename))

    # Three different types in each set.
    num_types = param_per_set = 3

    # Test that all parameters are valid.
    for param_set in params_array:
        assert len(param_set) == param_per_set
        for param in param_set:
            # Assert param has valid name
            assert param.get_param() in params.keys()
            upper = params[param.get_param()][1]
            lower = params[param.get_param()][0]
            # Assert param has valid value.
            assert lower <= param.get_val() <= upper

    # Test that all parameters are unique.
    for i in range(len(params_array)):
        for j in range(i + 1, len(params_array)):
            unique = False
            set1 = params_array[i]
            set2 = params_array[j]
            obs_types = 0
            for k in range(num_types):
                for l in range(k, num_types):
                    # Are they of the same type?
                    if set1[k].get_param() == set2[l].get_param():
                        obs_types += 1
                        if set1[k].get_val() != set2[l].get_val():
                            unique = True
            assert obs_types == num_types
            assert unique

    # Test that a sufficient number of parameters have been generated.
    assert len(params_array) == 3 * 3 * 6

