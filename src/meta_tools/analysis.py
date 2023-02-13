import statistics
from src.meta_tools.compute_results import *
from src.meta_tools.analysis_tools import RIGOUR, Result

# A script intended to analyse the parsed results of a TF run.
RESULTS_PATH = Path(
    'C:\\Users\\bfern\\OneDrive - University of Toronto\\Desktop\\Analysis Datasets')
DESKTOP = Path('C:\\Users\\bfern\\OneDrive - University of Toronto\\Desktop')


# "Unoptimised Long Primers - Replicates\\Results"
# "Unoptimised Long Primers - Rigour\\Unoptimised Long Primers Results"
#
#
#
#

def main():
    # === Config region ===
    # The range of parameters used to generate and/or analyse results files.
    params = {RIGOUR: (-10, 10)}
    # Directory containing results files.
    dirpath = DESKTOP / 'Demo\\Results'

    # One method will generate empty results files and analyse, the other will
    # just analyse. WARNING: re-creating empty results files overrides any
    # existing results files.
    # results = make_and_get_results(params, dirpath)
    # results = get_all_results(dirpath)
    results = custom_get_results(params, dirpath)

    # Do whatever you'd like with the results. See Result class for available
    # data.
    analyse(results, params)


def analyse(raw_results: List[Result], params: Dict[str, Tuple[int, int]]):
    for result in raw_results:
        result.prune()
    for result in raw_results:
        print(get_average_str(result))


def custom_get_results(params: Dict[str, Tuple[int, int]], dirpath: Path) \
        -> List[Result]:
    """Gets results from <dirpath> using non-standard methods."""
    results = []
    for filename in ['rig-10.txt', 'rig10.txt']:
        file_path = dirpath / filename
        file_params = parse_parameters(filename)
        file_homodimers, file_heterodimers = parse_tf_output(file_path)
        results.append(Result(file_params, file_homodimers, file_heterodimers))
    return results


def analyse_by_param(raw_results: List[Result], param: str, V: bool = False
                     ) -> None:
    """Prints results by rigour. Assumes raw_results contains one result for
    each rigour increment. Will print full content of every result iff V,
    otherwise prints summary."""
    for result in raw_results:
        result.prune()

    param_to_result = {}
    for result in raw_results:
        param_to_result[result.get_param_val(param)] = result

    param_vals = list(param_to_result.keys())
    param_vals.sort()

    for param_val in param_vals:
        result = param_to_result[param_val]

        print(get_average_str(result))
        if V:
            print(str(result))


def get_average_str(result: Result, do_max: bool = False,
                    newline: bool = False) -> str:
    """Returns a string formatted as
     'parameters:value - average homodimer comp: average heterodimer comp \n'
     With values extracted from result, returns max score instead of average iff
      do_max, appends newline iff newline."""
    if do_max:
        homo_val = max(result.get_homo_scores())
        hetero_val = max(result.get_hetero_scores())
    else:
        homo_val = statistics.mean(result.get_homo_scores())
        hetero_val = statistics.mean(result.get_hetero_scores())

    rtrn_str = ''
    delim = ', '
    num_param = len(result.get_parameters())
    for param_tup in enumerate(result.get_parameters()):
        param = param_tup[1]
        if param_tup[0] == num_param - 1:
            delim = ''
        rtrn_str += ''.join([str(param.get_param()), ':',
                             str(param.get_val()), delim])

    rtrn_str += " @ {homo: .2f} | {hetero: .2f}".format(homo=homo_val,
                                                        hetero=hetero_val)
    if newline:
        rtrn_str += '\n'

    return rtrn_str

    homo_avg = statistics.mean(result.get_homo_scores())
    hetero_avg = statistics.mean(result.get_hetero_scores())

    print("{param_str}: {homo:.1f}-{hetero:.2f}".format(param_str=param_str,
                                                        homo=homo_avg,
                                                        hetero=hetero_avg))


if __name__ == "__main__":
    main()
