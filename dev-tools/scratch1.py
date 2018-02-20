import argparse
import textwrap
import json

parameters_file = "/user/simon/development/soar/goodman/dev/param.json"


def load_json_params(json_file, param_group):
    with open(json_file) as json_params:
        params = json.load(json_params)
    redccd = []
    try:
        redccd_args = params['ccd_args']
        for arg in redccd_args:
            print("{:s} = {:s}".format(arg, redccd_args[arg]))
    except KeyError as key:
        print('Wrong Key {:s}'.format(key))

    redspec = []
    try:
        redspec_args = params['spec_args']
        for arg in redspec_args:
            print("{:s} = {:s}".format(arg, redspec_args[arg]))
            if redspec_args[arg] != "False":
                redspec.append(arg)
                redspec.append(redspec_args[arg])
    except KeyError as key:
        print('Wrong Key {:s}'.format(key))

    print(" ")
    print(redspec)
    print(" ")
    return redspec



arguments = load_json_params(json_file=parameters_file, param_group='spec_args')

# arguments = None
# arguments = ["--data-path", "./",
#              "--proc-path", "./",
#              "--search-pattern", "cfzsto",
#              "--output-prefix", "o",
#              "--extraction", "simple",
#              "--reference-files", "ref_comp/",
#              "--max-targets", "3"]


parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            '''Extracts goodman spectra and does wavelength calibration.'''))

parser.add_argument('--data-path',
                    action='store',
                    default='./',
                    type=str,
                    metavar='<Source Path>',
                    dest='source',
                    help='Path for location of raw data. Default <./>')

parser.add_argument('--proc-path',
                    action='store',
                    default='./',
                    type=str,
                    metavar='<Destination Path>',
                    dest='destination',
                    help='Path for destination of processed data. Default '
                         '<./>')

parser.add_argument('--search-pattern',
                    action='store',
                    default='cfzsto',
                    type=str,
                    metavar='<Search Pattern>',
                    dest='pattern',
                    help="Pattern for matching the goodman's reduced data.")

parser.add_argument('--output-prefix',
                    action='store',
                    default='g',
                    metavar='<Out Prefix>',
                    dest='output_prefix',
                    help="Prefix to add to calibrated spectrum.")

parser.add_argument('--extraction',
                    action='store',
                    default='simple',
                    type=str,
                    metavar='<Extraction Type>',
                    dest='extraction_type',
                    choices=['simple', 'optimal'],
                    help='Choose a which extraction to perform. Simple is a '
                         'sum across the spatial direction after the '
                         'background has been removed. Optimal is a more '
                         'advanced method that considers weights and profile'
                         'fitting.')

parser.add_argument('--reference-files',
                    action='store',
                    default='ref_comp/',
                    metavar='<Reference Dir>',
                    dest='reference_dir',
                    help="Directory of Reference files location")

parser.add_argument('--interactive',
                    action='store_true',
                    dest='interactive_ws',
                    help="Interactive wavelength solution."
                         "Disbled by default.")

parser.add_argument('--debug',
                    action='store_true',
                    dest='debug_mode',
                    help="Debugging Mode")

parser.add_argument('--log-to-file',
                    action='store_true',
                    dest='log_to_file',
                    help="Write log to a file")

parser.add_argument('--max-targets',
                    action='store',
                    dest='max_n_targets',
                    metavar='<max targets>',
                    type=int,
                    default=3,
                    help="Maximum number of targets to be found in a "
                         "single image. Default 3")

parser.add_argument('--save-plots',
                    action='store_true',
                    dest='save_plots',
                    help="Save all plots in a directory")

parser.add_argument('--plot-results',
                    action='store_true',
                    dest='plot_results',
                    help="Show wavelength calibrated spectrum at the end.")

args = parser.parse_args(args=arguments)

print(args)
