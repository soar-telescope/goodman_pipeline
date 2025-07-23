import argparse

from importlib.metadata import version

__version__ = version("goodman_pipeline")


def get_photometry_config_args(arguments=None):

    parser = argparse.ArgumentParser(
        description="Obtain photometry measurements")

    parser.add_argument('filename', nargs='?', help="Path to file to process.")

    parser.add_argument(
        '--aperture-radius',
        default=4.0,
        type=float,
        help="Aperture radius in pixels. (default=4.0)")

    # parser.add_argument(
    #     '--aperture-type',
    #     default='fixed',
    #     choices=['fixed', 'variable'],
    #     type=str,
    #     help="Aperture type to perform Photometry. (default 'fixed')")

    parser.add_argument(
        '--initial-fwhm',
        default=3.0,
        type=float,
        help="Initial FWHM in pixels for source detection. (default=3.0)")

    parser.add_argument(
        '--detection-threshold',
        default=5.0,
        type=float,
        help="Detection threshold above noise factor (default=5.0)")

    parser.add_argument(
        '--gaia-sources-limit',
        default=5000,
        type=int,
        help="Maximum number of Gaia sources to process. (default=5000)")

    parser.add_argument(
        '--imaging-filter-keyword',
        default='FILTER',
        type=str,
        help="FITS header keyword for imaging filter (default: 'FILTER').")
    #
    # parser.add_argument(
    #     '--ra-keyword',
    #     default='OBSRA',
    #     type=str,
    #     help="FITS header keyword for observed RA (default: 'OBSRA').")
    #
    parser.add_argument(
        '--gaia-photometry-column',
        default='',
        type=str,
        help="GAIA's filter corresponding column name.")

    parser.add_argument(
        '--overwrite',
        action='store_true',
        help="Allow overwriting of existing photometry results.")

    parser.add_argument(
        '--aperture-curve-of-growth',
        action='store_true',
        help="Run aperture curve of growth diagnostics and exit.")

    parser.add_argument(
        '--debug',
        action='store_true',
        help="Enable debug mode.")

    parser.add_argument(
        '-v', '--version',
        action='store_true')

    args = parser.parse_args(args=arguments)

    if args.version:
        parser.exit(status=0, message=__version__)


    if not args.filename:
        parser.print_help()
        parser.exit(status=0, message="\nPlease specify a filename to process.\n")

    return args
