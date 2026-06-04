import sys
import argparse
import logging

from importlib.metadata import version

log = logging.getLogger(__name__)

__version__ = version('goodman_pipeline')


def get_create_reference_lamp_args(arguments=None):

    parser = argparse.ArgumentParser(
        description="Creates reference lamp.\nPipeline Version: {:s}".format(__version__))
    parser.add_argument("--comparison-lamp", action="store", help="Comparison lamp file name. This is your NEW lamp.")
    parser.add_argument("--reference-lamp", action="store", default=None, help="Already calibrated comparison lamp file name.")
    parser.add_argument("--plots-theme", action="store", default="dark", choices=["light", "dark"], help="Choose a theme for plotting, default is dark.")
    parser.add_argument("--screen-size", action="store", default='large', choices=['small', 'medium', 'large'], help="Choose a screen size for sizing the plots, default is large.")
    parser.add_argument("--comp-intensity-start", action="store", default=None, help="Override y-axis start value for comparison lamp.")
    parser.add_argument("--comp-intensity-end", action="store", default=None, help="Override y-axis end value for comparison lamp.")
    parser.add_argument("--ref-wavelength-start", action="store", default=None, help="Override wavelength start value for reference lamp.")
    parser.add_argument("--ref-wavelength-end", action="store", default=None, help="Override wavelength end value for reference lamp.")
    parser.add_argument("--debug", action="store_true", default=False, help="Enable debug mode.")
    parser.add_argument("-v", "--version", action="version", version=__version__)
    args = parser.parse_args(args=arguments)

    if not args.comparison_lamp:
        log.error("Comparison lamp file name not specified.")
        parser.print_help()
        sys.exit("Please specify a comparison lamp file name.")

    return args
