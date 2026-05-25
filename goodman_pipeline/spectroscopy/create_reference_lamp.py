import argparse
import os
import sys

import astropy.units as u
import logging

from importlib.metadata import version

from astropy.nddata import CCDData
from matplotlib import pyplot as plt

from goodman_pipeline.core import (get_lines_in_lamp, get_spectral_characteristics)
from goodman_pipeline.wcs import WCS


__version__ = version('goodman_pipeline')

def get_args(arguments=None):
    log = logging.getLogger()

    parser = argparse.ArgumentParser(
        description="Creates reference lamp.\nPipeline Version: {:s}".format(__version__))
    parser.add_argument("--comparison-lamp", action="store", help="Comparison lamp file name")
    parser.add_argument("--reference-lamp", action="store", help="Already calibrated comparison lamp file name.")

    parser.add_argument("-v", "--version", action="version", version=__version__)

    args = parser.parse_args(args=arguments)
    if not args.comparison_lamp:
        log.error("Comparison lamp file name not specified.")
        parser.print_help()
        sys.exit("Please specify a comparison lamp file name.")

    return args

class CreateReferenceLamp:

    def __init__(self):
        self.log = logging.getLogger(__name__)
        self.args = None
        self.pixel_size = 15 * u.micrometer
        self.instrument_focal_length = 377.3 * u.mm
        self.wcs = WCS()
        self.pixel = []
        self.angstrom = []
        self.comparison_lamp = None
        self.reference_lamp = None
        self.fig, (self.ax_ref, self.ax_comp) = plt.subplots(nrows=2, ncols=1, figsize=(16, 10))

    def __call__(self, args=None):
        if args is None:
            self.args = get_args()
        else:
            self.args = args

        if not os.path.exists(self.args.comparison_lamp):
            self.log.error("Comparison lamp file not found.")
            sys.exit("Please specify a comparison lamp file name.")

        self.comparison_lamp = CCDData.read(self.args.comparison_lamp, unit=u.adu)
