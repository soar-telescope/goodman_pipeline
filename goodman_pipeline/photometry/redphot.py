from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from goodman_pipeline.images.data_classifier import DataClassifier
from goodman_pipeline.images.night_organizer import NightOrganizer
from goodman_pipeline.images.image_processor import ImageProcessor

#wrong:
from astropy.io           import fits
from .aperture_phot import AperturePhotometry
#from goodman_pipeline.photometry import

import os
import sys
import shutil
import argparse
import glob
import logging
import matplotlib

log = logging.getLogger(__name__)

try:
    matplotlib.use('Qt5Agg')
except ImportError as error:
    log.warning(error)

__version__ = __import__('goodman_pipeline').__version__


def get_args(arguments=None):
    """Get command line arguments.

    The list of arguments can be obtained by using the argument ``--help``.
    All the arguments start with two dashes and single-character arguments where
    avoided in order to eliminate confusion.

    Args:
        arguments (list): A list containing the arguments as elements.

    Returns:
        args (object): argparse instance. Contains all the arguments as
            attributes

    """
    parser = argparse.ArgumentParser(
        description="Goodman CCD Photometry - CCD photometry for Goodman "
                    "Imaging data."
                    "\nPipeline Version: {:s}".format(__version__))

    parser.add_argument('--auto-clean',
                        action='store_true',
                        dest='auto_clean',
                        help="Automatically clean reduced data directory")

    parser.add_argument('--debug',
                        action='store_true',
                        dest='debug_mode',
                        help="Show detailed information of the process.")

    parser.add_argument('--fwhm',
                        action='store',
                        metavar='<value>',
                        default='15',
                        dest='fwhm',
                        type=float,
                        help="FHWM for star finding. Default 15.")

    parser.add_argument('--sigma-threshold',
                        action='store',
                        default=5.,
                        dest='sigma_threshold',
                        metavar='<value>',
                        type=float,
                        help="Detection limit (in sigmas). Default 5")

    parser.add_argument('--source-id',
                        action='store',
                        default=0.,
                        dest='source_id',
                        metavar='<value>',
                        type=int,
                        help="Source ID for FWHM estimation. Default 5")

    parser.add_argument('--pixscale',
                        action='store',
                        default=10.,
                        dest='pixscale',
                        metavar='<value>',
                        type=float,
                        help="Pixel scale for FWHM estimation. Default 10")

    parser.add_argument('--stamp-radius',
                        action='store',
                        default=50.,
                        dest='stamp_radius',
                        metavar='<value>',
                        type=int,
                        help="Stamp radius for FWHM estimation. Default 50")


    parser.add_argument('--psf-model',
                        action='store',
                        default='moffat',
                        type=str,
                        metavar='<psf_model>',
                        dest='psf_model',
                        choices=['moffat', 'gaussian'],
                        help="Choose a model to fit the psf. Choices are: 'moffat',"
                             "'gaussian'. Default 'moffat'")

    parser.add_argument('--r-in',
                        action='store',
                        type=float,
                        metavar='<value>',
                        default=10.0,
                        help='Inner radius for aperture photometry. Default 10',
                        dest='r_in')

    parser.add_argument('--r-out',
                        action='store',
                        type=float,
                        metavar='<value>',
                        default=None,
                        help='Outer radius for aperture photometry (optional)',
                        dest='r_in')
    #This can be obtained from image
    parser.add_argument('--gain',
                        action='store',
                        type=float,
                        metavar='<value>',
                        default=1.48,
                        help="Gain for photometry. Default 1.48 "
                             "(eliminate this function ASAP)",
                        dest='gain')

    parser.add_argument('--plot', action='store_true', help='Plot results')

    parser.add_argument('--raw-path',
                        action='store',
                        metavar='<raw_path>',
                        default='./',
                        type=str,
                        help="Path to raw data.")

    parser.add_argument('--red-path',
                        action='store',
                        metavar='<red_path>',
                        type=str,
                        default='./RED',
                        help="Path to reduced data.")

    parser.add_argument('--version',
                        action='store_true',
                        dest='show_version',
                        help="Show current version of the Goodman Pipeline")

    args = parser.parse_args(args=arguments)

    return args


class MainApp(object):

    def __init__(self):
        """This method initializes the MainApp class

        The main task of this method is to call the get_args function that
        returns an argparse object.
        The arguments will be obtained here and they will be available for all
        execution of the program.
        Other attributes will be initialized as None.
        """
        self.log = logging.getLogger(__name__)
        self.args = get_args()
        self.log.debug('Initializing AperturePhotometry instance')

    def run(self):
        """Main execution method for the application.

        This method checks if the raw_path points to a directory. If so, it loops
        through all FITS files (with extension '.fits') within the directory and
        processes each file individually using the existing functionalities.
        If raw_path points to a single file, it behaves as before.

        Raises:
            Exception: If an error occurs during processing.
        """
        # Check if raw_path is a directory
        if os.path.isdir(self.args.raw_path):
            # Loop through all FITS files in the directory
            for fits_file in glob.glob(os.path.join(self.args.raw_path, "*.fits")):
                print(f"Processing file: {fits_file}")
                try:
                    ap = AperturePhotometry()
                    ap.process_single_file(fits_file, self.args)
                except Exception as e:
                    self.log.error(f"Error processing file {fits_file}: {e}")
        else:
            # Handle single FITS file (existing behavior)
            try:
                ap.process_single_file(self.args.raw_path, self.args)
                # ... (rest of your processing logic here)
                # ...
            except Exception as e:
                print(f"Error processing file: {e}")

if __name__ == '__main__':  # pragma: no cover
    main_app = MainApp()
    main_app()
