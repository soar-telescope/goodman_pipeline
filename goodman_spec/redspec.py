#!/usr/bin/env python2
# -*- coding: utf8 -*-
"""Pipeline for GOODMAN spectra Extraction.

This program finds reduced images, i.e. trimmed, bias subtracted, flat fielded,
etc. that match the <pattern> in the source folder, then classify them in two
groups: Science or Lamps. For science images, finds the spectrum or spectra and
traces it doing some fit.
Simon Torres 2016-06-28

"""
# TODO (simon): Change all astropy.io.fits to astropy.CCDData.read
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging
import os
import sys
import textwrap
import time
import warnings

import argparse
import ccdproc as ccd
import matplotlib
import numpy as np
import pandas as pd
matplotlib.use('GTK3Agg')

from ccdproc import CCDData
from matplotlib import pyplot as plt

from goodman_ccd.core import classify_spectroscopic_data
from .wavelength import process_spectroscopy_data


warnings.filterwarnings('ignore')
FORMAT = '%(levelname)s: %(asctime)s:%(module)s.%(funcName)s: %(message)s'
DATE_FORMAT = '%I:%M:%S%p'
LOG_FILENAME = 'goodman_spec.log'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt=DATE_FORMAT)
log = logging.getLogger('redspec')

__author__ = 'Simon Torres'
__date__ = '2016-06-28'
__version__ = "0.1"
__email__ = "storres@ctio.noao.edu"
__status__ = "Development"


def get_args(arguments=None):
    """Handles the argparse library and returns the arguments

    The possible arguments to be used are:

    --data-path: has to be the source directory, where the data (images) is.
        The location is `self.source` and the default value is ./
    --proc-path: is the destination where all new files/data will be
        placed. The location is `self.destiny` and the default value is ./
    --search-pattern: the pattern that matches the reduced data that
        will be processed. The location is `self.pattern` default value is
        fc\_
    --interactive: Interactive Wavelength Solution. Enabled by default
    --output-prefix: Prefix to use to name wavelength calibrated spectrum
    --extraction: Spectrum extraction type. Simple or optimal.
    --reference-files: Directory of reference files location
    --debug: Show plots for intermediate steps. Also prints more messages.
    --log-to-file: Write the log to a file.
    --save-plots: Save all plots ina directory
    --plot-results: Show a plot of the wavelength calibrated data at the
        end.


    Returns:
        An object that contains all the variables parsed through the argument
        system

    """
    leave = False
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
                        dest='destiny',
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
                        default='refdata/',
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

    parser.add_argument('--save-plots',
                        action='store_true',
                        dest='save_plots',
                        help="Save all plots in a directory")

    parser.add_argument('--plot-results',
                        action='store_true',
                        dest='plot_results',
                        help="Show wavelength calibrated spectrum at the end.")

    args = parser.parse_args(args=arguments)

    if args.log_to_file:
        log.info('Logging to file {:s}'.format(LOG_FILENAME))
        file_handler = logging.FileHandler(LOG_FILENAME)
        formatter = logging.Formatter(fmt=FORMAT, datefmt=DATE_FORMAT)
        file_handler.setFormatter(fmt=formatter)
        log.addHandler(file_handler)

    if args.debug_mode:
        log.info('Changing log level to DEBUG.')
        log.setLevel(level=logging.DEBUG)

    # get full path for reference files directory
    ref_full_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 args.reference_dir)
    if not os.path.isdir(ref_full_path):
        log.info("Reference files directory doesn't exist.")
        try:
            os.path.os.makedirs(ref_full_path)
            log.info('Reference Files Directory is: %s', ref_full_path)
            args.reference_dir = ref_full_path
        except OSError as err:
            log.error(err)
    else:
        args.reference_dir = ref_full_path

    if not os.path.isdir(args.source):
        leave = True
        log.error("Source Directory doesn't exist.")

    if not os.path.isdir(args.destiny):
        leave = True
        log.error("Destination folder doesn't exist.")
        try:
            os.path.os.makedirs(args.destiny)
            log.info('Destination folder created: %s', args.destiny)
        except OSError as err:
            log.error(err)

    if leave:
        parser.print_help()
        parser.exit("Leaving the Program.")
    # print_default_args(args)
    return args


class MainApp(object):
    """Defines and initialize all important variables for processing the data

    The MainApp class controls the way the night is organize to its further
    processing. It also sets the appropriate parameters that will allow for a
    smooth working in all the other modules.

    """
    def __init__(self, args=None):

        if args is None:
            self.args = get_args()
        else:
            self.args = args

        self.wavelength_solution_obj = None

    def __call__(self):

        # data_container instance of NightDataContainer defined in core
        data_container = classify_spectroscopic_data(
            path=self.args.source,
            search_pattern=self.args.pattern)

        # print('data_container.bias')
        # print(data_container.bias)
        # print('data_container.day_flats')
        # print(data_container.day_flats)
        # print('data_container.dome_flats')
        # print(data_container.dome_flats)
        # print('data_container.sky_flats')
        # print(data_container.sky_flats)
        # print('data_container.data_groups')
        # print(data_container.data_groups)
        # print('data_container.spec_groups')
        # print(data_container.spec_groups)

        self.wavelength_solution_obj = process_spectroscopy_data(
            data_container=data_container,
            args=self.args,
            extraction_type=self.args.extraction_type)


if __name__ == '__main__':
    MAIN_APP = MainApp()
    try:
        MAIN_APP()
    except KeyboardInterrupt:
        sys.exit(0)
