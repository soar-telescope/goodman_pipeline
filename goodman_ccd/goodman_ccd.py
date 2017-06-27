from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import glob
import logging
import os
import re
import time

import argparse
import numpy as np
import matplotlib
matplotlib.use('GTK3Agg')

from ccdproc import ImageFileCollection
from matplotlib import pyplot as plt

from .core import print_default_args
from .data_classifier import DataClassifier
from .night_organizer import NightOrganizer
from .image_processor import ImageProcessor

__author__ = 'David Sanmartim'
__date__ = '2016-07-15'
__version__ = "0.1"
__email__ = "dsanmartim@ctio.noao.edu"
__maintainer__ = "Simon Torres"

FORMAT = '%(levelname)s: %(asctime)s: %(module)s.%(funcName)s: %(message)s'
# DATE_FORMAT = '%m/%d/%Y %I:%M:%S%p'
DATE_FORMAT = '%I:%M:%S%p'
LOG_FILENAME = 'goodman_ccd.log'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt=DATE_FORMAT)
log = logging.getLogger('goodmanccd')


def get_args(arguments=None):
    """Get command line arguments.

    Returns:
        args (object): Arparse object. Contains all the arguments as attributes

    """
    # global log
    # Parsing Arguments ---
    parser = argparse.ArgumentParser(description="Goodman CCD Reduction - CCD"
                                                 "reductions for "
                                                 "Goodman spectroscopic data")

    parser.add_argument('-c', '--cosmic',
                        action='store_true',
                        dest='clean_cosmic',
                        help="Clean cosmic rays from science data.")

    # TODO (simon): Add argument to use calibration data from other day

    parser.add_argument('--ignore-bias',
                        action='store_true',
                        dest='ignore_bias',
                        help="Ignore bias correction")

    parser.add_argument('--auto-clean',
                        action='store_true',
                        dest='auto_clean',
                        help="Automatically clean reduced data directory")

    parser.add_argument('--saturation',
                        action='store',
                        default=55000.,
                        dest='saturation_limit',
                        metavar='<Value>',
                        help="Saturation limit. Default to 55.000 ADU (counts)")

    parser.add_argument('--raw-path',
                        action='store',
                        metavar='raw_path',
                        default='./',
                        type=str,
                        help="Path to raw data.")

    parser.add_argument('--red-path',
                        action='store',
                        metavar='red_path',
                        type=str,
                        default='./RED',
                        help="Path to reduced data.")

    parser.add_argument('--debug',
                        action='store_true',
                        dest='debug_mode',
                        help="Show detailed information of the process.")

    parser.add_argument('--log-to-file',
                        action='store_true',
                        dest='log_to_file',
                        help="Write log to a file.")

    parser.add_argument('--flat-normalize',
                        action='store',
                        default='simple',
                        type=str,
                        metavar='<Normalization Method>',
                        dest='flat_normalize',
                        choices=['mean', 'simple', 'full'],
                        help='Choose a method to normalize the master flat for'
                             'spectroscoy. Choices are: mean, simple (model) '
                             'and full (fits model to each line).')

    parser.add_argument('--flat-norm-order',
                        action='store',
                        default=15,
                        type=int,
                        metavar='<Order>',
                        dest='norm_order',
                        help='Defines the order of the model to be fitted.')

    parser.add_argument('--dcr-par-dir',
                        action='store',
                        default='files/',
                        metavar='<dcr.par directory>',
                        dest='dcr_par_dir',
                        help="Directory of default dcr.par file.")

    parser.add_argument('--keep-cosmic-files',
                        action='store_false',
                        dest='keep_cosmic_files',
                        help="After cleaning cosmic rays with dcr, do not "
                             "remove the input file and the cosmic rays file.")


    args = parser.parse_args(args=arguments)

    if args.log_to_file:
        log.info('Logging to file {:s}'.format(LOG_FILENAME))
        file_handler = logging.FileHandler(filename=LOG_FILENAME)
        file_handler.setLevel(level=logging.INFO)
        formatter = logging.Formatter(fmt=FORMAT, datefmt=DATE_FORMAT)
        file_handler.setFormatter(fmt=formatter)
        log.addHandler(file_handler)
    if args.debug_mode:
        log.info('Changing log level to DEBUG.')
        log.setLevel(level=logging.DEBUG)
    if os.path.isdir(args.raw_path):
        args.raw_path = os.path.abspath(args.raw_path)
        log.debug(os.path.abspath(args.raw_path))
    else:
        parser.print_help()
        parser.exit("Raw data folder doesn't exist")

    # updated full path for default dcr.par file. If it doesn't exist it will
    # create an empty one.
    dcr_par_full_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 args.dcr_par_dir)
    if not os.path.isdir(dcr_par_full_path) or args.dcr_par_dir != 'files/':
        log.info("dcr.par default location doesn't exist.")
        try:
            os.path.os.makedirs(dcr_par_full_path)
            log.info('Created dcr.par empty directory: %s', dcr_par_full_path)
            args.dcr_par_dir = dcr_par_full_path
        except OSError as err:
            log.error(err)
    else:
        args.dcr_par_dir = dcr_par_full_path
    # print_default_args(args)
    return args


class MainApp(object):

    def __init__(self, args=None):
        """This method initalizes the MainApp class

        The main task of this method is to call the get_args function that
        returns an argparse object.
        The arguments will be obtained here and they will be available for all
        execution of the program.
        Other attributes will be initilized as None.
        """
        if args is None:
            self.args = get_args()
        else:
            self.args = args
        self.data_container = None
        self.full_path = None
        # self.instrument = None
        # self.technique = None
    
    def __call__(self):
        """Call method for MainApp

        From the arguments this method finds the raw_path attribute and checks
        its contents for the existance of files containing the '.fits' string.
        If there is none it will assume every item is a different data directory
        and they will be treated independently. If there are '.fits' files the
        program will assume is a single data directory.
        Any subdirectory will be ignored.

        """

        folders = glob.glob(os.path.join(self.args.raw_path, '*'))
        if any('.fits' in item for item in folders):
            folders = [self.args.raw_path]
        for data_folder in folders:
            if not os.path.isdir(data_folder):
                continue

            self.args.raw_path = data_folder

            try:
                log.debug('Initializing DataClassifier Class')
                night_sorter = DataClassifier(self.args)
                log.debug('Calling night_sorter Instance of DataClassifier')
                night_sorter()
                # self.instrument = night_sorter.instrument
                # self.technique = night_sorter.technique
                # print(self.instrument)
                # print(self.technique)
            except AttributeError as error:
                log.error(error)
                log.error('Empty or Invalid data directory:'
                          '{:s}'.format(data_folder))
                continue

            # # check start
            # self.args.raw_path = data_folder
            if self.args.red_path == './RED' or len(folders) > 1:

                log.info('No special reduced data path defined. '
                         'Proceeding with defaults.')

                if self.args.raw_path not in self.args.red_path:
                    self.args.red_path = os.path.join(self.args.raw_path, 'RED')
                    # print(self.args.red_path)

            if os.path.isdir(self.args.red_path):
                if os.listdir(self.args.red_path) != []:
                    log.warning('Reduced Data Path is not empty')
                    if self.args.auto_clean:
                        for _file in os.listdir(self.args.red_path):
                            os.unlink(os.path.join(self.args.red_path, _file))

                        log.info('Cleaned Reduced data directory:'
                                 ' {:s}'.format(self.args.red_path))
                    else:
                        log.error('Please clean the reduced data folder or '
                                  'use --auto-clean')
                        break
                self.args.red_path = os.path.abspath(self.args.red_path)
                log.debug(os.path.abspath(self.args.red_path))
            else:
                try:
                    log.warning("Reduction folder doesn't exist.")
                    os.mkdir(os.path.abspath(self.args.red_path))
                    log.info('Created reduced data directory!')
                    log.info(os.path.abspath(self.args.red_path))
                except OSError as error:
                    log.error(error)
            # check ends

            # print(night_sorter.nights_dict)
            for night in night_sorter.nights_dict:
                # print(night_sorter.nights_dict[night])
                # night_organizer = \
                #     NightOrganizer(args=self.args,
                #                    night_dict=night_sorter.nights_dict[night])
                # nd = Night Dictionary
                nd = night_sorter.nights_dict[night]
                log.debug('Initializing NightOrganizer Class')
                night_organizer = NightOrganizer(full_path=nd['full_path'],
                                                 instrument=nd['instrument'],
                                                 technique=nd['technique'],
                                                 ignore_bias=
                                                 self.args.ignore_bias)

                log.debug('Calling night_organizer instance')
                self.data_container = night_organizer()
                if self.data_container is None or self.data_container is None:
                    log.error('Discarding night ' + str(night))
                    break
                process_images = ImageProcessor(self.args, self.data_container)
                process_images()


if __name__ == '__main__':
    main_app = MainApp()
    main_app()
