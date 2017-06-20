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
import sys
import os
import numpy as np
import time
import textwrap
import ccdproc as ccd
import pandas as pd
import argparse
import logging
import matplotlib
matplotlib.use('GTK3Agg')
import matplotlib.pyplot as plt
# from astropy import log
import warnings
from ccdproc import CCDData
from .process import Process, SciencePack
from .wavelength import WavelengthCalibration, process_spectroscopy_data
from goodman_ccd.core import (print_spacers,
                              ra_dec_to_deg,
                              convert_time,
                              print_default_args,
                              classify_spectroscopic_data)


warnings.filterwarnings('ignore')
FORMAT = '%(levelname)s: %(asctime)s:%(module)s: %(message)s'
DATE_FORMAT = '%m/%d/%Y %I:%M:%S%p'
LOG_FILENAME = 'goodman_spec.log'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt=DATE_FORMAT)
log = logging.getLogger('redspec')

__author__ = 'Simon Torres'
__date__ = '2016-06-28'
__version__ = "0.1"
__email__ = "storres@ctio.noao.edu"
__status__ = "Development"



def get_args():
    """Handles the argparse library and returns the arguments

    Returns:
        An object that contains all the variables parsed through the argument
        system
        The possible arguments to be returned are:

        -p or --data-path: has to be the source directory, where the data (images) is.
                the location is **self.source**
                default value is ./
        -d or --proc-path: is the destination where all new files/data will be placed.
                the location is self.destiny
                default value is ./
        -s or --search-pattern: the pattern that matches the reduced data that will be processed.
                the location is self.pattern
                default value is fc\_
        -m or --proc-mode: is one of the predefined observing modes and the options are:
                0: One or more lamps taken during the beginning or end of the night, i.e. single
                calibration to all data in that night
                1: One or more lamps around every science exposure.
                default value is 0
        -r or --reference-lamp: Name of reference lamp file for mode 0. If not present, the first one in the list
                will be selected
        -l or --lamp-file: Name of the ASCII file that contains the relation between science files and lamp
                files. An example is depicted below. Note that the lamps can be repeated.
                    #example of how the file should look
                    science_target_01.fits lamp_001.fits
                    science_target_02.fits lamp_001.fits
                    science_target_03.fits lamp_002.fits

                the location is self.lamp_file
                default value is lamps.txt

        -i or --non-interactive: Interactive Wavelength Solution. Enabled by default
        -o or --output-prefix: Prefix to use to name wavelength calibrated spectrum
        -R or --reference-files: Directory of reference files location
        --plots-enabled: Show plots for intermediate steps. For debugging only.

    """
    leave = False
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent(
                                         '''Extracts goodman spectra and does wavelength calibration.\n\n\
Supported Processing Modes are:
<0>: (Default) reads lamps taken at the begining or end of the night.\n\
<1>: one or more lamps around science exposure.
'''))
# \n\
# <2>: ASCII file describing which science target uses which lamp.\n\
# <3>: No lamps. Uses the sky lines

    parser.add_argument('-p', '--data-path',
                        action='store',
                        default='./',
                        type=str,
                        metavar='<Source Path>',
                        dest='source',
                        help='Path for location of raw data. Default <./>')

    parser.add_argument('-d', '--proc-path',
                        action='store',
                        default='./',
                        type=str,
                        metavar='<Destination Path>',
                        dest='destiny',
                        help='Path for destination of processed data. Default <./>')

    parser.add_argument('-s', '--search-pattern',
                        action='store',
                        default='cfzsto_',
                        type=str,
                        metavar='<Search Pattern>',
                        dest='pattern',
                        help="Pattern for matching the goodman's reduced data.")

    # parser.add_argument('-m', '--proc-mode',
    #                     action='store',
    #                     default=0,
    #                     type=int,
    #                     metavar='<Processing Mode>',
    #                     dest='procmode',
    #                     choices=[0, 1],
    #                     help='Defines the mode of matching lamps to science targets.')

    parser.add_argument('-r', '--reference-lamp',
                        action='store',
                        default='',
                        type=str,
                        metavar='<Reference Lamp>',
                        dest='lamp_all_night',
                        help="Name of reference lamp file for mode 0.\
                         If not present, the first one in the list will be selected")

    # parser.add_argument('-l', '--lamp-file',
    #                     action='store',
    #                     default='lamps.txt',
    #                     type=str,
    #                     metavar='<Lamp File>',
    #                     dest='lamp_file',
    #                     help="Name of an ASCII file describing which science"
    #                          "target uses which lamp. default <lamp.txt>")

    parser.add_argument('-o', '--output-prefix',
                        action='store',
                        default='g',
                        metavar='<Out Prefix>',
                        dest='output_prefix',
                        help="Prefix to add to calibrated spectrum.")

    parser.add_argument('-R', '--reference-files',
                        action='store',
                        default='refdata/',
                        metavar='<Reference Dir>',
                        dest='reference_dir',
                        help="Directory of Reference files location")

    parser.add_argument('-i', '--interactive',
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

    args = parser.parse_args()
    if args.debug_mode:
        log.info('Changing log level to DEBUG.')
        log.setLevel(level=logging.DEBUG)
    if args.log_to_file:
        log.info('Logging to file {:s}'.format(LOG_FILENAME))
        file_handler = logging.FileHandler(LOG_FILENAME)
        formatter = logging.Formatter(fmt=FORMAT, datefmt=DATE_FORMAT)
        file_handler.setFormatter(fmt=formatter)
        log.addHandler(file_handler)

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

    # if args.procmode == 2:
    #     # print(args.source + args.lamp_file)
    #     if not os.path.isfile(os.path.join(args.source, args.lamp_file)):
    #         if args.lamp_file == 'lamps.txt':
    #             leave = True
    #             log.error("Default <lamp file> doesn't exist.")
    #             log.error("Please define a <lamp file> using the flags -l or --lamp-file")
    #             log.error("or make sure you entered the right observing mode.")
    if leave:
        parser.print_help()
        parser.exit("Leaving the Program.")
    print_default_args(args)
    return args


class MainApp(object):
    """Defines and intialize all important variables for processing the data

    The MainApp class controls the way the night is organize to its further
    processing. It also sets the appropriate parameters that will allow for a
    smooth working in all the other modules.

    """
    def __init__(self, args=None):

        if args is None:
            self.args = get_args()
        else:
            self.args = args
        # self.image_collection = pd.DataFrame
        # self.night = self.set_night()
        # self.extracted_data = None
        # self.wsolution = None
        # self.calibration_lamp = None
        self.wavelength_solution_obj = None

    def __call__(self):

        # data_container instance of NightDataContainer defined in core
        data_container = classify_spectroscopic_data(
            path=self.args.source,
            search_pattern=self.args.pattern)

        # TODO (simon): add extraction to arguments.
        extracted, comps = process_spectroscopy_data(
            data_container=data_container,
            args=self.args,
            extraction_type='simple')


# class Night(object):
#     """Stores all data relevant to the night being processed
#
#     Note:
#         The night class stores the data relative to single observing night
#         therefore this software works on a per-night basis
#     """
#
#     def __init__(self, date, args):
#         """Initialize Night class
#
#         The night class will store filename of science images as well as lamp
#         images. It also stores the arguments and also the wavelength solution
#         for the night.
#
#         Args:
#             date (str): Date of the night being processed
#             args (object): argparse instance containing all the arguments the
#             program was started with
#         """
#         # source, destiny, pattern, mode, lamps
#         self.all = []
#         self.date = date
#         self.sci = []
#         self.lamp = []
#         self.args = args
#         self.sci_targets = []
#         self.telescope = False
#         self.night_wsolution = None
#         self.night_calibration_lamp = None
#         self.gratings = None
#
#     def add_sci(self, in_sci):
#         """Adds science object to list"""
#         for new_sci in in_sci:
#             if self.args.pattern == new_sci[0:len(self.args.pattern)]:
#                 self.sci.append(new_sci)
#                 self.all.append(new_sci)
#             else:
#                 log.info('Image %s rejected because does not match pattern!', new_sci)
#         # self.sci_targets = [[]] * len(self.sci)
#
#     def add_lamp(self, in_lamp):
#         """Adds lamp objects to list"""
#         for new_lamp in in_lamp:
#             if self.args.pattern == new_lamp[0:len(self.args.pattern)]:
#                 self.lamp.append(new_lamp)
#                 self.all.append(new_lamp)
#             else:
#                 log.info('Image %s rejected because does not match pattern!', new_lamp)
#
#     def add_sci_object(self, sci_obj):
#         """Appends a ScienceObject to the class attribute sci_targets"""
#         self.sci_targets.append(sci_obj)
#         log.info("Added science object %s", sci_obj.name)
#
#     def set_gratings(self, gratings):
#         """Adds an array of the names of all the gratings observed in the night"""
#         self.gratings = gratings
#
#     def set_night_calibration_lamp(self, calibration_lamp):
#         """Sets the filename of the calibration lamp as an attribute"""
#         self.night_calibration_lamp = calibration_lamp
#
#     def is_telescope(self):
#         """Sets the operation mode as Telescope Mode"""
#         self.telescope = True
#
#     def set_night_wsolution(self, wsolution):
#         """Sets the wavelength solution as a class attribute"""
#         if wsolution is not None:
#             self.night_wsolution = wsolution
#             log.info('Night Solution Defined')
#         else:
#             log.error("Wavelength solution still can't be defined")
#
#
# class ScienceObject(object):
#     """class that defines a science object attributes
#
#     Sci objects, for science, are the main targets which have lamps for
#     calibration. Their atritutes are: the name, the file name, the
#     observation time, right ascension and declination. The same information
#     can be obtained for the lamps but they are treated as lists
#     The list are created in a way that allows the use of elements' index
#     as the correlator between atributes of the lamp.
#
#     Attributes:
#         name (str): science object name
#         file_name (str): file name
#         obs_time (str): observing time in the format yyyy-mm-ddThh:mm:ss.ss for instance 2016-03-20T23:54:15.96
#         right_ascension (float): right ascension in degrees
#         declination (float): declination in degrees
#         lamp_count (int): lamps count
#         lamp_file (list): every element is a string with the file name of the lamp
#         lamp_type (list): every element is a string with the OBJECT value of the lamp i.e Cu, HgAr, etc
#         lamp_ra (list): every element is a float with lamp's right ascension in degrees
#         lamp_dec (list): every element is a float with lamp's declination in degrees
#
#     """
#
#     def __init__(self, name, file_name, obs_time, right_ascension, declination, grating):
#         self.name = name
#         self.file_name = file_name
#         self.obs_time = obs_time
#         self.right_ascension = right_ascension
#         self.declination = declination
#         self.lamp_count = 0
#         self.lamp_file = []
#         self.lamp_type = []
#         self.lamp_ra = []
#         self.lamp_dec = []
#         self.grating = grating
#         self.no_targets = 0
#
#     def add_lamp(self, new_lamp, new_type, right_ascension, declination):
#         """Adds a lamp to the science object
#
#         Args:
#             new_lamp (str): new lamp file name
#             new_type (str): new lamp type
#             right_ascension (float): right ascension in degrees
#             declination (float): declination in degrees
#
#         """
#         self.lamp_file.append(new_lamp)
#         self.lamp_type.append(new_type)
#         self.lamp_ra.append(right_ascension)
#         self.lamp_dec.append(declination)
#         self.lamp_count = int(len(self.lamp_file))
#
#     def update_no_targets(self, new_value=None, add_one=False):
#         """Update number of spectra in an image
#
#         An spectral image may contain multiple science target's spectra this
#         method is set to update its number as a class attribute. There are two
#         ways it can work. Add one to the existing number or set a new one.
#
#         Args:
#             new_value (int): New value for number of targets
#             add_one (bool): If True increase the count by one.
#         """
#         if add_one:
#             self.no_targets += 1
#             log.debug('Number of targets is %s', self.no_targets)
#         elif new_value is not None:
#             self.no_targets = new_value
#         else:
#             log.debug('Nothing to do: new_value: %s add_one: %s', new_value, str(add_one))
#
#     def print_all(self):
#         """Prints all the relevant attributes of the object
#
#         Note:
#             this method is mainly used for development purposes
#
#         """
#         print(' ')
#         log.info("Name: %s", self.name)
#         log.info("File: %s", self.file_name)
#         log.info("Obs-T: %s", self.obs_time)
#         if self.lamp_count > 0:
#             log.info("Lamp N: %s", self.lamp_count)
#             for i in range(self.lamp_count):
#                 log.info("Lamp %s: %s", (i + 1), self.lamp_file[i])
#                 log.info("Type %s: %s", (i + 1), self.lamp_type[i])
#         # time.sleep(10)


if __name__ == '__main__':
    MAIN_APP = MainApp()
    try:
        MAIN_APP()
    except KeyboardInterrupt:
        sys.exit(0)
