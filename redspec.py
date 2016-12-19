#!/usr/bin/python2.7
"""Pipeline for GOODMAN spectra Extraction.

This program finds reduced images, i.e. trimmed, bias subtracted, flat fielded, etc. that match the <pattern>
in the source folder, then classify them in two groups: Science or Lamps. For science images, finds the spectrum
or spectra and traces it doing some fit.
Simon Torres 2016-06-28

"""
from __future__ import print_function
import sys
import os
import numpy as np
import time
import textwrap
import ccdproc as ccd
import pandas as pd
import argparse
import logging
import warnings
from process import Process, SciencePack
from wavelength import WavelengthCalibration

warnings.filterwarnings('ignore')
FORMAT = '%(asctime)s:%(levelname)s:%(module)s: %(message)s'
DATE_FORMAT = '%m/%d/%Y %I:%M:%S%p'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt=DATE_FORMAT)
log = logging.getLogger('redspec')

__author__ = 'Simon Torres'
__date__ = '2016-06-28'
__version__ = "0.1"
__email__ = "storres@ctio.noao.edu"
__status__ = "Development"


class MainApp(object):
    """Defines and intialize all important variables for processing the data

    The MainApp class controls the way the night is organize to its further processing. It also sets the appropriate
    parameters that will allow for a smooth working in all the other modules.

    """
    def __init__(self):
        """Initalization of important parameters

        Initializes the list of images using ccdproc.ImageFileCollection and pandas the get the arguments that define
        the working of the pipeline using arpargse and instantiate a Night class, an object that will store relevant
        information of the observed night being processed.

        """
        self.image_collection = pd.DataFrame
        self.args = self.get_args()
        self.night = self.set_night()
        self.extracted_data = None
        self.wsolution = None
        self.calibration_lamp = None
        self.wavelength_solution_obj = None

    def __call__(self):
        """Call method for the MainApp class

        This is equivalent to a main() function where all the logic and controls are implemented.

        Raises:
            NotImplementedError: For observing modes 2 and 3

        """
        self.organize_full_night()
        # print(len(self.night.sci_targets))
        for i in range(len(self.night.sci_targets)):
            science_object = self.night.sci_targets[i]
            # print(science_object)
            # print(self.night.sci_targets)
            process = Process(science_object, self.args)
            if self.args.obsmode == 0:
                if self.wavelength_solution_obj is None:
                    self.extracted_data, self.night.sci_targets[i] = process()
                    if isinstance(self.extracted_data, SciencePack):
                        wavelength_calibration = WavelengthCalibration(self.extracted_data,
                                                                       self.night.sci_targets[i],
                                                                       self.args)

                        self.wavelength_solution_obj = wavelength_calibration()

                        # self.night.set_night_wsolution(process.get_wsolution())
                        # self.night.set_night_calibration_lamp(process.get_calibration_lamp())
                    else:
                        log.error('No data was extracted from this target.')
                else:
                    if self.wavelength_solution_obj.check_compatibility(process.header):
                        log.debug(self.night.night_wsolution)
                        self.extracted_data, self.night.sci_targets[i] = process(extract_lamps=True)
                        if isinstance(self.extracted_data, SciencePack):
                            wavelength_calibration = WavelengthCalibration(self.extracted_data,
                                                                           self.night.sci_targets[i],
                                                                           self.args)

                            wavelength_calibration(self.wavelength_solution_obj)
                        else:
                            log.error('No data was extracted from this target.')
                    else:
                        log.debug('Incompatibility of solution to new data')
                        time.sleep(3)
                        # TODO(simon): complete this part
            elif self.args.obsmode == 1:
                self.extracted_data, self.night.sci_targets[i] = process()
                if self.extracted_data != []:
                    wavelength_calibration = WavelengthCalibration(self.extracted_data,
                                                                   self.night.sci_targets[i],
                                                                   self.args)

                    self.wavelength_solution_obj = wavelength_calibration()
                else:
                    log.error('No data was extracted from this target.')
            elif self.args.obsmode == 2:
                raise NotImplementedError
            elif self.args.obsmode == 3:
                raise NotImplementedError
            # else:
                # process = Process(self.night.source, science_object, self.args, self.night.night_wsolution)

    @staticmethod
    def get_args():
        """Handles the argparse library and returns the arguments

        Returns:
            An object that contains all the variables parsed through the argument system
            The possible arguments to be returned are:

            -p or --data-path: has to be the source directory, where the data (images) is.
                    the location is **self.source**
                    default value is ./
            -d or --proc-path: is the destination where all new files/data will be placed.
                    the location is self.destiny
                    default value is ./
            -s or --search-pattern: the pattern that matches the reduced data that will be processed.
                    the location is self.pattern
                    default value is fc_
            -m or --obs-mode: is one of the predefined observing modes and the options are:
                    0: One or more lamps taken during the beginning or end of the night, i.e. single
                    calibration to all data in that night
                    1: One or more lamps around every science exposure.
                    2: An ASCII file will be read. This file contains a matching of science target files
                    to respective calibration lamp file that will be used for calibration.
                    3: No lamp used, use sky-lines instead.
                    the location is self.mode
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

        Raises:
            In the case when -m or --obs-mode is set to 3 will requiere the name of file parsed with the -l or
            --lamp-file parameter an IOError is raised

        """
        leave = False
        parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                         description=textwrap.dedent(
                                             '''Extracts goodman spectra and does wavelength calibration.\n\n\
Supported Observing modes are:
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
                            default='fzh_',
                            type=str,
                            metavar='<Search Pattern>',
                            dest='pattern',
                            help="Pattern for matching the goodman's reduced data.")

        parser.add_argument('-m', '--obs-mode',
                            action='store',
                            default=0,
                            type=int,
                            metavar='<Observing Mode>',
                            dest='obsmode',
                            choices=[0, 1, 2, 3],
                            help='Defines the mode of matching lamps to science targets.')

        parser.add_argument('-r', '--reference-lamp',
                            action='store',
                            default='',
                            type=str,
                            metavar='<Lamp File>',
                            dest='lamp_all_night',
                            help="Name of reference lamp file for mode 0.\
                             If not present, the first one in the list will be selected")

        parser.add_argument('-l', '--lamp-file',
                            action='store',
                            default='lamps.txt',
                            type=str,
                            metavar='<Lamp File>',
                            dest='lamp_file',
                            help="Name of an ASCII file describing which science target\
                                uses which lamp. default <lamp.txt>")

        # parser.add_argument('-t', '--telescope',
        #                     action='store_true',
        #                     default=False,
        #                     dest='telescope',
        #                     help="Enables the <Telescope> mode i.e. it run sequentially,\
        #                         designed to use while observing at the telescope. Catches\
        #                          new files arriving to the <source> folder. (NI!)")

        parser.add_argument('-i', '--non-interactive',
                            action='store_false',
                            default=True,
                            dest='interactive_ws',
                            help="Interactive wavelength solution. Enabled by default.")

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

        parser.add_argument('--plots-enabled',
                            action='store_true',
                            default=False,
                            dest='plots_enabled',
                            help="Show plots of intermediate steps. For debugging only.")

        args = parser.parse_args()

        # there must be a more elegant way to do this
        root_path = os.path.realpath(__file__).split('/')
        root_path[-1] = ''
        root_full_path = '/'.join(root_path)
        reference_full_path = root_full_path + args.reference_dir
        if not os.path.isdir(reference_full_path):
            log.info("Reference files directory doesn't exist.")
            try:
                os.path.os.makedirs(reference_full_path)
                log.info('Reference Files Directory is: %s', reference_full_path)
                args.reference_dir = reference_full_path
            except OSError as err:
                log.error(err)
        else:
            args.reference_dir = reference_full_path

        if not os.path.isdir(args.source):
            leave = True
            log.error("Source Directory doesn't exist.")
        else:
            if args.source[-1] != '/':
                args.source += '/'
        if not os.path.isdir(args.destiny):
            leave = True
            log.error("Destination folder doesn't exist.")
            try:
                os.path.os.makedirs(args.destiny)
                log.info('Destination folder created: %s', args.destiny)
            except:
                pass
        else:
            if args.destiny[-1] != '/':
                args.destiny += '/'
        if args.obsmode == 2:
            # print(args.source + args.lamp_file)
            if not os.path.isfile(args.source + args.lamp_file):
                if args.lamp_file == 'lamps.txt':
                    leave = True
                    log.error("Default <lamp file> doesn't exist.")
                    log.error("Please define a <lamp file> using the flags -l or --lamp-file")
                    log.error("or make sure you entered the right observing mode.")
        if leave:
            parser.print_help()
            parser.exit("Leaving the Program.")

        return args

    def set_night(self):
        """Defines and initialize the 'night' class

        Uses information parsed by arguments to construct a table with the values of the keys specified within the
        code itself. A night object stores specific values regarding the night that is going to be processed. If the
        program is not going to be used while observing at the telescope it creates two lists of images, one for
        science and another for lamp files. Although the telescope mode is not fully developed yet, in case of it being
        selected it will return the night object without the list of images.

        Returns:
            new_night (class): A class that stores critical data of the night that will be processed and can be
            parsed to other methods.

        """
        keys = ['date', 'date-obs', 'obstype', 'object', 'exptime', 'ra', 'dec', 'grating']
        try:
            image_collection = ccd.ImageFileCollection(self.args.source, keys)
            self.image_collection = image_collection.summary.to_pandas()
        except ValueError as error:
            log.error('The images contain duplicated keywords')
            log.error('ValueError: %s', error)
            sys.exit(0)
        # type(self.image_collection)

        date = self.image_collection.date[0]
        new_night = Night(date, self.args)

        gratings_array = self.image_collection.grating.unique()
        new_night.set_gratings(gratings=gratings_array)

        # if self.args.telescope:
        #     new_night.is_telescope()
        #     log.info("Telescope Mode is not implemented yet...")
        #     return new_night
        # else:
        new_night.add_sci(self.image_collection.file[self.image_collection.obstype == 'OBJECT'])
        new_night.add_lamp(self.image_collection.file[self.image_collection.obstype == 'COMP'])
        return new_night

    def organize_full_night(self):
        """Organize the data according to the Observing Mode

        There are four observing modes defined by numbers in Python's style. From 0 to 3:

        In mode 0 one lamp is used to calibrate all the science targets of the night. As of December 2016 it picks
        up the first calibration lamp and uses it to find the wavelength calibration it doesn't discriminate if the lamp
        has good quality. It can also be parsed as an argument using -r or --reference-lamp

        In mode 1 one or more lamps are linked with a science target by matching them using two parameters. Distance
        in the sky equal or lower than 1e-3 degrees and a time difference of 300 seconds this is without the exposure
        time itself. For the sky distance calculation a flat sky is assumed.

        In mode 2 a text file is defined which correlates the science target with one or more lamps. Comments can be
        used by using a octothorp or hash (#) followed by one space. Not implemented yet.

        In mode 3 no sky lamp is used, instead the science target's spectrum will be calibrated using sky lines. Not
        implemented yet.

        It does not return anything but creates a ScienceObject instance and stores it in the night class. ScienceObject
        is one of its attributes.

        Raises:
            NotImplementedError: For mode 2 and 3.

        """
        self.print_spacers("Processing night %s" % self.night.date)

        if self.args.obsmode == 0:
            self.obsmode_zero()

        if self.args.obsmode == 1:
            self.obsmode_one()

        if self.args.obsmode == 2:
            self.obsmode_two()

        if self.args.obsmode == 3:
            self.obsmode_three()

        # science_object.print_all()
        # self.print_spacers(name)
        # print(self.night.sci)
        # print(self.night.lamp)

    def obsmode_zero(self):
        """Observing/Processing mode 0

        In mode 0 one lamp is used to calibrate all the science targets of the night. As of September 2016 it picks
        up the first calibration lamp and uses it to find the wavelength calibration it doesn't discriminate if the lamp
        has good quality.

        Notes:
            Although you can parse one lamp as a whole night lamp it is not recommended since there might be
            different gratings which would rise the need to give one lamp per grating and in this case it better
            to let the software choose the first in the list and assume there will be no bad lamps.
        """
        log.info("Observation mode 0")
        log.debug("One lamp for all targets.")

        for target in self.night.sci:
            index = self.image_collection[self.image_collection['file'] == target].index.tolist()[0]
            name = self.image_collection.object.iloc[index]
            grating = self.image_collection.grating.iloc[index]
            # print('Grating', grating)
            obs_time = self.image_collection['date-obs'][index]
            right_ascension, declination = self.ra_dec_to_deg(self.image_collection.ra.iloc[index],
                                                              self.image_collection.dec.iloc[index])
            science_object = ScienceObject(name, target, obs_time, right_ascension, declination, grating)

            comp_files = self.image_collection[(self.image_collection['grating'] == grating)
                                               & (self.image_collection['obstype'] == 'COMP')].index.tolist()
            # Need to define a better method for selecting the lamp
            # Now is just picking the first in the list
            if self.args.lamp_all_night is not '':
                lamp = self.args.lamp_all_night
                lamp_index = self.image_collection[self.image_collection['file'] == lamp].index.tolist()[0]
            else:
                lamp_index = comp_files[0]
                lamp = self.image_collection.file.iloc[lamp_index]
            log.debug("Lamp File: %s", lamp)

            lamp_name = self.image_collection.object.iloc[lamp_index]
            lamp_grating = self.image_collection.grating.iloc[lamp_index]
            # lamp_obs_time = self.image_collection['date-obs'][lamp_index]
            lamp_ra, lamp_dec = self.ra_dec_to_deg(self.image_collection.ra.iloc[lamp_index],
                                                   self.image_collection.ra.iloc[lamp_index])

            if lamp_grating == grating:
                science_object.add_lamp(lamp, lamp_name, lamp_ra, lamp_dec)
            # else:
                # log.info('Gratings do not match. Looking for another one.')

                # print(comp_files)
                # lamp_index = comp_files[0]
                # lamp_name = self.image_collection.object.iloc[lamp_index]
                lamp_grating = self.image_collection.grating.iloc[lamp_index]
                # lamp_obs_time = self.image_collection['date-obs'][lamp_index]
                # lamp_ra, lamp_dec = self.ra_dec_to_deg(self.image_collection.ra.iloc[lamp_index],
                #                                        self.image_collection.ra.iloc[lamp_index])
                # science_object.add_lamp(lamp, lamp_name, lamp_ra, lamp_dec)
            # science_object.print_all()
            self.night.add_sci_object(science_object)

        return

    def obsmode_one(self):
        """Observing/Processing mode 1

        In mode 1 one or more lamps are linked with a science target by matching them using two parameters. Distance
        in the sky equal or lower than 1e-3 degrees and a time difference of 300 seconds this is without the exposure
        time itself. For the sky distance calculation a flat sky is assumed.
        """
        log.info("Observation mode 1")
        log.debug("One or more lamps around the target")
        for target in self.night.sci:
            # Get basic data of the target
            index = self.image_collection[self.image_collection['file'] == target].index.tolist()[0]
            name = self.image_collection.object.iloc[index]
            obs_time = self.image_collection['date-obs'][index]
            exptime = self.image_collection.exptime.iloc[index]
            grating = self.image_collection.grating.iloc[index]
            right_ascension, declination = self.ra_dec_to_deg(self.image_collection.ra.iloc[index],
                                                              self.image_collection.dec.iloc[index])
            # Reformat some data of the target for comparison
            target_time = self.convert_time(obs_time)
            # Define ScienceObject object
            science_object = ScienceObject(name, target, obs_time, right_ascension, declination, grating)
            # Loop trough lamps to find a match for target
            for lamp in self.night.lamp:
                lamp_index = self.image_collection[self.image_collection['file'] == lamp].index.tolist()[0]
                lamp_name = self.image_collection.object.iloc[lamp_index]
                lamp_time = self.convert_time(self.image_collection['date-obs'][lamp_index])
                lamp_exptime = self.image_collection.exptime.iloc[lamp_index]
                lamp_grating = self.image_collection.grating.iloc[lamp_index]
                lamp_ra, lamp_dec = self.ra_dec_to_deg(self.image_collection.ra.iloc[lamp_index],
                                                       self.image_collection.dec.iloc[lamp_index])
                # print(lamp, lamp_name, lamp_ra, lamp_dec)
                if lamp_grating == grating:
                    # Since we are not doing astrometry here we assume the sky is flat
                    sky_distance = np.sqrt((lamp_ra - right_ascension) ** 2 + (lamp_dec - declination) ** 2)
                    if sky_distance <= 1e-3:
                        log.debug("Lamps by distance")
                        time_dif = abs(target_time - lamp_time) - abs(exptime + lamp_exptime)
                        if time_dif <= 300:
                            science_object.add_lamp(lamp, lamp_name, lamp_ra, lamp_dec)
                            # print(target,lamp,time_dif,exptime,lamp_exptime,sep=' : ')
                        else:
                            log.warning("Lamp within sky distance but too large time difference %s. Ignored.", time_dif)
                else:
                    log.info('Gratings do not match')
            # science_object.print_all()
            self.night.add_sci_object(science_object)
        return

    def obsmode_two(self):
        """Observing/Processing mode 2

        In mode 2 a text file is defined which correlates the science target with one or more lamps. Comments can be
        used by using a octothorp or hash (#) followed by one space. Not implemented yet.
        """
        log.info("Observation mode 2")
        log.debug("A text file defines the relation of lamps and science targets")
        log.debug(self.night.lamps_file)
        lamps_file_full = self.night.source + self.night.lamps_file
        log.debug(lamps_file_full)

        read_file = open(lamps_file_full)
        read_file = read_file.readlines()
        for i in range(len(read_file)):
            if read_file[i][0] != '#':
                read_file[i] = read_file[i].split()
                print(read_file[i])

    @staticmethod
    def obsmode_three():
        """Observing/Processing Mode 3

        In mode 3 no sky lamp is used, instead the science target's spectrum will be calibrated using sky lines.
        """
        log.info("Observation mode 3")
        log.debug("No Lamps. Use sky lines")
        raise NotImplementedError

    # Bunch of small functions

    @staticmethod
    def convert_time(in_time):
        """Converts time to seconds since epoch

        Args:
            in_time (str): time obtained from header's keyword DATE-OBS

        Returns:
            time in seconds since epoch

        """
        return time.mktime(time.strptime(in_time, "%Y-%m-%dT%H:%M:%S.%f"))

    @staticmethod
    def ra_dec_to_deg(right_ascension, declination):
        """Converts right ascension and declination to degrees

        Args:
            right_ascension (str): right ascension
            declination (str): declination

        Returns:
            right ascension and declination in degrees

        """
        right_ascension = right_ascension.split(":")
        declination = declination.split(":")
        # RIGHT ASCENTION conversion
        right_ascension_deg = (float(right_ascension[0])
                               + (float(right_ascension[1])
                                  + (float(right_ascension[2]) / 60.)) / 60.) * (360. / 24.)
        # DECLINATION conversion
        # print(declination)
        # sign = float(declination[0]) / abs(float(declination[0]))
        if float(declination[0]) == abs(float(declination[0])):
            sign = 1
        else:
            sign = -1
        declination_deg = sign * (abs(float(declination[0]))
                                  + (float(declination[1])
                                     + (float(declination[2]) / 60.)) / 60.)
        return right_ascension_deg, declination_deg

    @staticmethod
    def print_spacers(message):
        """Miscelaneous function to print uniform spacers

        Prints a spacer of 80 columns with  and 3 rows height. The first and last rows contains the symbol "="
        repeated 80 times. The middle row contains the message centered and the extremes has one single "=" symbol.
        The only functionality of this is aesthetic.

        Args:
            message (str): a message to be printed

        Returns:
            a True boolean

        """

        columns = 80
        if len(message) % 2 == 1 and int(columns) % 2 != 1:
            message += " "
        bar_length = int(columns)
        spacer_bar = "=" * bar_length
        blanks = bar_length - 2
        space_length = int((blanks - len(message)) / 2)
        message_bar = "=" + " " * space_length + message + " " * space_length + "="
        print(spacer_bar)

        print(message_bar)
        print(spacer_bar)
        return True

    @staticmethod
    def print_progress(current, total):
        """Prints the percentage of a progress

        It works for FOR loops, requires to know the full length of the loop. Prints to the standard output.

        Args:
            current (int): Current value in the range of the loop.
            total (int): The length of the loop.

        Returns:
            Nothing

        """
        if current == total:
            sys.stdout.write("Progress {:.2%}\n".format(1.0 * current / total))
        else:
            sys.stdout.write("Progress {:.2%}\r".format(1.0 * current / total))
        sys.stdout.flush()
        return


class Night(object):
    """Stores all data relevant to the night being processed

    Note:
        The night class stores the data relative to single observing night
        therefore this software works on a per-night basis
    """

    def __init__(self, date, args):
        """Initialize Night class

        The night class will store filename of science images as well as lamp images. It also stores the arguments and
        also the wavelength solution for the night.

        Args:
            date (str): Date of the night being processed
            args (object): argparse instance containing all the arguments the program was started with
        """
        # source, destiny, pattern, mode, lamps
        self.all = []
        self.date = date
        self.sci = []
        self.lamp = []
        self.args = args
        # self.args.source,
        # self.args.destiny,
        # self.args.pattern,
        # self.args.obsmode,
        # self.args.lamp_file
        # self.source = args.source
        # self.destiny = args.destiny
        # self.pattern = args.pattern
        # self.obsmode = args.obsmode
        # self.lamps_file = args.lamps_file
        self.sci_targets = []
        self.telescope = False
        self.night_wsolution = None
        self.night_calibration_lamp = None
        self.gratings = None

    def add_sci(self, in_sci):
        """Adds science object to list"""
        for new_sci in in_sci:
            if self.args.pattern == new_sci[0:len(self.args.pattern)]:
                self.sci.append(new_sci)
                self.all.append(new_sci)
            else:
                log.info('Image %s rejected because does not match pattern!', new_sci)
        # self.sci_targets = [[]] * len(self.sci)

    def add_lamp(self, in_lamp):
        """Adds lamp objects to list"""
        for new_lamp in in_lamp:
            if self.args.pattern == new_lamp[0:len(self.args.pattern)]:
                self.lamp.append(new_lamp)
                self.all.append(new_lamp)
            else:
                log.info('Image %s rejected because does not match pattern!', new_lamp)

    def add_sci_object(self, sci_obj):
        """Appends a ScienceObject to the class attribute sci_targets"""
        self.sci_targets.append(sci_obj)
        log.info("Added science object %s", sci_obj.name)

    def set_gratings(self, gratings):
        """Adds an array of the names of all the gratings observed in the night"""
        self.gratings = gratings

    def set_night_calibration_lamp(self, calibration_lamp):
        """Sets the filename of the calibration lamp as an attribute"""
        self.night_calibration_lamp = calibration_lamp

    def is_telescope(self):
        """Sets the operation mode as Telescope Mode"""
        self.telescope = True

    def set_night_wsolution(self, wsolution):
        """Sets the wavelength solution as a class attribute"""
        if wsolution is not None:
            self.night_wsolution = wsolution
            log.info('Night Solution Defined')
        else:
            log.error("Wavelength solution still can't be defined")


class ScienceObject(object):
    """class that defines a science object attributes

    Sci objects, for science, are the main targets which have lamps for
    calibration. Their atritutes are: the name, the file name, the
    observation time, right ascension and declination. The same information
    can be obtained for the lamps but they are treated as lists
    The list are created in a way that allows the use of elements' index
    as the correlator between atributes of the lamp.

    Attributes:
        name (str): science object name
        file_name (str): file name
        obs_time (str): observing time in the format yyyy-mm-ddThh:mm:ss.ss for instance 2016-03-20T23:54:15.96
        right_ascension (float): right ascension in degrees
        declination (float): declination in degrees
        lamp_count (int): lamps count
        lamp_file (list): every element is a string with the file name of the lamp
        lamp_type (list): every element is a string with the OBJECT value of the lamp i.e Cu, HgAr, etc
        lamp_ra (list): every element is a float with lamp's right ascension in degrees
        lamp_dec (list): every element is a float with lamp's declination in degrees

    """

    def __init__(self, name, file_name, obs_time, right_ascension, declination, grating):
        self.name = name
        self.file_name = file_name
        self.obs_time = obs_time
        self.right_ascension = right_ascension
        self.declination = declination
        self.lamp_count = 0
        self.lamp_file = []
        self.lamp_type = []
        self.lamp_ra = []
        self.lamp_dec = []
        self.grating = grating
        self.no_targets = 0

    def add_lamp(self, new_lamp, new_type, right_ascension, declination):
        """Adds a lamp to the science object

        Args:
            new_lamp (str): new lamp file name
            new_type (str): new lamp type
            right_ascension (float): right ascension in degrees
            declination (float): declination in degrees

        """
        self.lamp_file.append(new_lamp)
        self.lamp_type.append(new_type)
        self.lamp_ra.append(right_ascension)
        self.lamp_dec.append(declination)
        self.lamp_count = int(len(self.lamp_file))

    def update_no_targets(self, new_value=None, add_one=False):
        """Update number of spectra in an image

        An spectral image may contain multiple science target's spectra this method is set to update its number
        as a class attribute. There are two ways it can work. Add one to the existing number or set a new one.

        Args:
            new_value (int): New value for number of targets
            add_one (bool): If True increase the count by one.
        """
        if add_one:
            self.no_targets += 1
            log.debug('Number of targets is %s', self.no_targets)
        elif new_value is not None:
            self.no_targets = new_value
        else:
            log.debug('Nothing to do: new_value: %s add_one: %s', new_value, str(add_one))

    def print_all(self):
        """Prints all the relevant attributes of the object

        Note:
            this method is mainly used for development purposes

        """
        print(' ')
        log.info("Name: %s", self.name)
        log.info("File: %s", self.file_name)
        log.info("Obs-T: %s", self.obs_time)
        if self.lamp_count > 0:
            log.info("Lamp N: %s", self.lamp_count)
            for i in range(self.lamp_count):
                log.info("Lamp %s: %s", (i + 1), self.lamp_file[i])
                log.info("Type %s: %s", (i + 1), self.lamp_type[i])
        # time.sleep(10)


if __name__ == '__main__':
    MAIN_APP = MainApp()
    try:
        MAIN_APP()
    except KeyboardInterrupt:
        sys.exit(0)
