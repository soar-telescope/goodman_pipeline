#!/usr/bin/python2.7
"""Pipeline for GOODMAN spectra Extraction.

This program finds reduced images, i.e. trimmed, bias subtracted, flat fielded, etc. that match the <pattern>
in the source folder, then classify them in two groups: Science or Lamps. For science images, finds the spectrum
or spectra and traces it doing some fit.
Simon Torres 2016-06-28

"""

import sys
import os
import numpy as np
import time
import textwrap
import ccdproc as ccd
import pandas as pd
import argparse
import logging as log
import warnings
from process import Process
from wavelengthCalibration import WavelengthCalibration

warnings.filterwarnings('ignore')
FORMAT = '%(levelname)s:%(filename)s:%(module)s: 	%(message)s'
log.basicConfig(level=log.DEBUG, format=FORMAT)

__author__ = 'Simon Torres'
__date__ = '2016-06-28'
__version__ = "0.1"
__email__ = "storres@ctio.noao.edu"
__status__ = "Development"


class MainApp:
    """Defines and intialize all important variables for processing the data

    Args:
        It doesn't take any arguments

    """
    def __init__(self):
        self.ic = pd.DataFrame
        self.args = self.get_args()
        self.night = self.set_night()
        self.extracted_data = None
        self.wsolution = None
        self.calibration_lamp = None
        self.wavelength_solution_obj = None

    def __call__(self, *args, **kwargs):
        self.organize_full_night()
        # print(len(self.night.sci_targets))
        for i in range(len(self.night.sci_targets)):
            science_object = self.night.sci_targets[i]
            # print(science_object)
            # print(self.night.sci_targets)
            process = Process(self.night.source, science_object, self.args)
            if self.args.obsmode == 0:
                if self.wavelength_solution_obj is None:
                    self.extracted_data, self.night.sci_targets[i] = process()
                    wavelength_calibration = WavelengthCalibration(self.night.source,
                                                                   self.extracted_data,
                                                                   self.night.sci_targets[i],
                                                                   self.args)
                    self.wavelength_solution_obj = wavelength_calibration()

                    # self.night.set_night_wsolution(process.get_wsolution())
                    # self.night.set_night_calibration_lamp(process.get_calibration_lamp())
                else:
                    log.debug(self.night.night_wsolution)
                    self.extracted_data, self.night.sci_targets[i] = process(extract_lamps=False)
                    wavelength_calibration = WavelengthCalibration(self.night.source,
                                                                   self.extracted_data,
                                                                   self.night.sci_targets[i],
                                                                   self.args)
                    wavelength_calibration(self.wavelength_solution_obj)
            elif self.args.obsmode == 1:
                self.extracted_data, self.night.sci_targets[i] = process()
                wavelength_calibration = WavelengthCalibration(self.night.source,
                                                               self.extracted_data,
                                                               self.night.sci_targets[i],
                                                               self.args)
                self.wavelength_solution_obj = wavelength_calibration()
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
            -l or --lamp-file: Name of the ASCII file that contains the relation between science files and lamp
                    files. An example is depicted below. Note that the lamps can be repeated.
                        #example of how the file should look
                        science_target_01.fits lamp_001.fits
                        science_target_02.fits lamp_001.fits
                        science_target_03.fits lamp_002.fits
                    the location is self.lamp_file
                    default value is lamps.txt

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
    <1>: one or more lamps around science exposure.\n\
    <2>: ASCII file describing which science target uses which lamp.\n\
    <3>: No lamps. Uses the sky lines
    '''))

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
                            default='fc_',
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

        parser.add_argument('-t', '--telescope',
                            action='store_true',
                            default=False,
                            dest='telescope',
                            help="Enables the <Telescope> mode i.e. it run sequentially,\
                                designed to use while observing at the telescope. Catches\
                                 new files arriving to the <source> folder. (NI!)")

        parser.add_argument('-i', '--interactive',
                            action='store_true',
                            default=False,
                            dest='interactive_ws',
                            help="Interactive wavelength solution. Disabled by default.")

        parser.add_argument('-o', '--output-prefix',
                            action='store',
                            default='g_',
                            dest='output_prefix',
                            help="Prefix to add to calibrated spectrum.")

        args = parser.parse_args()
        if not os.path.isdir(args.source):
            leave = True
            log.error("Source Directory doesn't exist.")
        else:
            if args.source[-1] != '/':
                args.source += '/'
        if not os.path.isdir(args.destiny):
            leave = True
            log.error("Destination folder doesn't exist.")
        else:
            if args.destiny[-1] != '/':
                args.destiny += '/'
        if args.obsmode == 2:
            print(args.source + args.lamp_file)
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
        keys = ['date', 'date-obs', 'obstype', 'object', 'exptime', 'ra', 'dec']
        image_collection = ccd.ImageFileCollection(self.args.source, keys)
        self.ic = image_collection.summary.to_pandas()
        # type(self.ic)
        date = self.ic.date[0]
        new_night = Night(date,
                          self.args.source,
                          self.args.destiny,
                          self.args.pattern,
                          self.args.obsmode,
                          self.args.lamp_file)
        if self.args.telescope:
            new_night.is_telescope()
            log.info("Telescope Mode is not implemented yet...")
            return new_night
        else:
            new_night.add_sci(self.ic.file[self.ic.obstype == 'OBJECT'])
            new_night.add_lamp(self.ic.file[self.ic.obstype == 'COMP'])
        print new_night
        return new_night

    def organize_full_night(self):
        """Organize the data according to the Observing Mode

        There are four observing modes defined by numbers in Python's style. From 0 to 3:

        In mode 0 one lamp is used to calibrate all the science targets of the night. As of September 2016 it picks
        up the first calibration lamp and uses it to find the wavelength calibration it doesn't discriminate if the lamp
        has good quality.

        In mode 1 one or more lamps are linked with a science target by matching them using two parameters. Distance
        in the sky equal or lower than 1e-3 degrees and a time difference of 300 seconds this is without the exposure
        time itself. For the sky distance calculation a flat sky is assumed.

        In mode 2 a text file is defined which correlates the science target with one or more lamps. Comments can be
        used by using a octothorp or hash (#) followed by one space. Not implemented yet.

        In mode 3 no sky lamp is used, instead the science target's spectrum will be calibrated using sky lines.

        Returns:
            Nothing. But creates a ScienceObject object and stores it in the night class. ScienceObject is one of its
            attributes.

        Raises:
            NotImplementedError: For mode 3.

        """
        self.print_spacers("Processing night %s" % self.night.date)
        if self.night.obsmode == 0:
            log.info("Observation mode 0")
            log.debug("One lamp for all targets.")
            """
            Need to define a better method for selecting the lamp
            Now is just picking the first in the list
            """
            if self.args.lamp_all_night is not '':
                lamp = self.args.lamp_all_night
            else:
                lamp = self.night.lamp[0]
            log.debug("Lamp File: " + lamp)
            lamp_index = self.ic[self.ic['file'] == lamp].index.tolist()[0]
            lamp_name = self.ic.object.iloc[lamp_index]
            # lamp_obs_time = self.ic['date-obs'][lamp_index]
            lamp_ra, lamp_dec = self.ra_dec_to_deg(self.ic.ra.iloc[lamp_index], self.ic.ra.iloc[lamp_index])
            for target in self.night.sci:
                index = self.ic[self.ic['file'] == target].index.tolist()[0]
                name = self.ic.object.iloc[index]
                obs_time = self.ic['date-obs'][index]
                ra, dec = self.ra_dec_to_deg(self.ic.ra.iloc[index], self.ic.dec.iloc[index])
                science_object = ScienceObject(name, target, obs_time, ra, dec)
                science_object.add_lamp(lamp, lamp_name, lamp_ra, lamp_dec)
                self.night.add_sci_object(science_object)
            return
        if self.night.obsmode == 1:
            log.info("Observation mode 1")
            log.debug("One or more lamps around the target")
            for target in self.night.sci:
                """Get basic data of the target"""
                index = self.ic[self.ic['file'] == target].index.tolist()[0]
                name = self.ic.object.iloc[index]
                obs_time = self.ic['date-obs'][index]
                exptime = self.ic.exptime.iloc[index]
                ra, dec = self.ra_dec_to_deg(self.ic.ra.iloc[index], self.ic.dec.iloc[index])
                """Reformat some data of the target for comparison"""
                target_time = self.convert_time(obs_time)
                """Define ScienceObject object"""
                science_object = ScienceObject(name, target, obs_time, ra, dec)
                """Loop trough lamps to find a match for target"""
                for lamp in self.night.lamp:
                    lamp_index = self.ic[self.ic['file'] == lamp].index.tolist()[0]
                    lamp_name = self.ic.object.iloc[lamp_index]
                    lamp_time = self.convert_time(self.ic['date-obs'][lamp_index])
                    lamp_exptime = self.ic.exptime.iloc[lamp_index]
                    lamp_ra, lamp_dec = self.ra_dec_to_deg(self.ic.ra.iloc[lamp_index], self.ic.dec.iloc[lamp_index])
                    # print(lamp, lamp_name, lamp_ra, lamp_dec)
                    """Since we are not doing astrometry here we assume the sky is flat"""
                    sky_distance = np.sqrt((lamp_ra - ra) ** 2 + (lamp_dec - dec) ** 2)
                    if sky_distance <= 1e-3:
                        log.debug("Lamps by distance")
                        time_dif = abs(target_time - lamp_time) - abs(exptime + lamp_exptime)
                        if time_dif <= 300:
                            science_object.add_lamp(lamp, lamp_name, lamp_ra, lamp_dec)
                            # print(target,lamp,time_dif,exptime,lamp_exptime,sep=' : ')
                        else:
                            log.warning("Lamp within sky distance but too large time difference. Ignored.")
                self.night.add_sci_object(science_object)
            return

        if self.night.obsmode == 2:
            log.info("Observation mode 2")
            log.debug("A text file defines the relation of lamps and science targets")
            log.debug(self.night.lamps_file)
            lamps_file_full = self.night.source + self.night.lamps_file
            log.debug(lamps_file_full)

            ff = open(lamps_file_full)
            ff = ff.readlines()
            for i in range(len(ff)):
                if ff[i][0] != '#':
                    ff[i] = ff[i].split()
                    print(ff[i])

        if self.night.obsmode == 3:
            log.info("Observation mode 3")
            log.debug("No Lamps. Use sky lines")
            raise NotImplementedError

        # science_object.print_all()
        # self.print_spacers(name)
        # print(self.night.sci)
        # print(self.night.lamp)

    """Bunch of small functions"""

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
    def ra_dec_to_deg(ra, dec):
        """Converts right ascension and declination to degrees

        Args:
            ra (str): right ascension
            dec (str): declination

        Returns:
            right ascension and declination in degrees

        """
        ra = ra.split(":")
        dec = dec.split(":")
        # RIGHT ASCENTION conversion
        ra_deg = (float(ra[0]) + (float(ra[1]) + (float(ra[2]) / 60.)) / 60.) * (360. / 24.)
        # DECLINATION conversion
        sign = float(dec[0]) / abs(float(dec[0]))
        dec_deg = sign * (abs(float(dec[0])) + (float(dec[1]) + (float(dec[2]) / 60.)) / 60.)
        return ra_deg, dec_deg

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
        bar = "=" * bar_length
        blanks = bar_length - 2
        space_length = int((blanks - len(message)) / 2)
        message_bar = "=" + " " * space_length + message + " " * space_length + "="
        print(bar)

        print(message_bar)
        print(bar)
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


class Night:
    """Stores all data relevant to the night being processed

    Note:
        The night class stores the data relative to single observing night
        therefore this software works on a per-night basis
    """

    def __init__(self, date, source, destiny, pattern, mode, lamps):
        self.all = []
        self.date = date
        self.sci = []
        self.lamp = []
        self.source = source
        self.destiny = destiny
        self.pattern = pattern
        self.obsmode = mode
        self.lamps_file = lamps
        self.sci_targets = []
        self.telescope = False
        self.night_wsolution = None
        self.night_calibration_lamp = None

    def add_sci(self, in_sci):
        """Adds science object to list"""
        for new_sci in in_sci:
            if self.pattern == new_sci[0:len(self.pattern)]:
                self.sci.append(new_sci)
                self.all.append(new_sci)
            else:
                log.info('Image %s rejected because does not match pattern!', new_sci)
        # self.sci_targets = [[]] * len(self.sci)

    def add_lamp(self, in_lamp):
        """Adds lamp objects to list"""
        for new_lamp in in_lamp:
            if self.pattern == new_lamp[0:len(self.pattern)]:
                self.lamp.append(new_lamp)
                self.all.append(new_lamp)
            else:
                log.info('Image %s rejected because does not match pattern!', new_lamp)


    def add_sci_object(self, sci_obj):
        """Appends a ScienceObject to the class attribute sci_targets"""
        self.sci_targets.append(sci_obj)
        log.info("Added science object %s" % sci_obj.name)

    def set_night_calibration_lamp(self, calibration_lamp):
        self.night_calibration_lamp = calibration_lamp

    def is_telescope(self):
        """Sets the operation mode as Telescope Mode"""
        self.telescope = True

    def set_night_wsolution(self, wsolution):
        if wsolution is not None:
            self.night_wsolution = wsolution
            log.info('Night Solution Defined')
        else:
            log.error("Wavelength solution still can't be defined")


class ScienceObject:
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
        ra (float): right ascension in degrees
        dec (float): declination in degrees
        lamp_count (int): lamps count
        lamp_file (list): every element is a string with the file name of the lamp
        lamp_type (list): every element is a string with the OBJECT value of the lamp i.e Cu, HgAr, etc
        lamp_ra (list): every element is a float with lamp's right ascension in degrees
        lamp_dec (list): every element is a float with lamp's declination in degrees

    """

    def __init__(self, name, file_name, obs_time, ra, dec):
        self.name = name
        self.file_name = file_name
        self.obs_time = obs_time
        self.ra = ra
        self.dec = dec
        self.lamp_count = 0
        self.lamp_file = []
        self.lamp_type = []
        self.lamp_ra = []
        self.lamp_dec = []

    def add_lamp(self, new_lamp, new_type, ra, dec):
        """Adds a lamp to the science object

        Args:
            new_lamp (str): new lamp file name
            new_type (str): new lamp type
            ra (float): right ascension in degrees
            dec (float): declination in degrees

        """
        self.lamp_file.append(new_lamp)
        self.lamp_type.append(new_type)
        self.lamp_ra.append(ra)
        self.lamp_dec.append(dec)
        self.lamp_count = int(len(self.lamp_file))

    def print_all(self):
        """Prints all the relevant attributes of the object

        Note:
            this method is mainly used for development purposes

        """
        log.info("Name: %s" % self.name)
        log.info("File: %s" % self.file_name)
        log.info("Obs-T: %s" % self.obs_time)
        if self.lamp_count > 0:
            log.info("Lamp N: %s" % self.lamp_count)
            for i in range(self.lamp_count):
                log.info("Lamp %s: %s" % ((i + 1), self.lamp_file[i]))
                log.info("Type %s: %s" % ((i + 1), self.lamp_type[i]))


if __name__ == '__main__':
    app = MainApp()
    app()
