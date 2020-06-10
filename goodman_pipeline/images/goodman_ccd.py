from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .data_classifier import DataClassifier
from .night_organizer import NightOrganizer
from .image_processor import ImageProcessor

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
    arg_log = logging.getLogger()

    parser = argparse.ArgumentParser(
        description="Goodman CCD Reduction - CCD reductions for Goodman "
                    "spectroscopic data."
                    "\nPipeline Version: {:s}".format(__version__))

    parser.add_argument('--auto-clean',
                        action='store_true',
                        dest='auto_clean',
                        help="Automatically clean reduced data directory")

    parser.add_argument('--cosmic',
                        action='store',
                        dest='clean_cosmic',
                        default='default',
                        choices=['default', 'dcr', 'lacosmic', 'none'],
                        metavar='<method>',
                        help="Clean cosmic rays from all data. Options are: "
                             "'default', 'dcr', 'lacosmic' or 'none'. "
                             "See manual for full description of dcr.")

    parser.add_argument('--combine',
                        action='store_true',
                        dest='combine',
                        help="Combine compatible data (experimental)")

    parser.add_argument('--dcr-par-dir',
                        action='store',
                        default='data/params',
                        metavar='<dcr.par_directory>',
                        dest='dcr_par_dir',
                        help="Directory of default dcr.par file")

    parser.add_argument('--debug',
                        action='store_true',
                        dest='debug_mode',
                        help="Show detailed information of the process.")

    # TODO (simon): Add argument to use calibration data from other day

    parser.add_argument('--flat-normalize',
                        action='store',
                        default='simple',
                        type=str,
                        metavar='<normalization_method>',
                        dest='flat_normalize',
                        choices=['mean', 'simple', 'full'],
                        help="Choose a method to normalize the master flat for "
                             "spectroscopy. Choices are: 'mean',"
                             "'simple' (model) and 'full' "
                             "(fits model to each line). Default 'simple'")

    parser.add_argument('--flat-norm-order',
                        action='store',
                        default=15,
                        type=int,
                        metavar='<order>',
                        dest='norm_order',
                        help='Defines the order of the model to be fitted. '
                             'Default to 15')

    parser.add_argument('--ignore-bias',
                        action='store_true',
                        dest='ignore_bias',
                        help="Ignore bias correction")

    parser.add_argument('--ignore-flats',
                        action='store_true',
                        dest='ignore_flats',
                        help="Ignore flat field correction")

    parser.add_argument('--skip-slit-trim',
                        action='store_true',
                        dest='skip_slit_trim',
                        help="Do not apply slit trimming")

    parser.add_argument('--keep-cosmic-files',
                        action='store_true',
                        dest='keep_cosmic_files',
                        help="After cleaning cosmic rays with 'dcr', do not "
                             "remove the input file and the cosmic rays file. "
                             "For 'lacosmic' it saves the mask to a fits file.")

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

    parser.add_argument('--saturation_threshold',
                        action='store',
                        default=1.,
                        dest='saturation_threshold',
                        metavar='<value>',
                        help="Maximum percent of pixels above saturation_threshold "
                             "threshold. Default 1 percent.")

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
        self.args = None
        self.log.debug('Initializing DataClassifier instance')
        self.data_classifier = DataClassifier()
        self.data_container = None
        self.full_path = None
        # self.instrument = None
        # self.technique = None
        self._pipeline_version = __version__

    def __call__(self, args=None):
        """Call method for MainApp

        From the arguments this method finds the raw_path attribute and checks
        its contents for the existence of files containing the '.fits' string.
        If there is none it will assume every item is a different data directory
        and they will be treated independently. If there are '.fits' files the
        program will assume is a single data directory.
        Any subdirectory will be ignored.

        Args:
            args (list): a list of arguments and values, this is useful when you
                want to import the class.
        """

        if args is None:
            self.args = get_args()
        else:
            self.args = args

        if self.args.show_version:
            print("Goodman HTS Pipeline {:s}".format(__version__))
            sys.exit(0)

        if not self._check_args():
            sys.exit()

        # Division point for the future implementation of *live reduction mode*

        try:
            self.log.debug('Calling data_classifier '
                           'Instance of DataClassifier')
            self.data_classifier(raw_path=self.args.raw_path)
        except AttributeError as error:

            if 'instconf' in error.args[0]:
                self.log.error("Card '{:s}' ".format('instconf'.upper()) +
                               "not found inside headers. "
                               "This keyword contains which Goodman "
                               "Camera was used for this observation. "
                               "Please add it manually (Blue/Red) and "
                               "run again. Leaving the program now.")
                sys.exit(1)

            if 'wavmode' in error.args[0]:
                self.log.error(
                    "Card '{:s}' ".format('wavmode'.upper()) +
                    "not found inside headers. This keyword contains what "
                    "is the Goodman Wavelength Configuration that was used"
                    " for this observation. Please add it manually (see "
                    "https://github.com/soar-telescope/goodman_pipeline/blob/development/goodman_modes.md) "
                    "and run again. Leaving the program now.")
                sys.exit(1)

            self.log.error(error)
            self.log.error('Empty or Invalid data directory:'
                           '{:s}'.format(self.args.raw_path))

        # print(self.data_classifier.nights_dict)
        for night in self.data_classifier.nights_dict:
            nd = self.data_classifier.nights_dict[night]
            self.log.debug('Initializing night organizer procedure')
            night_organizer = NightOrganizer(
                full_path=nd['full_path'],
                instrument=nd['instrument'],
                technique=nd['technique'],
                ignore_bias=self.args.ignore_bias,
                ignore_flats=self.args.ignore_flats)

            self.log.debug('Calling night organizer procedure')
            try:
                data_container_list = night_organizer()
            except IOError as error:
                self.log.critical(error)
                sys.exit(1)
            for self.data_container in data_container_list:
                # print(self.data_container)
                if self.data_container is None or \
                        self.data_container.is_empty:
                    self.log.info("Data container is empty")
                    self.log.error('Discarding night {:s}'
                                   ''.format(str(night)))
                else:
                    self.log.debug("Initializing image processing "
                                   "procedure")
                    process_images = ImageProcessor(
                        args=self.args,
                        data_container=self.data_container)
                    self.log.debug("Calling image processing procedure.")
                    process_images()

    def _check_args(self):
        """Perform checks to arguments

        Returns:
            False in case anything fails, True if everything succeeds.

        """
        self.log.debug("Starting argument checks.")
        # check if raw folder exist
        if os.path.isdir(self.args.raw_path):
            self.args.raw_path = os.path.abspath(self.args.raw_path)
            self.log.debug("Raw data folder \"{:s}\" exists."
                           "".format(os.path.abspath(self.args.raw_path)))
        else:
            log.critical("Raw data folder \"{:s}\" doesn't "
                         "exist".format(self.args.raw_path))
            return False

        raw_folder_content = glob.glob(os.path.join(self.args.raw_path,
                                                    '*fits'))

        if any(['.fits' in _item for _item in raw_folder_content]):
            self.log.info("Found {:d} fits files in {:s}."
                          "".format(len(raw_folder_content),
                                    self.args.raw_path))
        else:
            self.log.critical("Raw data folder \"{:s}\" does not contain any "
                              "'*.fits' files.".format(self.args.raw_path))
            return False

        # check start
        if self.args.red_path == './RED':

            self.log.info('No special reduced data path defined. '
                          'Proceeding with defaults.')

            if self.args.raw_path not in self.args.red_path:
                self.args.red_path = os.path.join(self.args.raw_path, 'RED')
                self.log.debug("Folder for reduced data defined to: {:s}"
                               "".format(self.args.red_path))

        if os.path.isdir(self.args.red_path):
            try:
                _directory_content = os.listdir(self.args.red_path)
            except PermissionError as error:
                self.log.debug(error)
                self.log.critical("Unable to read on directory {:s}"
                                  "".format(self.args.red_path))
                sys.exit()

            if _directory_content != []:
                self.log.warning('Folder for reduced data is not empty')
                if self.args.auto_clean:
                    self.log.info("--auto-clean is set")
                    for _file in os.listdir(self.args.red_path):
                        try:
                            _to_delete = os.path.join(self.args.red_path, _file)

                            os.unlink(_to_delete)

                            self.log.debug("Cleanup: Deleting {:s}"
                                           "".format(_to_delete))
                        except OSError as error:
                            self.log.error(
                                'OSError: {:s}'.format(str(error)))
                            self.log.warning('Removing Directory '
                                             '{:s}'.format(_file))

                            try:
                                shutil.rmtree(os.path.join(self.args.red_path,
                                                           _file))
                            except PermissionError as error:
                                self.log.debug(error)
                                self.log.critical("Unable to delete files on "
                                                  "directory {:s}"
                                                  "".format(self.args.red_path))
                                self.log.info("Please check permissions")
                                return False

                    self.log.info('Cleaned Reduced data directory:'
                                  ' {:s}'.format(self.args.red_path))
                else:
                    self.log.error('Please define another folder using '
                                   '--red-path or or use --auto-clean to '
                                   'automatically wipe the current one.')
                    return False
            else:
                self.log.info("Folder for reduced data is empty (good).")
                self.args.red_path = os.path.abspath(self.args.red_path)
                self.log.debug(os.path.abspath(self.args.red_path))
        else:
            try:
                self.log.warning("Reduction folder doesn't exist.")
                os.mkdir(os.path.abspath(self.args.red_path))
                self.log.info('Created directory for reduced data.')
                self.log.info("New directory: {:s}".format(
                    os.path.abspath(self.args.red_path)))
            except OSError as error:
                self.log.error(error)
                # check ends

        # updated full path for default dcr.par file. If it doesn't exist it will
        # create an empty one.
        # print(sys.modules['goodman_pipeline'].__file__)
        if not os.path.isabs(self.args.dcr_par_dir):
            dcr_par_full_path = os.path.join(
                os.path.dirname(sys.modules['goodman_pipeline'].__file__),
                self.args.dcr_par_dir)
        else:
            dcr_par_full_path = self.args.dcr_par_dir
        if not os.path.isdir(dcr_par_full_path) and \
                        self.args.dcr_par_dir != 'data/params':
            self.log.info("dcr.par location {:s} doesn't exist."
                          "".format(dcr_par_full_path))
            try:
                os.path.os.makedirs(dcr_par_full_path)
                self.log.info('Created dcr.par empty directory: '
                              '{:s}'.format(dcr_par_full_path))
                self.args.dcr_par_dir = dcr_par_full_path
            except OSError as err:
                self.log.error(err)
        else:
            self.args.dcr_par_dir = dcr_par_full_path
            # print_default_args(args)
        return True


if __name__ == '__main__':  # pragma: no cover
    main_app = MainApp()
    main_app()
