from __future__ import print_function
import os
from ccdproc import ImageFileCollection
import matplotlib.pyplot as plt
import time
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import re
import glob
import logging
import argparse
from data_classifier import DataClassifier
from night_organizer import NightOrganizer
from image_processor import ImageProcessor

FORMAT = '%(levelname)s: %(asctime)s: %(module)s: %(message)s'
DATE_FORMAT = '%m/%d/%Y %I:%M:%S%p'
logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt=DATE_FORMAT)
log = logging.getLogger('goodmanccd')


class MainApp(object):
    def __init__(self):
        """This method initalizes the MainApp class

        The main task of this method is to call the get_args function that returns an argparse object.
        The arguments will be obtained here and they will be available for all execution of the program.
        Other attributes will be initilized as None.
        """
        self.args = self.get_args()
        self.data_container = None
        self.full_path = None
        self.instrument = None
        self.technique = None
    
    def __call__(self):
        """Call method for MainApp

        From the arguments this method finds the raw_path attribute and checks its contents for the existance of
        files containing the '.fits' string. If there is none it will assume every item is a different data directory and they
        will be treated independently. If there are '.fits' files the program will assume is a single data directory.
        Any subdirectory will be ignored.

        """
        folders = glob.glob(re.sub('//', '/', '/'.join(self.args.raw_path.split('/') + ['*'])))
        # print(re.sub('//', '/', '/'.join(self.args.raw_path.split('/') + ['*'])))
        if any('.fits' in item for item in folders):
            folders = [self.args.raw_path]
        for data_folder in folders:
            # check start
            self.args.raw_path = data_folder
            if self.args.red_path == './RED' or len(folders) > 1:
                log.info('No special reduced data path defined. Proceeding with defaults.')
                if self.args.raw_path not in self.args.red_path:
                    self.args.red_path = re.sub('//', '/', '/'.join(self.args.raw_path.split('/') + ['RED']))
                    # print(self.args.red_path)
            if os.path.isdir(self.args.red_path):
                if os.listdir(self.args.red_path) != []:
                    log.warning('Reduced Data Path is not empty')
                    if self.args.auto_clean:
                        for _file in os.listdir(self.args.red_path):
                            os.unlink(self.args.red_path + '/' + _file)
                        log.info('Cleaned Reduced data directory: ' + self.args.red_path)
                    else:
                        log.error('Please clean the reduced data folder or use --auto-clean')
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
            night_sorter = DataClassifier(self.args)
            night_sorter()
            self.instrument = night_sorter.instrument
            self.technique = night_sorter.technique
            # print(night_sorter.nights_dict)
            for night in night_sorter.nights_dict:
                # print(night_sorter.nights_dict[night])
                night_organizer = NightOrganizer(args=self.args, night_dict=night_sorter.nights_dict[night])
                self.data_container = night_organizer()
                if self.data_container is False or self.data_container is None:
                    log.error('Discarding night ' + str(night))
                    break
                process_images = ImageProcessor(self.args, self.data_container)
                process_images()

    @staticmethod
    def get_args():
        """Get command line arguments.

        Returns:
            args (object): Arparse object. Contains all the arguments as attributes

        """
        # Parsing Arguments ---
        parser = argparse.ArgumentParser(description="PyGoodman CCD Reduction - CCD reductions for "
                                                     "Goodman spectroscopic data")

        # parser.add_argument('-c', '--clean',
        #                     action='store_true',
        #                     help="Clean cosmic rays from science data.")

        # # removed because is not working properly
        # parser.add_argument('-s', '--slit', action='store_true',
        #                     help="Find slit edge to make an additional trimming (Maintainer: Not recommended for now).")

        # TODO (simon): Add argument to use calibration data from other day

        # remove saturated data
        parser.add_argument('--remove-saturated',
                            action='store_true',
                            dest='remove_saturated',
                            help="Remove images above saturation level")

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
                            help="Path to raw data (e.g. /home/jamesbond/soardata/).")

        parser.add_argument('--red-path',
                            action='store',
                            metavar='red_path',
                            type=str,
                            default='./RED',
                            help="Full path to reduced data (e.g /home/jamesbond/soardata/RED/).")

        args = parser.parse_args()
        if os.path.isdir(args.raw_path):
            args.raw_path = os.path.abspath(args.raw_path)
            log.debug(os.path.abspath(args.raw_path))
        else:
            parser.print_help()
            parser.exit("Raw data folder doesn't exist")

        return args


if __name__ == '__main__':
    main_app = MainApp()
    main_app()
