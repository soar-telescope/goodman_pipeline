from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import logging
import os
import glob
import re
import random

import numpy as np

from ccdproc import ImageFileCollection
from astropy.io import fits

from .core import fix_duplicated_keywords

log = logging.getLogger('goodmanccd.dataclassifier')


class DataClassifier(object):
    """Class Definition

    Data classifier is intended to define the camera that is being used and the
    technique in use. This will be used later to make important decisions
    regarding the process to be used.

    """

    def __init__(self, args):
        """Initialization method for the DataClassifier class

        The general arguments of the program are parsed and become part of the
        class attributes. The rest of attributes are initialized as None.

        Args:
            args (object): Argparse object
        """
        self.args = args
        self.nights_dict = None
        self.instrument = None
        self.image_collection = None
        self.objects_collection = None
        self.technique = None

    def __call__(self):
        """Call method for the DataClassifier class

        This method call specific method that define all the attributes of the
        class. The intention is to define the instrument and technique in use.

        """
        self.nights_dict = {}
        log.debug('Raw path: ' + self.args.raw_path)
        self.get_instrument(self.args.raw_path)
        log.info('Instrument: ' + self.instrument + ' Camera')
        self.get_obs_technique()
        log.info('Observing Technique: ' + self.technique)
        if self.instrument is not None and self.technique is not None:
            # folder name is used as key for the dictionary
            night = self.args.raw_path.split('/')[-1]
            self.nights_dict[night] = {'full_path': self.args.raw_path,
                                       'instrument': self.instrument,
                                       'technique': self.technique}
        else:
            log.error('Failed to determine Instrument or Technique '
                      'for the night: {:s}'.format(self.args.raw_path))

    def get_instrument(self, night_folder):
        """Identify Goodman's Camera

        Goodman has two camera, a blue one and another optimized for redder
        wavelength. Their headers are differents so this methods uses that to
        discover which camera the data belong to. The red camera has an specific
        keyword that says which camera is but the blue does not.

        Args:
            night_folder (str): The full path for the raw data location

        The result is stored as an attribute of the class.

        """
        while True:
            try:
                ifc = ImageFileCollection(night_folder)
                self.image_collection = ifc.summary.to_pandas()

                self.objects_collection = self.image_collection[
                    self.image_collection.obstype != 'ZERO']

                # print(len(self.objects_collection))
                if len(self.objects_collection) > 0:

                    indexes = self.objects_collection.index.tolist()
                    index = random.choice(indexes)
                    # print('Index: ' + str(index))
                    try:
                        # print(self.objects_collection.instconf[index])
                        self.instrument = \
                            self.objects_collection.instconf[index]

                    except AttributeError as error:
                        log.error(error)
                        # print(self.objects_collection.file[index])
                        self.instrument = 'Blue'
                else:
                    log.error('There is no useful data in this folder.')
            except ValueError as error:
                if 'Inconsistent data column lengths' in str(error):

                    log.error('There are duplicated keywords in the headers. '
                              'Fix it first!')

                    fix_duplicated_keywords(night_folder)
                    continue
                else:
                    log.error('Unknown Error: ' + str(error))
            break

    def get_obs_technique(self):
        """Identify if the data is Imaging or Spectroscopy

        Besides the fact there are two cameras there are two observational
        techniques. Imaging and Spectroscopy. The red camera has an specific
        keyword that contains that information but the blue does not.

        The result is stored as an attribute of the class.

        """

        if self.instrument == 'Red':
            wavmodes = self.objects_collection.wavmode.unique()
            if len(wavmodes) == 1 and wavmodes[0] == 'Imaging':
                self.technique = 'Imaging'
                log.info('Detected Imaging Data from RED Camera')
            elif 'Imaging' in wavmodes and len(wavmodes) > 1:
                log.error('There are mixed observation techniques this night. '
                          'Please classify your data')
                self.technique = 'Unknown'
            else:
                self.technique = 'Spectroscopy'
                log.info('Detected Spectroscopy Data from RED Camera')
        elif self.instrument == 'Blue':
            self.remove_conflictive_keywords()
            gratings = self.objects_collection.grating.unique()
            cam_targ = self.objects_collection.cam_targ.unique()
            # print(gratings)
            if ['<NO GRATING>'] not in gratings:
                self.technique = 'Spectroscopy'
                log.info('Detected Spectroscopy Data from BLUE Camera')
            elif ['<NO GRATING>'] in gratings and int(np.mean(cam_targ)) == 0:
                self.technique = 'Imaging'
                log.info('Detected Imaging Data from BLUE Camera')
            else:
                log.error('It was not possible to determine observing '
                          'technique')
                self.technique = 'Unknown'

    def remove_conflictive_keywords(self):
        """Removes problematic keywords

        The blue camera has a set of keywords whose comments contain non-ascii
        characters, in particular the degree symbol. Those keyords are not
        needed in any stage of the data reduction therefore they are removed.
        The data will be overwritten with the keywords removed. The user will
        need to have backups of raw data.

        """
        log.info('Removing conflictive keywords in Blue Camera Headers')
        log.warning('Files will be overwritten')
        for blue_file in self.image_collection.file.tolist():
            full_path = os.path.join(self.args.raw_path, blue_file)
            log.debug('Processing file {:s}'.format(blue_file))
            try:
                data, header = fits.getdata(full_path,
                                            header=True,
                                            ignore_missing_end=True)

                keys_to_remove = ['PARAM0',
                                  'PARAM61',
                                  'PARAM62',
                                  'PARAM63',
                                  'NAXIS3']

                if data.ndim == 3:
                    header['NAXIS'] = 2
                    data = data[0]

                    log.debug('Modified file to be 2D instead of 3D '
                              '(problematic)')

                for keyword in keys_to_remove:
                    header.remove(keyword)

                    log.debug('Removed conflictive keyword '
                              '{:s}'.format(keyword))

                log.debug('Updated headers')

                fits.writeto(full_path,
                             data,
                             header,
                             clobber=True)

            except KeyError as error:
                log.debug(error)


if __name__ == '__main__':
    pass


