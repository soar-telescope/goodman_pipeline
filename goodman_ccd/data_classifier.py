from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import logging
import numpy as np
import random
from ccdproc import ImageFileCollection
from .core import fix_duplicated_keywords, remove_conflictive_keywords

log = logging.getLogger('goodmanccd.dataclassifier')

# TODO (simon): Move methods as functions to goodman_ccd.core.py


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
        no_bias_collection = self.image_collection[
            self.image_collection.obstype != 'BIAS']
        self.get_obs_technique(image_collection=no_bias_collection)
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

        Goodman has two camera, *Blue* and *Red*. They are, as the name suggest
        optimized for bluer and redder wavelength respectively. Their headers
        are different so this methods uses their differences to discover which
        camera the data belong to. The red camera has an specific keyword that
        says which camera is but the blue does not.
        The result is stored as an attribute of the class.

        Notes:
            As of April 2017 the blue camera computer was upgraded and as a
            result the headers where updated too. But we need to keep this
            feature for *backward compatibility*

        Args:
            night_folder (str): The full path for the raw data location

        """
        while True:
            try:
                ifc = ImageFileCollection(night_folder)
                self.image_collection = ifc.summary.to_pandas()

                self.objects_collection = self.image_collection[
                    self.image_collection.obstype != 'BIAS']

                if len(self.objects_collection) > 0:

                    indexes = self.objects_collection.index.tolist()
                    index = random.choice(indexes)

                    try:

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

    def get_obs_technique(self, image_collection):
        """Identify if the data is Imaging or Spectroscopy

        Besides the fact there are two cameras there are two observational
        techniques. Imaging and Spectroscopy. The red camera has an specific
        keyword that contains that information but the blue does not.

        The result is stored as an attribute of the class.

        """

        if self.instrument == 'Red':
            wavmodes = image_collection.wavmode.unique()
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
            file_list = self.image_collection.file.tolist()
            remove_conflictive_keywords(path=self.args.raw_path,
                                        file_list=file_list)
            # gratings = image_collection.grating.unique()
            cam_targ = image_collection.cam_targ.unique()

            if int(np.mean(cam_targ)) != 0:
                self.technique = 'Spectroscopy'
                log.info('Detected Spectroscopy Data from BLUE Camera')
            elif int(np.mean(cam_targ)) == 0:
                self.technique = 'Imaging'
                log.info('Detected Imaging Data from BLUE Camera')
            else:
                log.error('It was not possible to determine observing '
                          'technique')
                self.technique = 'Unknown'


if __name__ == '__main__':
    pass


