from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import glob
import os
import logging
from pipeline.core import SpectroscopicMode


from ccdproc import CCDData


class FixHeaders(object):

    def __init__(self):
        self.__path = None
        self.__pattern = None
        self._log = logging.getLogger(__name__)
        self._file_list = None
        self._spectroscopic_mode = SpectroscopicMode()

    def __call__(self, path, pattern='*.fits'):

        self.path = path
        self.pattern = pattern

        self._file_list = glob.glob(os.path.join(self.path, self.pattern))

        self._fix_keywords()

        self._fix_data_shape()

        self._add_new_keywords()

    @property
    def path(self):
        return self.__path

    @path.setter
    def path(self, value):
        self.__path = value

    @property
    def pattern(self):
        return self.__pattern

    @pattern.setter
    def pattern(self, value):
        self.__pattern = value

    def _fix_keywords(self):
        """Reads and write files to correct bad keywords using ccdproc tools"""
        for input_file in self._file_list:
            self._log.info(input_file)
            # read and write file to fix non-fits compliant keys. This should
            # fix most of the issues
            ccd = CCDData.read(input_file, unit='adu')
            ccd.write(input_file, overwrite=True)

    def _fix_data_shape(self):
        """Fix data shape

        Old Goodman's Blue camera data has three dimensional data, this has to
        be corrected in order to use the Goodman Pipeline.

        To remove one dimension from the data do as follows.

        >>>> from ccdproc import CCDData
        >>>> input_file = '/path/to/some/file.fits'
        >>>> ccd = CCDData.read(input_file, unit='adu')
        >>>> ccd.data = ccd.data[0]
        >>>> ccd.write(input_file, overwrite=True)
        """
        for input_file in self._file_list:
            self._log.info(input_file)

            # fix data shape
            ccd = CCDData.read(input_file, unit='adu')
            if int(ccd.header['NAXIS']) > 2:
                self._log.info("Fixing data shape")
                ccd.data = ccd.data[0]
                ccd.write(input_file, overwrite=True)
            else:
                self._log.debug('Shape is alright')

    def _add_new_keywords(self):
        """Performs the action of setting a new keyword to the header

        Reads a file from a list defined somewhere else on the process and
        adds all the new keywords

        >>>>  ccd = CCDData.read(input_file, unit='adu')
        >>>>  for keyword in new_keywords.keys():
        >>>>      ccd.header.set(keyword, value=new_keywords[keyword], comment='')
        >>>>  ccd.write(input_file, overwrite=True)

        """

        for input_file in self._file_list:
            ccd = CCDData.read(input_file, unit='adu')
            new_keywords = self._define_keyword_values(ccd=ccd)
            for keyword in new_keywords.keys():
                ccd.header.set(keyword,
                               value=new_keywords[keyword],
                               comment='')
                print(ccd.header[keyword])
            ccd.write(input_file, overwrite=True)

    def _define_keyword_values(self, ccd):
        """Defines new keywords and values for new keywords

        Notes:
            New keywords are some of the new keywords that the Red camera
            came with that the Blue didn't have. Some of them are:
            'INSTCONF',
            'WAVMODE',
            'ROI',
            'LAMP_HGA',
            'LAMP_NE',
            'LAMP_AR',
            'LAMP_FE',
            'LAMP_CU',
            'LAMP_QUA',
            'LAMP_QPE',
            'LAMP_BUL',
            'LAMP_DOM',
            'LAMP_DPE'

        """
        keyword_dict = {}

        lamp_keywords = ['LAMP_HGA',
                         'LAMP_NE',
                         'LAMP_AR',
                         'LAMP_FE',
                         'LAMP_CU',
                         'LAMP_QUA',
                         'LAMP_QPE',
                         'LAMP_BUL',
                         'LAMP_DOM',
                         'LAMP_DPE']

        if 'INSTCONF' not in ccd.header.keys():
            if len(ccd.header['PARAM*']) > 0:
                keyword_dict['INSTCONF'] = 'Blue'
            else:
                keyword_dict['INSTCONF'] = 'Red'

        if 'WAVMODE' not in ccd.header.keys():
            keyword_dict['WAVMODE'] = self._spectroscopic_mode(
                header=ccd.header)

        if 'ROI' not in ccd.header.keys():
            keyword_dict['ROI'] = 'user-defined'

        if not all([lamp_key in ccd.header.keys() for lamp_key in lamp_keywords]):
            if ccd.header['OBSTYPE'] != 'COMP':
                for lamp_key in lamp_keywords:
                    keyword_dict[lamp_key] = 'FALSE'
            else:
                self._log.info('Trying to retrieve lamp status from '
                               'OBJECT value')
                self._log.warning('Setting all lamp keywords to FALSE')
                for lamp_key in lamp_keywords:
                    keyword_dict[lamp_key] = 'FALSE'
        return keyword_dict


if __name__ == '__main__':
    fix_headers = FixHeaders()
    fix_headers(path='/user/simon/data/soar/work/william_data/2016-03-07')