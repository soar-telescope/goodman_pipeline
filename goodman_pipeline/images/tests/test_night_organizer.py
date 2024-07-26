from __future__ import absolute_import

import numpy as np
import os

import shutil

from astropy.io import fits
from ccdproc import CCDData
from unittest import TestCase, skip

from ..night_organizer import NightOrganizer


def create_fake_data(technique, instrument, path):
    if os.path.isdir(path):
        card_values = [
            {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:23:24.285',
             'obsdec': '-39:12:53.954'},
            {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:24:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:24:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:24:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:25:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:25:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:25:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:26:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:27:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:28:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:23:34.285',
             'obsdec': ' 39:13:53.954'},
            {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:23:34.285',
             'obsdec': ' 39:13:53.954'},
            {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:23:34.285',
             'obsdec': ' 39:13:53.954'},
            {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:23:34.285',
             'obsdec': ' 39:13:53.954'},
            {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:23:34.285',
             'obsdec': ' 39:13:53.954'},
            {'obstype': 'COMP', 'object': 'HgArNe', 'obsra': '16:22:34.285',
             'obsdec': '-39:23:53.954'},
            {'obstype': 'COMP', 'object': 'HgArNe', 'obsra': '16:24:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'COMP', 'object': 'HgArNe', 'obsra': '16:24:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'COMP', 'object': 'HgArNe', 'obsra': '16:24:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'COMP', 'object': 'HgArNe', 'obsra': '16:25:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'COMP', 'object': 'HgArNe', 'obsra': '16:25:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'COMP', 'object': 'HgArNe', 'obsra': '16:25:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'COMP', 'object': 'HgArNe', 'obsra': '16:26:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'COMP', 'object': 'HgArNe', 'obsra': '16:27:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'COMP', 'object': 'HgArNe', 'obsra': '16:28:34.285',
             'obsdec': '-39:13:53.954'},
            {'obstype': 'COMP', 'object': 'HgArNe', 'obsra': '16:23:34.285',
             'obsdec': '39:13:53.954'},
            {'obstype': 'COMP', 'object': 'HgArNe', 'obsra': '16:23:34.285',
             'obsdec': '39:13:53.954'},
            {'obstype': 'COMP', 'object': 'HgArNe', 'obsra': '16:23:34.285',
             'obsdec': '39:13:53.954'},
            {'obstype': 'COMP', 'object': 'HgArNe', 'obsra': '16:23:34.285',
             'obsdec': '39:13:53.954'},
            {'obstype': 'COMP', 'object': 'HgArNe', 'obsra': '16:23:34.285',
             'obsdec': '39:13:53.954'}]
        for i in range(len(card_values)):
            ccd = CCDData(data=np.ones((3, 3)),
                          meta=fits.Header(),
                          unit='adu')

            ccd.header.set('DATE', value='2019-03-22', comment='nc')
            ccd.header.set('SLIT', value='1.0" long slit', comment='nc')
            ccd.header.set('DATE-OBS', value='2019-03-22T09:59:33.654',
                           comment='nc')
            ccd.header.set('OBSTYPE', value=card_values[i]['obstype'],
                           comment='nc')
            ccd.header.set('OBJECT', value=card_values[i]['object'],
                           comment='nc')
            ccd.header.set('EXPTIME', value='10', comment='nc')
            ccd.header.set('OBSRA', value=card_values[i]['obsra'], comment='nc')
            ccd.header.set('OBSDEC', value=card_values[i]['obsdec'],
                           comment='nc')
            ccd.header.set('GRATING', value='SYZY_400', comment='nc')
            ccd.header.set('CAM_TARG', value='16.1', comment='nc')
            ccd.header.set('GRT_TARG', value='7.5', comment='nc')
            ccd.header.set('FILTER', value='<NO FILTER>', comment='nc')
            ccd.header.set('FILTER2', value='GG455', comment='nc')
            ccd.header.set('GAIN', value=1.48, comment='nc')
            ccd.header.set('RDNOISE', value=3.89, comment='nc')
            ccd.header.set('ROI', value='Spectroscopic 2x2', comment='nc')
            ccd.header.set('INSTCONF', value=instrument, comment='nc')
            ccd.header.set('WAVMODE', value=technique, comment='nc')

            ccd.write(
                os.path.join(path, 'test_file_{:03d}.fits'.format(i)))


class NightOrganizerTest(TestCase):

    def setUp(self):
        self.full_path = os.path.join(
            os.path.dirname(__file__),
            '../../data/test_data/night-organizer-test')

        if not os.path.isdir(self.full_path):
            os.mkdir(self.full_path)

        self.instrument = 'Red'
        self.technique = 'Spectroscopy'
        self.ignore_bias = False
        self.ignore_flat = False

    def tearDown(self):
        if os.path.isdir(self.full_path):
            shutil.rmtree(self.full_path)

    def test_night_organizer_red_spectroscopy_night(self):
        create_fake_data(technique=self.technique,
                         instrument=self.instrument,
                         path=self.full_path)

        night_organizer = NightOrganizer(
            full_path=self.full_path,
            instrument=self.instrument,
            technique=self.technique,
            ignore_bias=self.ignore_bias,
            ignore_flats=self.ignore_flat)

        result = night_organizer()
        self.assertIsInstance(result, list)

    def test_night_organizer_imaging_night(self):
        pass
