from __future__ import absolute_import

import numpy as np
import os

import shutil

from astropy.io import fits
from ccdproc import CCDData
from unittest import TestCase, skip

from ..data_classifier import DataClassifier


class DataClassifierTests(TestCase):

    def setUp(self):
        self.raw_path = os.path.join(
            os.path.dirname(__file__),
            '../../data/test_data/classify-data')

        if not os.path.isdir(self.raw_path):
            os.mkdir(self.raw_path)

        self.create_fake_spectroscopic_data()

        self.data_classifier = DataClassifier()


    def create_fake_spectroscopic_data(self):
        if os.path.isdir(self.raw_path):
            card_values = [
                {'obstype': 'OBJECT', 'object': 'NGC2070',
                 'obsra': '16:23:24.285', 'obsdec': '-39:12:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070',
                 'obsra': '16:24:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070',
                 'obsra': '16:24:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070',
                 'obsra': '16:24:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070',
                 'obsra': '16:25:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070',
                 'obsra': '16:25:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070',
                 'obsra': '16:25:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070',
                 'obsra': '16:26:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070',
                 'obsra': '16:27:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070',
                 'obsra': '16:28:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070',
                 'obsra': '16:23:34.285', 'obsdec': ' 39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070',
                 'obsra': '16:23:34.285', 'obsdec': ' 39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070',
                 'obsra': '16:23:34.285', 'obsdec': ' 39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070',
                 'obsra': '16:23:34.285', 'obsdec': ' 39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070',
                 'obsra': '16:23:34.285', 'obsdec': ' 39:13:53.954'},
                {'obstype': 'COMP', 'object': 'HgArNe',
                 'obsra': '16:22:34.285', 'obsdec': '-39:23:53.954'},
                {'obstype': 'COMP', 'object': 'HgArNe',
                 'obsra': '16:24:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP', 'object': 'HgArNe',
                 'obsra': '16:24:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP', 'object': 'HgArNe',
                 'obsra': '16:24:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP', 'object': 'HgArNe',
                 'obsra': '16:25:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP', 'object': 'HgArNe',
                 'obsra': '16:25:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP', 'object': 'HgArNe',
                 'obsra': '16:25:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP', 'object': 'HgArNe',
                 'obsra': '16:26:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP', 'object': 'HgArNe',
                 'obsra': '16:27:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP', 'object': 'HgArNe',
                 'obsra': '16:28:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP', 'object': 'HgArNe',
                 'obsra': '16:23:34.285', 'obsdec': '39:13:53.954'},
                {'obstype': 'COMP', 'object': 'HgArNe',
                 'obsra': '16:23:34.285', 'obsdec': '39:13:53.954'},
                {'obstype': 'COMP', 'object': 'HgArNe',
                 'obsra': '16:23:34.285', 'obsdec': '39:13:53.954'},
                {'obstype': 'COMP', 'object': 'HgArNe',
                 'obsra': '16:23:34.285', 'obsdec': '39:13:53.954'},
                {'obstype': 'COMP', 'object': 'HgArNe',
                 'obsra': '16:23:34.285', 'obsdec': '39:13:53.954'}]
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
                ccd.header.set('OBSRA', value=card_values[i]['obsra'],
                               comment='nc')
                ccd.header.set('OBSDEC', value=card_values[i]['obsdec'],
                               comment='nc')
                ccd.header.set('GRATING', value='SYZY_400', comment='nc')
                ccd.header.set('CAM_TARG', value='16.1', comment='nc')
                ccd.header.set('GRT_TARG', value='7.5', comment='nc')
                ccd.header.set('FILTER', value='<NO FILTER>', comment='nc')
                ccd.header.set('FILTER2', value='GG455', comment='nc')
                ccd.header.set('GAIN', value='1.48', comment='nc')
                ccd.header.set('RDNOISE', value='3.89', comment='nc')

                ccd.header.set('INSTCONF', value='Red', comment='nc')
                ccd.header.set('WAVMODE', value='Spectroscopy', comment='nc')

                ccd.write(os.path.join(self.raw_path,
                                       'test_file_{:03d}.fits'.format(i)))



    def tearDown(self):
        if os.path.isdir(self.raw_path):
            shutil.rmtree(self.raw_path)

    def test___repr__undefined(self):
        with self.assertRaises(TypeError):
            self.data_classifier.__repr__()

    def test___repr__(self):
        self.data_classifier(raw_path=self.raw_path)
        result = self.data_classifier.__repr__()
        self.assertIn('Raw Path: {:s}'.format(self.raw_path), result)
        self.assertIn('Instrument: Red Camera', result)
        self.assertIn('Observing Technique: Spectroscopy', result)

    def test_data_classifier_expected_usage(self):

        assert isinstance(self.data_classifier, DataClassifier)
        self.data_classifier(raw_path=self.raw_path)

        self.assertEqual('Red', self.data_classifier.instrument)
        self.assertEqual('Spectroscopy', self.data_classifier.technique)

    def test_data_classifier_all_imaging(self):

        for _file in os.listdir(self.raw_path):
            raw_path_full = os.path.join(self.raw_path, _file)
            recovered_ccd = CCDData.read(raw_path_full, unit='adu')
            # recovered_ccd.header['INSTCONF'] = 'Blue'
            recovered_ccd.header['WAVMODE'] = 'Imaging'

            recovered_ccd.write(raw_path_full, overwrite=True)

        self.data_classifier(raw_path=self.raw_path)
        self.assertEqual('Imaging', self.data_classifier.technique)

    def test_data_classifier_mixed_technique(self):
        sample_file = os.listdir(self.raw_path)[0]
        raw_path_full = os.path.join(self.raw_path, sample_file)
        recovered_ccd = CCDData.read(raw_path_full, unit='adu')
        # recovered_ccd.header['INSTCONF'] = 'Blue'
        recovered_ccd.header['WAVMODE'] = 'Imaging'

        recovered_ccd.write(raw_path_full, overwrite=True)

        self.data_classifier(raw_path=self.raw_path)

        self.assertEqual('Spectroscopy', self.data_classifier.technique)

    def test_data_classifier_mixed_instconf(self):
        sample_file = os.listdir(self.raw_path)[0]
        raw_path_full = os.path.join(self.raw_path, sample_file)
        recovered_ccd = CCDData.read(raw_path_full, unit='adu')
        recovered_ccd.header['INSTCONF'] = 'Blue'
        # recovered_ccd.header['WAVMODE'] = 'Imaging'
        recovered_ccd.write(raw_path_full, overwrite=True)


        with self.assertRaises(SystemExit):
            self.data_classifier(raw_path=self.raw_path)

