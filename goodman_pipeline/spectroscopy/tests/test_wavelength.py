from __future__ import absolute_import

import numpy as np
import os
import re

from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
from astropy.io import fits
from astropy.modeling import models, Model
from ccdproc import CCDData
from unittest import TestCase, skip
from ..wavelength import (WavelengthCalibration)

from ..redspec import get_args
from ...core import add_wcs_keys
from ...core import ReferenceData, NoMatchFound


class WavelengthCalibrationTests(TestCase):

    def setUp(self):
        self.file_list = []
        argument_list = ['--data-path', os.getcwd(),
                         '--proc-path', os.getcwd(),
                         '--search-pattern', 'cfzsto',
                         '--output-prefix', 'w',
                         '--extraction', 'fractional',
                         '--reference-files', 'data/ref_comp',
                         '--max-targets', '3',
                         ]
        arguments = get_args(argument_list)
        self.wc = WavelengthCalibration()

        self.ccd = CCDData(data=np.random.random_sample(200),
                           meta=fits.Header(),
                           unit='adu')
        self.ccd = add_wcs_keys(ccd=self.ccd)
        self.ccd.header.set('SLIT',
                            value='1.0_LONG_SLIT',
                            comment="slit [arcsec]")
        self.ccd.header.set('GSP_FNAM',
                            value='some_name.fits',
                            comment='Name of the current file')
        self.ccd.header.set('OBSTYPE',
                            value='SPECTRUM',
                            comment='Obstype')
        self.ccd.header.set('OBJECT',
                            value='An X Object',
                            comment='Some random object name')
        self.ccd.header.set('GSP_FLAT',
                            value='some_flat_file.fits',
                            comment='The name of the flat')
        self.ccd.header.set('CCDSUM',
                            value='1 1',
                            comment='Binning')
        self.ccd.header.set('WAVMODE',
                            value='400 M1',
                            comment='wavmode')

        self.lamp = self.ccd.copy()
        self.lamp.header.set('OBSTYPE',
                             value='COMP',
                             comment='Comparison lamp obstype')
        self.lamp.header.set('OBJECT', value='HgArNe')

    def tearDown(self):
        for _file in self.file_list:
            if os.path.isfile(_file):
                os.unlink(_file)

    def test__automatic_wavelength_solution_no_match(self):
        self.wc.reference_data = ReferenceData(
            reference_dir='goodman_pipeline/data/ref_comp')

        self.lamp.header.set('OBJECT', value='PbNa')
        self.wc.lamp = self.lamp
        self.assertRaises(NoMatchFound,
                          self.wc._automatic_wavelength_solution,
                          os.getcwd())

    def test__save_wavelength_calibrated(self):
        self.wc.sci_target_file = 'target_sci_file.fits'
        fname = self.wc._save_wavelength_calibrated(ccd=self.lamp,
                                                    original_filename='file_name.fits',
                                                    save_data_to=os.getcwd(),
                                                    lamp=True)
        self.file_list.append(fname)
        self.assertEqual(fname, os.path.join(os.getcwd(), 'wfile_name.fits'))

        fname = self.wc._save_wavelength_calibrated(ccd=self.lamp,
                                                    index=1,
                                                    original_filename='file_name.fits',
                                                    save_data_to=os.getcwd())
        self.file_list.append(fname)
        self.assertEqual(fname, os.path.join(os.getcwd(), 'wfile_name_ws_1.fits'))

    def test__save_science_data(self):
        wavelength_solution = models.Chebyshev1D(degree=3)
        wavelength_solution.c0.value = 4419.161693945127
        wavelength_solution.c1.value = 1.321103785944705
        wavelength_solution.c2.value = -2.9766005683232e-06
        wavelength_solution.c3.value = -4.864180906701e-10
        fname = self.wc._save_science_data(
            ccd=self.ccd,
            wavelength_solution=wavelength_solution,
            save_to=os.getcwd(),
            index=None,
            plot_results=False,
            save_plots=False,
            plots=False)
        self.file_list.append(fname)

        fname_2 = self.wc._save_science_data(
            ccd=self.ccd,
            wavelength_solution=wavelength_solution,
            save_to=os.getcwd(),
            index=1,
            plot_results=False,
            save_plots=False,
            plots=False)
        self.file_list.append(fname_2)
        expected_name = os.getcwd() + '/w' + re.sub('.fits', '', os.path.basename(self.ccd.header['GSP_FNAM'])) + "_ws_{:d}".format(1) + ".fits"

        self.assertEqual(fname, os.getcwd() + '/w' + os.path.basename(self.ccd.header['GSP_FNAM']))
        self.assertEqual(fname_2, expected_name)

    def test___call___method_wrong_ccd(self):
        self.assertRaises(AssertionError, self.wc, [], [], '', '')

    def test___call___method_wrong_comp_list(self):
        self.assertRaises(AssertionError, self.wc, self.ccd, 'comp_list', '', '')

    def test___call___method_no_comparison_lamps_json(self):
        json_output = self.wc(ccd=self.ccd,
                              comp_list=[],
                              save_data_to='',
                              reference_data='',
                              json_output=True)

        self.assertEqual(json_output['error'], 'Unable to process without reference lamps')
        self.assertEqual(json_output['warning'], '')
        self.assertEqual(json_output['wavelength_solution'], [])

    def test___call___method_no_comparison_lamps(self):
        output = self.wc(ccd=self.ccd,
                         comp_list=[],
                         save_data_to='',
                         reference_data='',
                         json_output=False)

        self.assertIsNone(output)

    def test___call___method_one_comparison_lamps(self):
        json_output = self.wc(ccd=self.ccd,
                              comp_list=[self.lamp],
                              save_data_to='',
                              reference_data='goodman_pipeline/data/ref_comp',
                              json_output=True)

        self.assertEqual(json_output['error'], 'Unable to obtain wavelength solution')
        self.assertEqual(json_output['warning'], '')
        self.assertEqual(json_output['wavelength_solution'], [])



