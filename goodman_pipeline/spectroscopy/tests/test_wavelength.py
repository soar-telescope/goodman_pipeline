from __future__ import absolute_import

import numpy as np
import os

from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
from astropy.io import fits
from astropy.modeling import models, Model
from ccdproc import CCDData
from unittest import TestCase, skip
from ..wavelength import (WavelengthCalibration,
                          WavelengthSolution)

from ..redspec import get_args
from ...core import add_wcs_keys


class WavelengthCalibrationTests(TestCase):

    def setUp(self):
        argument_list = ['--data-path', os.getcwd(),
                         '--proc-path', os.getcwd(),
                         '--search-pattern', 'cfzsto',
                         '--output-prefix', 'w',
                         '--extraction', 'fractional',
                         '--reference-files', 'data/ref_comp',
                         '--max-targets', '3',
                         ]
        arguments = get_args(argument_list)
        self.wc = WavelengthCalibration(args=arguments)

        self.ccd = CCDData(data=np.random.random_sample(200),
                           meta=fits.Header(),
                           unit='adu')
        self.ccd = add_wcs_keys(ccd=self.ccd)
        self.ccd.header.set('SLIT',
                            value='1.0" long slit',
                            comment="slit [arcsec]")

    def test_add_wavelength_solution(self):
        self.wc.rms_error = 0.1
        self.wc.n_points = 100
        self.wc.n_rejections = 2

        self.wc.calibration_lamp = 'non-existent.fits'

        crval1 = 3977.948
        npix = 4060
        cdelt = 0.9910068

        x_axis = np.linspace(crval1,
                             crval1 + cdelt * npix,
                             npix)

        self.ccd = self.wc.add_wavelength_solution(ccd=self.ccd,
                                                   x_axis=x_axis)

        self.assertEqual(self.ccd.header['CTYPE1'], 'LINEAR')
        self.assertEqual(self.ccd.header['CRVAL1'], crval1)
        self.assertEqual(self.ccd.header['CRPIX1'], 1)
        self.assertAlmostEqual(self.ccd.header['CDELT1'], cdelt, places=3)
        self.assertEqual(self.ccd.header['DCLOG1'],
                         'REFSPEC1 = {:s}'.format(self.wc.calibration_lamp))
        self.assertEqual(self.ccd.header['GSP_WRMS'], self.wc.rms_error)
        self.assertEqual(self.ccd.header['GSP_WPOI'], self.wc.n_points)
        self.assertEqual(self.ccd.header['GSP_WREJ'], self.wc.n_rejections)

    @skip
    def test_automatic_wavelength_solution(self):
        pass

    @skip
    def test__cross_correlation(self):
        self.wc.lamp = self.ccd.copy()
        self.wc.serial_binning = 1

        x_axis = np.arange(0, 4060, 1)

        reference = np.zeros(4060)
        gaussian = models.Gaussian1D(stddev=2)

        for i in sorted(np.random.choice(x_axis, 30)):
            gaussian.mean.value = i
            reference += gaussian(x_axis)

        offset = np.random.choice(range(1, 15), 1)[0]
        
        for slit in [1, 2, 3, 4, 5]:

            new_array = np.append(reference[offset:], np.zeros(offset))

            if slit > 3:
                box_kernel = Box1DKernel(width=slit / 0.15)
                new_array = convolve(new_array, box_kernel)

            self.assertEqual(len(reference), len(new_array))

            self.wc.lamp.header['SLIT'] = '{:d}.0" long slit'.format(slit)

            correlation_value = self.wc._cross_correlation(reference=reference,
                                                           new_array=new_array)
            self.assertEqual(correlation_value, offset)

    def test__evaluate_solution(self):

        differences = np.array([0.5] * 10)

        clipped_differences = np.ma.masked_array(differences,
                                                 mask=[0,
                                                       0,
                                                       1,
                                                       0,
                                                       0,
                                                       1,
                                                       0,
                                                       0,
                                                       1,
                                                       0])

        rms_error, n_points, n_rej = self.wc._evaluate_solution(
            clipped_differences=clipped_differences)

        self.assertEqual(rms_error, 0.5)
        self.assertEqual(n_points, 10)
        self.assertEqual(n_rej, 3)

    @skip
    def test__get_lines_in_lamp(self):
        pass

    def test__get_spectral_characteristics(self):
        self.ccd.header.set('GRATING', 'SYZY_400')
        self.ccd.header.set('GRT_ANG', 7.5)
        self.ccd.header.set('CAM_ANG', 16.1)
        self.ccd.header.set('CCDSUM', '1 1')
        self.wc.lamp = self.ccd.copy()

        spec_charact = self.wc._get_spectral_characteristics()

        self.assertIsInstance(spec_charact, dict)
        self.assertEqual(len(spec_charact), 7)

    def test_get_wsolution(self):
        self.assertIsNone(self.wc.get_wsolution())

        self.wc.wsolution = models.Chebyshev1D(degree=2)
        self.wc.wsolution.c0.value = 3977.9485
        self.wc.wsolution.c1.value = 1.00153387
        self.wc.wsolution.c2.value = -1.2891437
        self.assertIsInstance(self.wc.get_wsolution(), Model)