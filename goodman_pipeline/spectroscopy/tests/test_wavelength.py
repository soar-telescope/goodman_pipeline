from __future__ import absolute_import

import numpy as np
import os

from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
from astropy.io import fits
from astropy.modeling import models, Model
from ccdproc import CCDData
from unittest import TestCase, skip
from ..wavelength import (WavelengthCalibration)

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
                            value='1.0_LONG_SLIT',
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