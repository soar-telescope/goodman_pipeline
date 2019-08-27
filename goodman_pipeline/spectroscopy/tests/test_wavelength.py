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
        self.wc = WavelengthCalibration()

        self.ccd = CCDData(data=np.random.random_sample(200),
                           meta=fits.Header(),
                           unit='adu')
        self.ccd = add_wcs_keys(ccd=self.ccd)
        self.ccd.header.set('SLIT',
                            value='1.0_LONG_SLIT',
                            comment="slit [arcsec]")



    @skip
    def test_automatic_wavelength_solution(self):
        pass