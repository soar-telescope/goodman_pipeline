from __future__ import absolute_import

from unittest import TestCase, skip
from ..wcs import WCS
import numpy as np
import os
import re
import sys
from astropy.io import fits
from astropy.modeling import (models, fitting, Model)

from ccdproc import CCDData


class TestWCS(TestCase):
    
    def setUp(self):
        self.wcs = WCS()

    def test_wcs__call__(self):
        self.assertRaisesRegex(SystemExit, '1', self.wcs)
        self.assertRaises(SystemExit, self.wcs)
        
    def test_pm_none(self):
        self.assertRaises(NotImplementedError, self.wcs._none)

    def test_pm_linear_solution(self):
        self.wcs.wcs_dict['crval'] = 3514.5662540243
        self.wcs.wcs_dict['crpix'] = 1.0
        self.wcs.wcs_dict['cdelt'] = 0.653432383822811
        self.wcs._linear_solution()
        self.assertIsInstance(self.wcs.model, Model)
        self.assertEqual(self.wcs.model.__class__.__name__, 'Linear1D')

    def test_pm_log_linear(self):
        self.assertRaises(NotImplementedError, self.wcs._log_linear)

    def test_pm_chebyshev(self):
        self.wcs.wcs_dict['order'] = 3
        self.wcs.wcs_dict['pmin'] = 1616.37
        self.wcs.wcs_dict['pmax'] = 3259.98
        self.wcs.wcs_dict['fpar'] = [5115.64008185559,
                                     535.515983711607,
                                     -0.779265625182385]
        self.wcs._chebyshev()
        self.assertIsInstance(self.wcs.model, Model)
        self.assertEqual(self.wcs.model.__class__.__name__, 'Chebyshev1D')

    def test_pm_non_linear_legendre(self):
        self.wcs.wcs_dict['order'] = 3
        self.wcs.wcs_dict['pmin'] = 1616.37
        self.wcs.wcs_dict['pmax'] = 3259.98
        self.wcs.wcs_dict['fpar'] = [5115.64008185559,
                                     535.515983711607,
                                     -0.779265625182385]
        self.wcs._non_linear_legendre()
        self.assertIsInstance(self.wcs.model, Model)
        self.assertEqual(self.wcs.model.__class__.__name__, 'Legendre1D')

    def test_pm_non_linear_lspline(self):
        self.assertRaises(NotImplementedError, self.wcs._non_linear_lspline)

    def test_pm_non_linear_cspline(self):
        self.assertRaises(NotImplementedError, self.wcs._non_linear_cspline)

