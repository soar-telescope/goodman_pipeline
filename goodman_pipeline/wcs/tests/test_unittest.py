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

    def test_binning_property_getter(self):
        self.assertEqual(self.wcs.binning, 1)

    def test_binning_property_setter_linear_model(self):
        self.wcs.wcs_dict['crval'] = 3514.5662540243
        self.wcs.wcs_dict['crpix'] = 1.0
        self.wcs.wcs_dict['cdelt'] = 0.653432383822811
        self.wcs._linear_solution()

        self.wcs.binning = 2

        self.assertEqual(self.wcs.binning, 2)

    def test_binning_property_setter_non_linear_model(self):
        self.wcs.wcs_dict['order'] = 3
        self.wcs.wcs_dict['pmin'] = 1616.37
        self.wcs.wcs_dict['pmax'] = 3259.98
        self.wcs.wcs_dict['fpar'] = [5115.64008185559,
                                     535.515983711607,
                                     -0.779265625182385]
        self.wcs._chebyshev()

        self.wcs.binning = 2

        self.assertEqual(self.wcs.binning, 2)

    def test_binning_property_setter_undefined_model(self):
        with self.assertRaises(NotImplementedError):
            self.wcs.binning = 2

    def test_binning_property_setter_invalid_input(self):
        with self.assertRaises(NotImplementedError):
            self.wcs.binning = 0

    def test_pm_fitter_undefined_model_and_fitter(self):
        pixel = list(range(100))
        angstrom = list(range(100))
        self.assertRaises(RuntimeError, self.wcs._fitter, pixel, angstrom)
        self.assertRaisesRegex(RuntimeError,
                               "Undefined model and fitter",
                               self.wcs._fitter,
                               pixel, angstrom)
        # self.wcs._fitter(physical=pixel, wavelength=angstrom)

    def test_pm_fitter_not_enough_points(self):
        pixel = [1,2]
        angstrom = [8000, 8001]

        self.wcs.model_name = 'chebyshev'
        self.wcs.degree = 3
        self.wcs._model_constructor()

        result = self.wcs._fitter(physical=pixel, wavelength=angstrom)
        self.assertIsNone(result)


    def test_pm_set_math_model__none(self):
        self.wcs.wcs_dict['dtype'] = -1
        self.assertRaises(NotImplementedError, self.wcs._set_math_model)

    def test_pm_set_math_model__log_linear(self):
        self.wcs.wcs_dict['dtype'] = 1
        # self.wcs.wcs_dict['crval'] = 1
        # self.wcs.wcs_dict['cdelt'] = 1
        # self.wcs.wcs_dict['crpix'] = 1
        self.assertRaises(NotImplementedError, self.wcs._set_math_model)

    def test_pm_set_math_model__pixel_coordinates(self):
        self.wcs.wcs_dict['dtype'] = 2
        self.wcs.wcs_dict['ftype'] = 5
        self.assertRaises(NotImplementedError, self.wcs._set_math_model)

    def test_pm_set_math_model__sampled_coordinate_array(self):
        self.wcs.wcs_dict['dtype'] = 2
        self.wcs.wcs_dict['ftype'] = 6
        self.assertRaises(NotImplementedError, self.wcs._set_math_model)

    def test_pm_set_math_model__wrong_dtype(self):
        self.wcs.wcs_dict['dtype'] = 3
        # self.wcs.wcs_dict['ftype'] = 6
        self.assertRaises(SyntaxError, self.wcs._set_math_model)
        self.assertRaisesRegex(SyntaxError,
                               'dtype {:d} is not defined in the '
                               'standard'.format(self.wcs.wcs_dict['dtype']))

    def test_pm_set_math_model__wrong_ftype(self):
        self.wcs.wcs_dict['dtype'] = 2
        self.wcs.wcs_dict['ftype'] = 7
        self.assertRaises(SyntaxError, self.wcs._set_math_model)
        self.assertRaisesRegex(SyntaxError,
                               'ftype {:d} is not defined in the '
                               'standard'.format(self.wcs.wcs_dict['ftype']))
        
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

