from __future__ import absolute_import

from unittest import TestCase, skip
from ..wcs import WCS
import numpy as np
import os
import re
import sys
from astropy.io import fits
from astropy.modeling import (models, fitting, Model)
import matplotlib.pyplot as plt

from ccdproc import CCDData


class TestWCSBase(TestCase):

    def setUp(self):
        self.data_path = os.path.join(
            os.path.dirname(sys.modules['goodman_pipeline'].__file__),
            'data/test_data/wcs_data')

        self.wcs = WCS()

    @staticmethod
    def _recover_lines(ccd):
        lines_pixel = []
        lines_angstrom = []
        pixel_keywords = ccd.header['GSP_P*']
        for pixel_key in pixel_keywords:
            if re.match(r'GSP_P\d{3}', pixel_key) is not None:
                angstrom_key = re.sub('GSP_P', 'GSP_A', pixel_key)
                if int(ccd.header[angstrom_key]) != 0:
                    lines_pixel.append(float(ccd.header[pixel_key]))
                    lines_angstrom.append(float(ccd.header[angstrom_key]))
        return lines_pixel, lines_angstrom


class TestWCS(TestWCSBase):

    # def test_wcs__call__(self):
    #     self.assertRaisesRegex(SystemExit, '1', self.wcs)
    #     self.assertRaises(SystemExit, self.wcs)

    def test_fit_chebyshev(self):
        test_file = os.path.join(self.data_path,
                                 'goodman_comp_400M1_HgArNe.fits')
        ccd = CCDData.read(test_file, unit='adu')
        pixel, angstrom = self._recover_lines(ccd=ccd)
        model = self.wcs.fit(physical=pixel, wavelength=angstrom)
        self.assertIsInstance(model, Model)

        self.assertEqual(model.__class__.__name__, ccd.header['GSP_FUNC'])
        self.assertEqual(model.degree, ccd.header['GSP_ORDR'])
        for i in range(model.degree + 1):
            self.assertAlmostEqual(model.__getattribute__('c{:d}'.format(i)).value,
                             ccd.header['GSP_C{:03d}'.format(i)])

    def test_fit_linear(self):
        test_file = os.path.join(self.data_path,
                                 'goodman_comp_400M1_HgArNe.fits')
        ccd = CCDData.read(test_file, unit='adu')
        pixel, angstrom = self._recover_lines(ccd=ccd)
        model = self.wcs.fit(physical=pixel,
                             wavelength=angstrom,
                             model_name='linear')
        self.assertIsInstance(model, Model)

    def test_fit_invalid(self):
        test_file = os.path.join(self.data_path,
                                 'goodman_comp_400M1_HgArNe.fits')
        ccd = CCDData.read(test_file, unit='adu')
        pixel, angstrom = self._recover_lines(ccd=ccd)

        self.assertRaisesRegex(NotImplementedError,
                               'The model invalid is not implemented',
                               self.wcs.fit,
                               pixel,
                               angstrom,
                               'invalid')

        self.assertRaises(NotImplementedError,
                          self.wcs.fit,
                          pixel,
                          angstrom,
                          'invalid')

    def test_fit__unable_to_fit(self):
        pixel = [0, 1, 2, 3]
        angstrom = [20, 30, 40]
        # self.assertRaisesRegex(ValueError,
        #                        'x and y should have the same shape',
        #                        self.wcs.fit, pixel, angstrom)
        self.assertRaises(ValueError, self.wcs.fit, pixel, angstrom)

    def test_read__linear(self):
        test_file = os.path.join(self.data_path,
                                 'linear_fits_solution.fits')
        self.assertTrue(os.path.isfile(test_file))

        ccd = CCDData.read(test_file, unit='adu')

        result = self.wcs.read(ccd=ccd)

        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)
        self.assertIsInstance(self.wcs.get_model(), Model)

    def test_read__log_linear(self):
        test_file = os.path.join(self.data_path,
                                 'log-linear_fits_solution.fits')
        self.assertTrue(os.path.isfile(test_file))

        ccd = CCDData.read(test_file, unit='adu')
        #
        # result = self.wcs.read(ccd=ccd)
        #
        # self.assertIsInstance(result, list)
        # self.assertEqual(len(result), 2)
        # self.assertIsInstance(self.wcs.get_model(), Model)
        self.assertRaises(NotImplementedError, self.wcs.read, ccd)

    def test_read__non_linear_chebyshev(self):
        test_file = os.path.join(self.data_path,
                                 'non-linear_fits_solution_cheb.fits')
        self.assertTrue(os.path.isfile(test_file))

        ccd = CCDData.read(test_file, unit='adu')

        result = self.wcs.read(ccd=ccd)
        self.assertIsInstance(self.wcs.model, Model)
        self.assertEqual(self.wcs.model.__class__.__name__, 'Chebyshev1D')

    def test_read__non_linear_legendre(self):
        test_file = os.path.join(self.data_path,
                                 'non-linear_fits_solution_legendre.fits')
        self.assertTrue(os.path.isfile(test_file))

        ccd = CCDData.read(test_file, unit='adu')

        result = self.wcs.read(ccd=ccd)
        self.assertIsInstance(self.wcs.model, Model)
        self.assertEqual(self.wcs.model.__class__.__name__, 'Legendre1D')

    def test_read__non_linear_lspline(self):
        test_file = os.path.join(self.data_path,
                                 'non-linear_fits_solution_linear-spline.fits')
        self.assertTrue(os.path.isfile(test_file))

        ccd = CCDData.read(test_file, unit='adu')
        # self.wcs.read(ccd=ccd)
        self.assertRaises(NotImplementedError, self.wcs.read, ccd)
        self.assertRaisesRegex(NotImplementedError,
                               'Linear spline is not implemented',
                               self.wcs.read, ccd)

    def test_read__non_linear_cspline(self):
        test_file = os.path.join(self.data_path,
                                 'non-linear_fits_solution_cubic-spline.fits')
        self.assertTrue(os.path.isfile(test_file))

        ccd = CCDData.read(test_file, unit='adu')
        self.assertRaises(NotImplementedError, self.wcs.read, ccd)
        self.assertRaisesRegex(NotImplementedError,
                               'Cubic spline is not implemented',
                               self.wcs.read, ccd)

    def test_write_fits_wcs(self):
        self.assertRaises(NotImplementedError, self.wcs.write_fits_wcs,
                          None,
                          None)

    def test_read__invalid(self):
        test_file = os.path.join(self.data_path,
                                 'linear_fits_solution.fits')
        self.assertTrue(os.path.isfile(test_file))

        ccd = CCDData.read(test_file, unit='adu')
        ccd.wcs.wcs.ctype[0] = 'INVALID'

        self.assertRaisesRegex(NotImplementedError,
                               'CTYPE INVALID is not recognized',
                               self.wcs.read,
                               ccd)
        self.assertRaises(NotImplementedError, self.wcs.read, ccd)

    def test_write_gsp_wcs(self):
        test_file = os.path.join(self.data_path,
                                 'goodman_comp_400M1_HgArNe.fits')
        ccd = CCDData.read(test_file, unit='adu')
        pixel, angstrom = self._recover_lines(ccd=ccd)
        model = self.wcs.fit(physical=pixel, wavelength=angstrom)
        self.assertIsInstance(model, Model)

        blank_ccd = CCDData(data=np.ones(ccd.data.shape),
                          meta=fits.Header(),
                          unit='adu')
        blank_ccd.header.set('GSP_WREJ', value=None, comment='empty')

        new_ccd = self.wcs.write_gsp_wcs(ccd=blank_ccd, model=model)

        self.assertEqual(new_ccd.header['GSP_FUNC'], ccd.header['GSP_FUNC'])
        self.assertEqual(new_ccd.header['GSP_ORDR'], ccd.header['GSP_ORDR'])
        self.assertEqual(new_ccd.header['GSP_NPIX'], ccd.header['GSP_NPIX'])
        for i in range(model.degree + 1):
            self.assertAlmostEqual(new_ccd.header['GSP_C{:03d}'.format(i)],
                             ccd.header['GSP_C{:03d}'.format(i)])

    def test_read_gsp_wcs(self):
        test_file = os.path.join(self.data_path,
                                 'goodman_comp_400M1_HgArNe.fits')

        self.assertTrue(os.path.isfile(test_file))
        ccd = CCDData.read(test_file, unit='adu')
        result = self.wcs.read_gsp_wcs(ccd=ccd)

        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)
        self.assertIsInstance(self.wcs.get_model(), Model)

    def test_get_model_is_None(self):
        self.wcs.model = None
        self.assertIsNone(self.wcs.get_model())

    def test_get_model_is_not_None(self):
        self.wcs.model = models.Chebyshev1D(degree=3)
        self.assertIsInstance(self.wcs.get_model(), Model)

    def test_pm_none(self):
        # test_file = os.path.join(self.data_path,
        #                          'non-linear_fits_solution_cheb.fits')
        # self.assertTrue(os.path.isfile(test_file))
        #
        # ccd = CCDData.read(test_file, unit='adu')
        #
        # WAT2_001 = 'wtype = multispec spec1 = "1 1 2 1. 1.5114461210693 4096 0. 834.39 864'
        # WAT2_002 = '.39 1. 0. 1 3 1616.37 3259.98 5115.64008185559 535.515983711607 -0.7'
        # WAT2_003 = '79265625182385"'
        #
        # dtype = -1
        self.assertRaises(NotImplementedError, self.wcs._none)
