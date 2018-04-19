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
        self.data_path = os.path.join(
            os.path.dirname(sys.modules['pipeline'].__file__),
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

    def test_wcs__call__(self):
        self.assertRaisesRegex(SystemExit, '1', self.wcs)
        self.assertRaises(SystemExit, self.wcs)

    def test_fit(self):
        test_file = os.path.join(self.data_path,
                                 'goodman_comp_400M1_HgArNe.fits')
        ccd = CCDData.read(test_file, unit='adu')
        pixel, angstrom = self._recover_lines(ccd=ccd)
        model = self.wcs.fit(physical=pixel, wavelength=angstrom)
        self.assertIsInstance(model, Model)

        self.assertEqual(model.__class__.__name__, ccd.header['GSP_FUNC'])
        self.assertEqual(model.degree, ccd.header['GSP_ORDR'])
        for i in range(model.degree + 1):
            self.assertAlmostEqual(model.__getattr__('c{:d}'.format(i)).value,
                             ccd.header['GSP_C{:03d}'.format(i)])

    def test_read__linear(self):
        test_file = os.path.join(self.data_path,
                                 'linear_fits_solution.fits')
        self.assertTrue(os.path.isfile(test_file))

        ccd = CCDData.read(test_file, unit='adu')

        result = self.wcs.read(ccd=ccd)

        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)
        self.assertIsInstance(self.wcs.get_model(), Model)

    def test_read__non_linear(self):
        test_file = os.path.join(self.data_path,
                                 'non-linear_fits_solution_cheb.fits')
        self.assertTrue(os.path.isfile(test_file))

        ccd = CCDData.read(test_file, unit='adu')

        self.assertRaises(NotImplementedError,  self.wcs.read, ccd)

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


    # @skip
    # def test_pm_model_constructor(self):
    #     self.fail()
    #
    # @skip
    # def test_pm_fitter(self):
    #     self.fail()
    #
    # @skip
    # def test_pm_read_non_linear(self):
    #     self.fail()
    #
    # @skip
    # def test_pm_read_linear(self):
    #     self.fail()
    #
    # @skip
    # def test_pm_set_math_model(self):
    #     self.fail()
    #
    # @skip
    # def test_pm_none(self):
    #     self.fail()
    #
    # @skip
    # def test_pm_linear_solution(self):
    #     self.fail()
    #
    # @skip
    # def test_pm_log_linear(self):
    #     self.fail()
    #
    # @skip
    # def test_pm_chebyshev(self):
    #     self.fail()
    #
    # @skip
    # def test_pm_non_linear_legendre(self):
    #     self.fail()
    #
    # @skip
    # def test_pm_non_linear_lspline(self):
    #     self.fail()
    #
    # @skip
    # def test_pm_non_linear_cspline(self):
    #     self.fail()

