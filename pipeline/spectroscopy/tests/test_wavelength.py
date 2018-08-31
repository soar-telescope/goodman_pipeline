from __future__ import absolute_import

import numpy as np
import os


from astropy.io import fits
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

        self.ccd = CCDData(data=np.ones((100, 100)),
                           meta=fits.Header(),
                           unit='adu')
        self.ccd = add_wcs_keys(ccd=self.ccd)

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

    def test_automatic_wavelength_solution(self):
        pass

    def test_bin_reference_data(self):
        wavelength = np.linspace(3000, 7000, 4000)
        intensity = np.random.random_sample(4000)

        for i in range(1, 4):
            self.wc.serial_binning = i

            new_wavelength, new_intensity = self.wc._bin_reference_data(
                wavelength=wavelength,
                intensity=intensity)

            self.assertEqual(len(wavelength), len(intensity))
            self.assertEqual(len(new_wavelength), len(new_intensity))
            self.assertEqual(len(new_wavelength), np.floor(len(wavelength) / i))


def test_process_spectroscopy_data():
    pass


def test_wavelength_calibration___init__():
    pass


def test_wavelength_calibration___call__():
    pass


def test_wavelength_calibration_get_wsolution():
    pass


def test_wavelength_calibration_get_calibration_lamp():
    pass


def test_wavelength_calibration_get_lines_in_lamp():
    pass


def test_wavelength_calibration_get_best_filling_value():
    pass


def test_wavelength_calibration_recenter_lines():
    pass


def test_wavelength_calibration_recenter_broad_lines():
    pass


def test_wavelength_calibration_get_spectral_characteristics():
    pass


def test_wavelength_calibration_interpolate():
    pass


def test_wavelength_calibration_recenter_line_by_data():
    pass


def test_wavelength_calibration_predicted_wavelength():
    pass


def test_wavelength_calibration_automatic_wavelength_solution():
    pass


def test_wavelength_calibration_interactive_wavelength_solution():
    pass


def test_wavelength_calibration_cross_correlation():
    pass


def test_wavelength_calibration_on_click():
    pass


def test_wavelength_calibration_key_pressed():
    pass


def test_wavelength_calibration_register_mark():
    pass


def test_wavelength_calibration_find_more_lines():
    pass


def test_wavelength_calibration_update_marks_plot():
    pass


def test_wavelength_calibration_plot_raw_over_reference():
    pass


def test_wavelength_calibration_evaluate_solution():
    pass


def test_wavelength_calibration_fit_pixel_to_wavelength():
    pass


def test_wavelength_calibration_linearize_spectrum():
    pass


def test_wavelength_calibration_add_wavelength_solution():
    pass


def test_wavelength_calibration_display_onscreen_message():
    pass


def test_wavelength_calibration_display_help_text():
    pass


def test_wavelength_solution___init__():
    pass


def test_wavelength_solution_set_spectral_features():
    pass


def test_wavelength_solution_check_compatibility():
    pass


def test_wavelength_solution_set_solution_name():
    pass

# old tests


def test_wavelength_calibration_instantiate():
    from ..wavelength import WavelengthCalibration
    from .test_redspec import test_get_args
    import os

    test_args = test_get_args()

    wavelength_calibration = WavelengthCalibration(args=test_args)

    assert isinstance(wavelength_calibration, WavelengthCalibration)
    assert os.path.exists(test_args.reference_dir)


def test_get_wsolution_none():
    pass


def test_lambda_pipeline_to_iraf():
    pass


def test_lamp_correlation():
    """Cross correlate lamp

    :return:
    """
    pass


def test_rms_pipeline_to_iraf():
    pass
