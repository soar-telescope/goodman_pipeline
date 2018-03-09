from __future__ import absolute_import

from ..wavelength import (WavelengthCalibration,
                          WavelengthSolution)

# TODO (simon): Figure out what tests are required and/or if some major
# refactoring is required.


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
