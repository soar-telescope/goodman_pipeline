from __future__ import absolute_import


def test_wavelength_calibration_isntantiate():
    from ..goodman.spectroscopy.wavelength import WavelengthCalibration
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
