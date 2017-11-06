from __future__ import absolute_import

def test_class_initialization():
    from ..spectroscopy.wavelength import WavelengthCalibration
    from .test_redspec import test_get_args
    import os

    test_args = test_get_args()

    wavelength_calibration = WavelengthCalibration(args=test_args)

    assert isinstance(wavelength_calibration, WavelengthCalibration)
    assert os.path.exists(test_args.reference_dir)


def test_get_wsolution_none():
    pass