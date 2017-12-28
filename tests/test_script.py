from __future__ import absolute_import


def test_class_instantiate():
    from ..pipeline.spectroscopy import MainApp
    main_app = MainApp()
    assert isinstance(main_app, MainApp)
    assert main_app.args is None
    assert main_app.wavelength_solution_obj is None
