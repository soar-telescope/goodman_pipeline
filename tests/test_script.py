from __future__ import absolute_import


def test_class_instantiate():
    from goodman.spectroscopy import MainApp
    main_app = MainApp()
    assert isinstance(main_app, MainApp)