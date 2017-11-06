from __future__ import absolute_import


def test_class_instantiate():
    from ..goodman import spectroscopy
    main_app = spectroscopy.MainApp()
    assert isinstance(main_app, spectroscopy.MainApp)