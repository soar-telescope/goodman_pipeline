from __future__ import absolute_import

import argparse
import os
from unittest import TestCase

from ...spectroscopy.redspec import (get_args, MainApp)


class TestArguments(TestCase):

    def setUp(self):
        self.dummy_path = 'dummy_path'
        os.mkdir(self.dummy_path)

    def tearDown(self):
        if os.path.exists(self.dummy_path):
            os.rmdir(self.dummy_path)

    def test_input_path_is_relative(self):
        arguments = ['--data-path', self.dummy_path]
        args = get_args(arguments)
        self.assertTrue(os.path.isabs(args.source))

    def test_input_path_is_absolute(self):
        absolute_path = os.path.join(os.getcwd(), self.dummy_path)
        arguments = ['--data-path', absolute_path]
        args = get_args(arguments)
        self.assertTrue(os.path.isabs(args.source))

    def test_convert_from_relative_to_absolute(self):

        self.assertFalse(os.path.isabs(self.dummy_path))

        arguments = ['--data-path', self.dummy_path]
        args = get_args(arguments)
        self.assertTrue(os.path.isabs(args.source))

    def test_absolute_path_does_not_exists(self):
        if os.path.exists(self.dummy_path):
            os.rmdir(self.dummy_path)
        arguments = ['--data-path', self.dummy_path]
        self.assertRaises(SystemExit, get_args, arguments)

    def test_absolute_path_exists(self):
        arguments = ['--data-path', self.dummy_path]
        args = get_args(arguments)

        self.assertTrue(os.path.exists(args.source))
        self.assertTrue(os.path.isdir(args.source))

        # If yes, carry on

            # If not, convert it to become absolute

        # Check if the source folder exists

        #  If it exists, carry on

        # If source folder does not exists, print message and leave program
        # error_message = 'Input folder "{} does not exists."'


def test_get_args():
    from ...spectroscopy.redspec import get_args
    import argparse
    arguments = ['--data-path', './',
                 '--proc-path', './',
                 '--search-pattern', 'test-pattern',
                 '--output-prefix', 'w',
                 '--extraction', 'fractional']

    args = get_args(arguments)

    assert isinstance(args, argparse.Namespace)
    assert args.pattern == 'test-pattern'
    return args


class TestMainApp(TestCase):

    def setUp(self):
        self.main_app = MainApp()

    def test_instantiation_without_args(self):
        self.assertIsInstance(self.main_app, MainApp)
        self.assertIsNone(self.main_app.args)
        self.assertIsNone(self.main_app.wavelength_solution_obj)
        self.assertIsNone(self.main_app.wavelength_calibration)
        self.assertIsNone(self.main_app.reference)

    def test_main_app(self):
        pass


if __name__ == '__main__':
    test_get_args()
