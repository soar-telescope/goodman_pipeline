from __future__ import absolute_import

import argparse
import os
from unittest import TestCase

from ...spectroscopy.redspec import (get_args, MainApp)


class TestMainApp(TestCase):

    def setUp(self):
        self.main_app = MainApp()
        self.dummy_path = 'dummy_path'
        os.mkdir(self.dummy_path)

    def tearDown(self):
        if os.path.exists(self.dummy_path):
            os.rmdir(self.dummy_path)

    def test_instantiation_without_args(self):
        self.assertIsInstance(self.main_app, MainApp)
        self.assertIsNone(self.main_app.args)
        self.assertIsNone(self.main_app.wavelength_solution_obj)
        self.assertIsNone(self.main_app.wavelength_calibration)
        self.assertIsNone(self.main_app.reference)

    def test_input_path_is_relative(self):
        arguments = ['--data-path', self.dummy_path]
        args = get_args(arguments=arguments)
        self.main_app.args = args
        self.main_app._check_args()
        self.assertTrue(os.path.isabs(self.main_app.args.source))

    def test_input_path_is_absolute(self):
        absolute_path = os.path.join(os.getcwd(), self.dummy_path)
        arguments = ['--data-path', absolute_path]
        args = get_args(arguments)
        self.assertTrue(os.path.isabs(args.source))

    def test_convert_from_relative_to_absolute(self):

        self.assertFalse(os.path.isabs(self.dummy_path))

        arguments = ['--data-path', self.dummy_path]
        args = get_args(arguments)
        self.main_app.args = args
        self.main_app._check_args()
        self.assertTrue(os.path.isabs(self.main_app.args.source))

    def test_absolute_path_does_not_exists(self):
        if os.path.exists(self.dummy_path):
            os.rmdir(self.dummy_path)
        arguments = ['--data-path', self.dummy_path]
        args = get_args(arguments=arguments)
        # self.main_app.args = args
        self.assertRaises(SystemExit, self.main_app, args)

    def test_absolute_path_exists(self):
        arguments = ['--data-path', self.dummy_path]
        args = get_args(arguments)

        self.assertTrue(os.path.exists(args.source))
        self.assertTrue(os.path.isdir(args.source))

    def test_main_app(self):
        pass
