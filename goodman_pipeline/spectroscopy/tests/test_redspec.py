from __future__ import absolute_import

import argparse
import os
import shutil
from unittest import TestCase

from ...spectroscopy.redspec import (get_args, ReduceSpectroscopy)


class TestArguments(TestCase):

    def setUp(self):
        self.dummy_path = 'dummy_path'
        os.mkdir(self.dummy_path)

    def tearDown(self):
        for _path in [self.dummy_path,
                      'goodman_pipeline/testing/reference_files',
                      'goodman_pipeline/testing',
                      'testing/processed_files',
                      'testing']:
            print(_path)
            if os.path.exists(_path):
                print("Deleting {}".format(_path))
                os.rmdir(_path)

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

    def test_get_args_output_path_does_not_exist(self):
        arguments = ['--data-path', 'does_not_exist']
        self.assertRaises(SystemExit, get_args, arguments)

    def test_reference_dir_does_not_exists(self):
        arguments = ['--reference-files', 'testing/reference_files']

        args = get_args(arguments=arguments)
        self.assertTrue(os.path.exists(
            os.path.join(os.path.dirname(__file__),
                         '../../testing/reference_files')))
        self.assertTrue(os.path.isabs(args.reference_dir))

    def test_reference_dir_os_error(self):
        arguments = ['--reference-files', '/testing/reference_files']

        args = get_args(arguments=arguments)
        print(args)
        self.assertFalse(os.path.exists('/testing/reference_files'))

    def test_destination_path_not_abs(self):
        arguments = ['--proc-path', 'testing/processed_files']
        args = get_args(arguments=arguments)
        self.assertTrue(os.path.isabs(args.destination))
        self.assertTrue(os.path.exists(args.destination))
        self.assertEqual(args.destination,
                         os.path.join(os.getcwd(),
                                      'testing/processed_files'))

    def test_destination_path_sysexit(self):
        arguments = ['--proc-path', '/testing/processed_files']
        self.assertRaises(SystemExit, get_args, arguments)


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
    assert os.path.normpath(args.source) == os.getcwd()
    assert os.path.normpath(args.destination) == os.getcwd()
    assert args.pattern == 'test-pattern'
    assert args.output_prefix == 'w'
    assert args.extraction_type == 'fractional'


class TestReduceSpectroscopy(TestCase):

    def setUp(self):
        self.main_app = ReduceSpectroscopy()

    def test_instantiation_without_args(self):
        self.assertIsInstance(self.main_app, ReduceSpectroscopy)
        self.assertIsNone(self.main_app.args)
        self.assertIsNone(self.main_app.wavelength_solution_obj)
        self.assertIsNone(self.main_app.wavelength_calibration)
        self.assertIsNone(self.main_app.reference)

    def test___call___no_args(self):
        self.assertRaises(SystemExit, self.main_app)

    def test___call___with_valid_arguments(self):
        arguments = ['--data-path', './',
                     '--proc-path', './',
                     '--search-pattern', 'test-pattern',
                     '--output-prefix', 'w',
                     '--extraction', 'fractional']
        args = get_args(arguments=arguments)
        self.assertRaises(SystemExit, self.main_app, args)
