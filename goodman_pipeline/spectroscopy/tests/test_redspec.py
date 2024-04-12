from __future__ import absolute_import

import os
from unittest import TestCase

from ...spectroscopy.redspec import (get_args, RedSpec)


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

    def test_get_args_output_path_does_not_exist(self):
        arguments = ['--data-path', 'does_not_exist']
        self.assertRaises(SystemExit, get_args, arguments)

    def test_reference_dir_does_not_exists(self):
        arguments = ['--reference-files', 'testing/reference_files']

        args = get_args(arguments=arguments)
        self.assertTrue(os.path.exists(
            os.path.join(os.getcwd(),
                         'goodman_pipeline/testing/reference_files')))
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
    assert args.pattern == 'test-pattern'
    return args


class TestRedSpec(TestCase):

    def setUp(self):
        self.red_spec = RedSpec()

    def test_instantiation_without_args(self):
        self.assertIsInstance(self.red_spec, RedSpec)
        self.assertIsNone(self.red_spec.args)
        self.assertIsNone(self.red_spec.wavelength_solution_obj)
        self.assertIsNone(self.red_spec.wavelength_calibration)
        self.assertIsNone(self.red_spec.reference)

    def test___call___no_args(self):
        self.assertRaises(SystemExit, self.red_spec)

    def test___call___with_valid_arguments(self):
        arguments = ['--data-path', './',
                     '--proc-path', './',
                     '--search-pattern', 'test-pattern',
                     '--output-prefix', 'w',
                     '--extraction', 'fractional']
        args = get_args(arguments=arguments)
        self.assertRaises(SystemExit, self.main_app, args)


if __name__ == '__main__':
    test_get_args()
