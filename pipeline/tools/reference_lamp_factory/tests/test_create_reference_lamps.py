from __future__ import absolute_import

from unittest import TestCase, skip
import logging
import os

from ..create_reference_lamps import ReferenceLibraryFactory

logging.disable(logging.CRITICAL)


class TestReferenceLibraryFactory(TestCase):

    def setUp(self):
        self.cwd = os.getcwd()
        arguments = ['--from', self.cwd,
                     '--save-to', self.cwd]
        self.factory = ReferenceLibraryFactory(arguments=arguments)

    def tearDown(self):
        pass

    def test_instantiation(self):
        self.assertIsInstance(self.factory, ReferenceLibraryFactory)
        self.assertIsNotNone(self.factory.args)

    def test_get_args(self):
        self.assertEqual(self.factory.args.path, self.cwd)
        self.assertEqual(self.factory.args.save_to, self.cwd)


class TestExceptions(TestCase):

    # def test_get_args(self):
    #     # rf = ReferenceLibraryFactory(arguments=None)
    #     self.assertRaisesRegex(SystemExit,
    #                            '2',
    #                            ReferenceLibraryFactory,
    #                            None)
    #     self.assertRaises(SystemExit, ReferenceLibraryFactory, None)

    def test_get_args_path_no_exist(self):
        arguments = ['--from', os.path.join(os.getcwd(), 'fake-dir'),
                     '--save-to', os.getcwd()]
        self.assertRaisesRegex(SystemExit,
                               'Input folder',
                               ReferenceLibraryFactory,
                               arguments)
        self.assertRaises(SystemExit, ReferenceLibraryFactory, arguments)

    def test_get_args_save_to_no_exist(self):
        arguments = ['--from', os.getcwd(),
                     '--save-to', os.path.join(os.getcwd(), 'fake-dir')]

        self.assertRaisesRegex(SystemExit,
                               'Output folder',
                               ReferenceLibraryFactory,
                               arguments)
        self.assertRaises(SystemExit, ReferenceLibraryFactory, arguments)
