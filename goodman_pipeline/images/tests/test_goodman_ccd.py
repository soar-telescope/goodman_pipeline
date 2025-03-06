from __future__ import absolute_import
from unittest import TestCase, skip

from ..goodman_ccd import get_args, ReduceCCD


class MainAppTest(TestCase):

    def setUp(self):
        self.main_app = ReduceCCD()

    def test___call__(self):
        self.assertRaises(SystemExit, self.main_app)

    def test___call___show_version(self):
        arguments = ['--version']
        args = get_args(arguments=arguments)
        self.assertRaises(SystemExit, self.main_app, args)
