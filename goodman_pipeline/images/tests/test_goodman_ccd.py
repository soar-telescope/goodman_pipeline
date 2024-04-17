from __future__ import absolute_import
from unittest import TestCase, skip

from ..goodman_ccd import get_args, RedCCD


class MainAppTest(TestCase):

    def setUp(self):
        self.red_ccd = RedCCD()

    def test___call__(self):
        self.assertRaises(SystemExit, self.red_ccd)

    def test___call___show_version(self):
        arguments = ['--version']
        args = get_args(arguments=arguments)
        self.assertRaises(SystemExit, self.red_ccd, args)
