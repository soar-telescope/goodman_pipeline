
__author__ = 'Bruno Quint'

import os
import unittest

from goodman_pipeline.core import check_version

__version__ = __import__('goodman_pipeline').__version__


class TestVersionChecker(unittest.TestCase):

    def test_get_last(self):
        try:
            v = check_version.get_last()
            self.assertRegex(v, '^(\*|\d+(\.\d+){0,2}(\.\*)?)$')
            # self.assertEqual(v, __version__)
        except ConnectionRefusedError:  # pragma: no cover
            pass

    def test_get_last_no_token(self):
        print(os.environ['GITHUB_ACCESS_TOKEN'])
        try:
            del os.environ['GITHUB_ACCESS_TOKEN']
            v = check_version.get_last()
            self.assertRegex(v, '^(\*|\d+(\.\d+){0,2}(\.\*)?)$')
            # self.assertEqual(v, __version__)
        except ConnectionRefusedError:  # pragma: no cover
            pass

    def test_am_i_updated(self):
        try:
            self.assertTrue(check_version.am_i_updated(__version__))
            self.assertFalse(check_version.am_i_updated('v0.0.0'))
        except ConnectionRefusedError:  # pragma: no cover
            pass
