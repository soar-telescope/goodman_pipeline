
__author__ = 'Bruno Quint'

import unittest

from goodman_pipeline.core import check_version

__version__ = __import__('goodman_pipeline').__version__


class TestVersionChecker(unittest.TestCase):

    def test_get_last(self):
        try:
            v = check_version.get_last()
            self.assertRegex(v, '^(\*|\d+(\.\d+){0,2}(\.\*)?)$')
            # self.assertEqual(v, __version__)
        except ConnectionRefusedError:
            pass

    def test_am_i_updated(self):
        try:
            self.assertTrue(check_version.am_i_updated(__version__))
            self.assertFalse(check_version.am_i_updated('v0.0.0'))
        except ConnectionRefusedError:
            pass


if __name__ == '__main__':
    unittest.main()
