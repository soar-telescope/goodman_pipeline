
__author__ = 'Bruno Quint'

import unittest

from ...tools import version

__version__ = __import__('pipeline').__version__


class TestVersionChecker(unittest.TestCase):

    def test_get_last(self):
        try:
            v = version.get_last()
            self.assertRegex(v, '^(\*|\d+(\.\d+){0,2}(\.\*)?)$')
            # self.assertEqual(v, __version__)
        except ConnectionRefusedError:
            pass

    def test_am_i_updated(self):
        self.assertTrue(version.am_i_updated(__version__))
        self.assertFalse(version.am_i_updated('v0.0.0'))


if __name__ == '__main__':
    unittest.main()
