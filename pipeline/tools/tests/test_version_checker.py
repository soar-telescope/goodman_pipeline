
__author__ = 'Bruno Quint'

import unittest

from ... import info
from ...tools import version


class TestVersionChecker(unittest.TestCase):

    def test_get_last(self):
         v = version.get_last()
         self.assertEqual(v, info.__version__)

    def test_am_i_updated(self):
        self.assertTrue(version.am_i_updated(info.__version__))
        self.assertFalse(version.am_i_updated('v0.0.0'))


if __name__ == '__main__':
    unittest.main()