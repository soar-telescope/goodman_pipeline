
__author__ = 'Bruno Quint'

import unittest

from ... import info
from ...tools import version


class TestVersionChecker(unittest.TestCase):

    def test_get_last(self):

        api, feature, bug = version.get_last()

        self.assertEqual(api, info.api)
        self.assertEqual(feature, info.feature)
        self.assertEqual(bug, info.bug)

    def test_check_last_positive(self):
        with self.assertLogs() as cm:
            version.check_last(info.__version__)

    def test_check_last_negative(self):
        with self.assertLogs() as cm:
            version.check_last('v0.0.0')

if __name__ == '__main__':
    unittest.main()