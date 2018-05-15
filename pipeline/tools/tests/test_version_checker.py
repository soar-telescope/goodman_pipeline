
__author__ = 'Bruno Quint'

import unittest

from ...tools import version


class TestVersionChecker(unittest.TestCase):

    def test_get_last(self):

        api, feature, bug = version.get_last()

        self.assertEqual(api, version.api)
        self.assertEqual(feature, version.feature)
        self.assertEqual(bug, version.bug)

    def test_check_last(self):

        temp_api = version.api
        version.api = 0

        with self.assertLogs() as cm:
            version.check_last()

        self.assertEqual(3, len(cm.records))

        version.api = temp_api


if __name__ == '__main__':
    unittest.main()