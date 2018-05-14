
__author__ = 'Bruno Quint'

import unittest

from ...tools import version


class TestVersionChecker(unittest.TestCase):

    def setUp(self):
        self.url = \
            "https://github.com/soar-telescope/goodman/releases/tag/1.0.0"

    def test_version_checker(self):

        class Version:
            api = 1
            feature = 0
            bug = 0

        api, feature, bug = version.check_last(self.url)

        self.assertEqual(api, Version.api)
        self.assertEqual(feature, Version.feature)
        self.assertEqual(bug, version.bug)


if __name__ == '__main__':
    unittest.main()