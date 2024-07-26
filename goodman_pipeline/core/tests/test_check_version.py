
__author__ = 'Bruno Quint'

import os
import unittest
import requests

from importlib.metadata import version

from ..core import check_version

__version__ = version('goodman_pipeline')


class TestVersionChecker(unittest.TestCase):

    def test_get_last(self):
        try:
            v = check_version.get_last()
            self.assertRegex(v.base_version, '^(\*|\d+(\.\d+){0,2}(\.\*)?)$')
            # self.assertEqual(v, __version__)
        except ConnectionRefusedError:  # pragma: no cover
            pass
        except requests.exceptions.ConnectionError:  # pragma: no cover
            pass

    def test_get_last_no_token(self):
        try:
            v = check_version.get_last(github_api_token='NONEXISTANTVAR')
            self.assertRegex(v.base_version, '^(\*|\d+(\.\d+){0,2}(\.\*)?)$')
            # self.assertEqual(v, __version__)
        except ConnectionRefusedError:  # pragma: no cover
            pass
        except requests.exceptions.ConnectionError:  # pragma: no cover
            pass
        except KeyError:  # pragma: no cover
            pass

    def test_get_last_token(self):
        os.environ['FAKETOKEN'] = 'ThisIsNotARealToken'
        try:
            self.assertRaises(ConnectionRefusedError,
                              check_version.get_last,
                              'FAKETOKEN')
        except requests.exceptions.ConnectionError:
            pass


    def test_am_i_updated(self):
        try:
            self.assertTrue(check_version.am_i_updated(__version__))
            self.assertFalse(check_version.am_i_updated('v0.0.0'))
        except ConnectionRefusedError:  # pragma: no cover
            pass
        except requests.exceptions.ConnectionError:  # pragma: no cover
            pass
