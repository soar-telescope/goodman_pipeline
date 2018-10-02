from __future__ import absolute_import

from unittest import TestCase, skip
from ccdproc.core import slice_from_string

class TestCcdprocMethods(TestCase):

    def test_slice_from_string(self):
        string = '[1:100,15:85]'
        result = tuple([slice(14, 85, None), slice(0, 100, None)])
        python_slice = slice_from_string(string=string, fits_convention=True)
        self.assertEqual(result, python_slice)
        # self.fail()
