from __future__ import absolute_import

from ..data_classifier import DataClassifier


def test_data_classifier():
    data_classifier = DataClassifier()
    assert isinstance(data_classifier, DataClassifier)

