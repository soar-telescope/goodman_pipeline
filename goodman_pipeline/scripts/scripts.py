#!/usr/bin/env python3
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from goodman_pipeline.core import setup_logging
from goodman_pipeline.images import ReduceCCD
from goodman_pipeline.spectroscopy import ReduceSpectroscopy

import sys

if '-h' not in sys.argv and \
                '--help' not in sys.argv and \
                '--version' not in sys.argv:  # pragma: no cover
    setup_logging()

def redccd():  # pragma: no cover
    reduce_ccd = ReduceCCD()
    reduce_ccd()


def redspec():  # pragma: no cover
    reduce_spectroscopy = ReduceSpectroscopy()
    reduce_spectroscopy()
