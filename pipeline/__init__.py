from __future__ import absolute_import

from . import spectroscopy
from . import images
from . import info
from . import core

from .core import setup_logging

import sys

if '-h' not in sys.argv and '--help' not in sys.argv:
    setup_logging()
