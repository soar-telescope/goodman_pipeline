from __future__ import absolute_import

from . import spectroscopy
from . import images
from . import info
from . import core

import sys

if '-h' not in sys.argv and '--help' not in sys.argv:
    core.setup_logging()
