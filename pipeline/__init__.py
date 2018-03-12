from __future__ import absolute_import

from . import spectroscopy
from . import images
from . import info
from . import core

import logging
import sys

if '--debug' in sys.argv:
    FORMAT = '%(levelname)s: %(asctime)s:%(module)s.%(funcName)s: %(message)s'
else:
    FORMAT = '%(levelname)s: %(asctime)s: %(message)s'
DATE_FORMAT = '%I:%M:%S%p'


logging.basicConfig(level=logging.INFO,
                    format=FORMAT,
                    datefmt=DATE_FORMAT)

log = logging.getLogger(__name__)
