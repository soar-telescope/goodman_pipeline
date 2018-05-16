from __future__ import absolute_import

from . import spectroscopy
from . import images
from . import info
from . import core
from . import tools

from .core import setup_logging

import os
import sys

# set version information
# Get metadata from setup.cfg
try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser
conf = ConfigParser()

conf.read([os.path.join(os.path.dirname(__file__), '..', 'setup.cfg')])

metadata = dict(conf.items('metadata'))

__version__ = metadata['version']

__license__ = metadata['license']

__description__ = metadata['description']

__long_description__ = metadata['long_description']

__author__ = metadata['author']

__author_email__ = metadata['author_email']


if '-h' not in sys.argv and '--help' not in sys.argv:
    setup_logging()
