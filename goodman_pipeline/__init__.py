from __future__ import absolute_import

from importlib.metadata import version

from . import spectroscopy
from . import images
from . import core

from .core import setup_logging

__version__ = version('goodman_pipeline')
