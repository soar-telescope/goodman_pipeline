from __future__ import absolute_import

from . import spectroscopy
from . import images
from . import info
from . import core

import logging
import sys
import datetime


def setup_logging():

    LOG_FILENAME = 'goodman_log.txt'

    logging_level = logging.INFO
    if '--debug' in sys.argv:
        FORMAT = '[%(asctime)s][%(levelname)8s][%(module)s.%(funcName)s:%(lineno)d]: %(message)s'
        logging_level = logging.DEBUG
    else:
        FORMAT = '[%(asctime)s][%(levelname)8s]: %(message)s'
        logging_level = logging.INFO
    DATE_FORMAT = '%I:%M:%S%p'

    formatter = logging.Formatter(fmt=FORMAT, datefmt=DATE_FORMAT)
    logging.basicConfig(level=logging_level,
                        format=FORMAT,
                        datefmt=DATE_FORMAT)

    log = logging.getLogger(__name__)

    file_handler = logging.FileHandler(filename=LOG_FILENAME)
    file_handler.setFormatter(fmt=formatter)
    file_handler.setLevel(level=logging_level)
    log.addHandler(file_handler)

    log.info("Starting Goodman HTS Pipeline Log")
    log.info("Local Time    : {:}".format(
        datetime.datetime.now()))
    log.info("Universal Time: {:}".format(
        datetime.datetime.utcnow()))


if '-h' not in sys.argv and '--help' not in sys.argv:
    setup_logging()
