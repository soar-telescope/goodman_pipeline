import sys

from goodman_pipeline.core import setup_logging
from goodman_pipeline.images import ReduceCCD


if '-h' not in sys.argv and '--help' not in sys.argv and '--version' not in sys.argv:  # pragma: no cover
    setup_logging()


def redccd():  # pragma: no cover
    reduce_ccd = ReduceCCD()
    reduce_ccd()
