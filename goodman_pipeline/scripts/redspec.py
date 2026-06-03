import sys

from goodman_pipeline.core import setup_logging
from goodman_pipeline.spectroscopy import ReduceSpectroscopy, CreateReferenceLamp

if '-h' not in sys.argv and '--help' not in sys.argv and '--version' not in sys.argv:  # pragma: no cover
    setup_logging()

def redspec():  # pragma: no cover
    reduce_spectroscopy = ReduceSpectroscopy()
    reduce_spectroscopy()


def create_reference_lamp():
    create_ref_lamp = CreateReferenceLamp()
    create_ref_lamp()
