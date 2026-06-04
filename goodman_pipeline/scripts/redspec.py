import sys

from goodman_pipeline.core import setup_logging
from goodman_pipeline.spectroscopy import ReduceSpectroscopy, CreateReferenceLamp
from goodman_pipeline.spectroscopy.utils import get_create_reference_lamp_args

if '-h' not in sys.argv and '--help' not in sys.argv and '--version' not in sys.argv:  # pragma: no cover
    setup_logging()


def redspec():  # pragma: no cover
    reduce_spectroscopy = ReduceSpectroscopy()
    reduce_spectroscopy()


def create_reference_lamp():
    args = get_create_reference_lamp_args()
    create_ref_lamp = CreateReferenceLamp(
        comparison_lamp_full_path=args.comparison_lamp,
        reference_lamp_full_path=args.reference_lamp,
        comparison_lamp_intensity_start=args.comp_intensity_start,
        comparison_lamp_intensity_end=args.comp_intensity_end,
        reference_lamp_wavelength_start=args.ref_wavelength_start,
        reference_lamp_wavelength_end=args.ref_wavelength_end,
        plots_theme=args.plots_theme,
        screen_size=args.screen_size,
        debug=args.debug)
    create_ref_lamp()
