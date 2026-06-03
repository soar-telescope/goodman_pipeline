import sys

from goodman_pipeline.core import setup_logging
from goodman_pipeline.photometry.utils import get_photometry_config_args
from goodman_pipeline.photometry import Photometry

if '-h' not in sys.argv and '--help' not in sys.argv and '--version' not in sys.argv:  # pragma: no cover
    setup_logging()


def redphotometry():
    args = get_photometry_config_args()

    photometry = Photometry(aperture_radius=args.aperture_radius,
                            # aperture_type=args.aperture_type,
                            detection_threshold=args.detection_threshold,
                            initial_fwhm=args.initial_fwhm,
                            gaia_sources_limit=args.gaia_sources_limit,
                            gaia_photometry_column=args.gaia_photometry_column,
                            imaging_filter_keyword=args.imaging_filter_keyword,
                            aperture_curve_of_growth=args.aperture_curve_of_growth,
                            disable_mask_creation=args.disable_mask_creation,
                            plots=args.plots,
                            overwrite=args.overwrite,
                            debug=args.debug)

    photometry(filename=args.filename, flat_image_filename=args.flat)
