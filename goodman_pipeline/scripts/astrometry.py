import sys

from goodman_pipeline.core import setup_logging
from goodman_pipeline.astrometry import Astrometry
from goodman_pipeline.astrometry.utils import get_astrometry_config_args

if '-h' not in sys.argv and '--help' not in sys.argv and '--version' not in sys.argv:  # pragma: no cover
    setup_logging()


def redastrometry():
    args = get_astrometry_config_args()
    astrometry = Astrometry(pixel_scale=args.pixel_scale,
                            pixel_scale_tolerance=args.pixel_scale_tolerance,
                            scale_units=args.scale_units,
                            detection_threshold=args.detection_threshold,
                            initial_fwhm=args.initial_fwhm,
                            downsample_factor=args.downsample_factor,
                            binning_keyword=args.binning_keyword,
                            ra_keyword=args.ra_keyword,
                            dec_keyword=args.dec_keyword,
                            index_directory=args.index_directory,
                            ignore_goodman_vignetting=args.ignore_goodman_vignetting,
                            overwrite=args.overwrite,
                            plots=args.plots,
                            debug=args.debug,
                            verbose=args.verbose)

    astrometry(filename=args.filename, flat_image_filename=args.flat)
