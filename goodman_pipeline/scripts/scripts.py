#!/usr/bin/env python3
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from goodman_pipeline.core import setup_logging
from goodman_pipeline.images import ReduceCCD
from goodman_pipeline.spectroscopy import ReduceSpectroscopy

from goodman_pipeline.astrometry import Astrometry
from goodman_pipeline.astrometry.utils import get_astrometry_config_args
from goodman_pipeline.photometry.utils import get_photometry_config_args

from goodman_pipeline.photometry import Photometry

import sys

if '-h' not in sys.argv and \
                '--help' not in sys.argv and \
                '--version' not in sys.argv:  # pragma: no cover
    setup_logging()

def redccd():  # pragma: no cover
    reduce_ccd = ReduceCCD()
    reduce_ccd()


def redspec():  # pragma: no cover
    reduce_spectroscopy = ReduceSpectroscopy()
    reduce_spectroscopy()


def redastrometry():
    args = get_astrometry_config_args()
    astrometry = Astrometry(pixel_scale=args.pixel_scale,
                            pixel_scale_tolerance=args.pixel_scale_tolerance,
                            scale_units=args.scale_units,
                            downsample_factor=args.downsample_factor,
                            binning_keyword=args.binning_keyword,
                            ra_keyword=args.ra_keyword,
                            dec_keyword=args.dec_keyword,
                            index_directory=args.index_directory,
                            overwrite=args.overwrite,
                            debug=args.debug,
                            verbose=args.verbose)

    astrometry(filename=args.filename, return_json=False)


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

    results = photometry(filename=args.filename)
