#!/usr/bin/env python3
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from goodman_pipeline.core import setup_logging
from goodman_pipeline.images import ReduceCCD
from goodman_pipeline.spectroscopy import ReduceSpectroscopy

from goodman_pipeline.astrometry import Astrometry
from goodman_pipeline.astrometry.utils import get_astrometry_config_args

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
                            index_directory=args.index_directory,
                            overwrite=args.overwrite,
                            debug=args.debug,
                            verbose=args.verbose)

    astrometry(filename=args.filename, return_json=False)


def redphotometry():
    # args = get_photometry_config_args()
    # photometry = Photometry()
    # photometry(filename=args.filename)
    pass
