#Felipe's template code libraries
from astropy              import units as u
from astropy              import wcs
from astropy.io           import fits
from astropy.coordinates  import SkyCoord
from astropy.modeling     import models, fitting
from astropy.stats        import sigma_clipped_stats, SigmaClip

from photutils            import DAOStarFinder, CircularAperture, aperture_photometry, CircularAnnulus
from photutils.background import Background2D, MedianBackground
from photutils.utils      import calc_total_error

from ..core import (phot_do_aperture_phot,
                    phot_display,
                    phot_display_table,
                    phot_get_pixelscale,
                    phot_get_fwhm,
                    phot_table_xy2sky,
                    phot_get_background2d,
                    phot_get_background_statistics,
                    phot_find_stars)

import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
#import seaborn
#---------------------------------

class AperturePhotometry(object):

    def __init__(self):
        # Initialize your class attributes if needed
        pass

    def process_single_file(self, file_path, args):
        """
        This function performs photometry analysis on a single FITS file.

        Args:
            file_path (str): Path to the FITS file.
            args (argparse.Namespace): Namespace object containing command-line arguments.

        Returns:
            None
        """

        # Open the FITS file
        img, hdr = fits.getdata(file_path, header=True)

        # Call helper functions for processing
        stars = phot_find_stars(img, fwhm=args.fwhm, sigma_thresh=args.sigma_threshold)

        print("PRE FWHM FIT")

        fwhm_fit = phot_get_fwhm(img, stars, source_id=args.source_id, pixscale=args.pixscale,
                                 stamp_radius=args.stamp_radius, psf_model=args.psf_model,
                                 plot_results=args.plot)
        print("POST FWHM FIT")

        if fwhm_fit is None:
            print("Warning: FWHM estimation failed.")
        print("* FWHM fit: ",fwhm_fit)

        bkg = phot_get_bkg2d(img, plot_results=args.plot)
        positions = (stars['xcentroid'], stars['ycentroid'])
        phot_table = phot_do_aperture_phot(img, positions, bkg=bkg, gain=args.gain, r_in=args.r_in,
                                           r_out=args.r_out)

        # Optional: Display results tables
        if args.plot:
            displ_table(stars)
            displ_table(phot_table)


if __name__ == '__main__':  # pragma: no cover
    sys.exit('This can not be run on its own.')