import logging
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np

from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.stats import mad_std
from astropy.visualization import ZScaleInterval
from astropy.wcs import WCS
from astropy import units as u

from astroquery.gaia import Gaia

from photutils.aperture import aperture_photometry, CircularAperture
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder

from scipy.spatial import cKDTree

from ..core import detect_point_sources, get_goodman_vigneting_mask, validate_fits_file_or_read

log = logging.getLogger(__name__)


class Photometry(object):

    def __init__(self,
                 aperture_radius: float = 4.0,
                 aperture_type: str = 'fixed',
                 detection_threshold: float = 5.0,
                 initial_fwhm: float = 3.0,
                 gaia_sources_limit: int = 5000,
                 gaia_photometry_column: str = '',
                 imaging_filter_keyword: str = 'FILTER',
                 aperture_curve_of_growth: bool = False,
                 plots: bool = False,
                 overwrite: bool = False,
                 debug: bool = False):
        self.filename = None
        self.image_data = None
        self.image_header = None
        self.background_subtracted_data = None
        self.wcs = None
        self.aperture_radius = aperture_radius
        self.aperture_type = aperture_type
        self.detection_threshold = detection_threshold
        self.initial_fwhm = initial_fwhm
        self.gaia_sources_limit = gaia_sources_limit
        self.gaia_photometry_column = gaia_photometry_column
        self.imaging_filter_keyword = imaging_filter_keyword
        self.sources = None
        self.gaia_table = None
        self.gaia_coords = None
        self.photometry_table = None
        self.photometry_table_name = ''
        self.filter_name: str = ''
        self.aperture_curve_of_growth = aperture_curve_of_growth
        self.plots = plots
        self.overwrite = overwrite
        self.debug = debug

        self.zero_point_median = 0
        self.zero_point_std = 0

        log.setLevel(logging.DEBUG if self.debug else logging.INFO)

    def __call__(self, filename):
        self.filename = filename

        self._initial_checks()

        self._subtract_background()

        self._create_mask()

        self.sources = detect_point_sources(data=self.background_subtracted_data,
                                            mask=self.mask,
                                            initial_fwhm=self.initial_fwhm,
                                            detection_threshold=self.detection_threshold,
                                            plots=self.plots)

        if self.aperture_curve_of_growth:
            log.info("Running aperture curve of growth diagnostics.")
            self._estimate_best_aperture_size()
            sys.exit(0)

        if self.aperture_type in ['fixed', 'variable']:
            self._aperture_photometry()

        self._save_photometry_table()

        self._get_gaia_sources()

        self._get_photometric_zeropoint()



        log.info("END")
        return {
            "photometry_table": self.photometry_table_name,
            "measured_sources": len(self.sources),
            "gaia_photometry_column": self.gaia_photometry_column,
            "zero_point_median": self.zero_point_median,
            "zero_point_std": self.zero_point_std,
        }

    def _initial_checks(self):
        self.image_data, self.image_header = validate_fits_file_or_read(self.filename)

        if self.imaging_filter_keyword in self.image_header:
            self.filter_name = self.image_header[self.imaging_filter_keyword]
        else:
            log.error(f"Keyword {self.imaging_filter_keyword} not found in {self.filename}'s header.")
            sys.exit(1)

        log.info("Checking for celestial WCS in the file's header.")
        try:
            self.wcs = WCS(self.image_header)
            if not self.wcs.has_celestial:
                log.error("No Celestial WCS found")
                sys.exit(1)
            else:
                log.info("Celestial WCS found")
        except Exception as e:
            log.error(f"WCS parsing failed: {e}")
            sys.exit(1)
        else:
            log.debug("WCS parsing successful")

        log.debug("Validating aperture type")
        if self.aperture_type not in ['fixed', 'variable']:
            log.error(f"Invalid aperture type: {self.aperture_type}")
            sys.exit(1)

    def _get_gaia_sources(self):
        if self.wcs:
            ny, nx = self.background_subtracted_data.shape
            corners = SkyCoord.from_pixel(xp=[0, nx - 1, nx - 1, 0], yp=[0, 0, ny - 1, ny - 1], wcs=self.wcs, origin=0)
            center = SkyCoord(ra=np.mean(corners.ra.deg) * u.deg, dec=np.mean(corners.dec.deg) * u.deg, frame='icrs')

            radius = max(corners.separation(center)).to(u.deg) * 0.8

            log.info(f"Found radius for querying GAIA: {radius.arcmin} arcmin / {radius.to(u.deg).value} deg")

            log.info(f"Querying Gaia DR3 around RA={center.ra.deg:.5f}, Dec={center.dec.deg:.5f}, radius={radius:.2f}")

            try:
                Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
                query = f"""
                            SELECT TOP {self.gaia_sources_limit}
                                source_id, ra, dec, phot_g_mean_mag, phot_rp_mean_mag, phot_bp_mean_mag
                            FROM gaiadr3.gaia_source
                            WHERE
                                1 = CONTAINS(
                                    POINT('ICRS', ra, dec),
                                    CIRCLE('ICRS', {center.ra.deg}, {center.dec.deg}, {radius.to(u.deg).value})
                                )
                                AND phot_g_mean_mag < 19
                            """

                job = Gaia.launch_job_async(query)
                self.gaia_table = job.get_results()
                log.info(f"Retrieved {len(self.gaia_table)} Gaia sources")

                if len(self.gaia_table) > 0:
                    self.gaia_coords = SkyCoord(ra=self.gaia_table['ra'], dec=self.gaia_table['dec'], unit='deg',
                                                frame='icrs')
                    if self.plots:
                        interval = ZScaleInterval()
                        vmin, vmax = interval.get_limits(self.background_subtracted_data)

                        fig, ax = plt.subplots(
                            subplot_kw={'projection': self.wcs} if self.wcs else {}, figsize=(16, 12)
                        )
                        im = ax.imshow(self.background_subtracted_data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)

                        if self.wcs:
                            ax.set_xlabel('Right Ascension (J2000)')
                            ax.set_ylabel('Declination (J2000)')
                        else:
                            ax.set_xlabel('X Pixel')
                            ax.set_ylabel('Y Pixel')

                        ax.set_title(os.path.basename(self.filename))
                        plt.colorbar(im, ax=ax, label='Pixel value')
                        x_pix, y_pix = self.wcs.world_to_pixel(self.gaia_coords, )
                        #             print(x_pix, y_pix)
                        ax.plot(x_pix, y_pix, 'o', markersize=5, markerfacecolor='none', markeredgecolor='cyan',
                                label='Gaia DR2')
                        ax.legend(loc='upper right')
                        ax.set_xlim(0, nx)
                        ax.set_ylim(0, ny)

                        plt.tight_layout()
                        plt.show()
                else:
                    log.error("Unable to continue due to GAIA table is empty.")

            except Exception as e:
                log.error(f"Gaia query failed: {e}")
        else:
            log.error("Can't query GAIA without having a WCS solution")

    def _detect_sources(self):
        log.info("Performing Fixed Aperture Photometry")
        log.debug("Estimating background for image data")
        background = Background2D(
            data=self.image_data,
            box_size=(64, 64),
            filter_size=(3, 3),
            bkg_estimator=MedianBackground())
        log.info("Subtracting background to image data.")
        self.background_subtracted_data = self.image_data - background.background

        log.debug("Estimating image's noise")
        noise = mad_std(self.background_subtracted_data)

        log.info(
            f"Running DAOStarFinder with fwhm={self.initial_fwhm} and threshold={self.detection_threshold} * noise")
        daofind = DAOStarFinder(fwhm=self.initial_fwhm, threshold=self.detection_threshold * noise)

        self.sources = daofind(self.background_subtracted_data)

        log.info(f"Detected {len(self.sources)} sources.")

    def _estimate_best_aperture_size(self):
        nx, ny = self.background_subtracted_data.shape
        center_x = nx / 2
        center_y = ny / 2

        dx = self.sources['xcentroid'] - center_x
        dy = self.sources['ycentroid'] - center_y
        distance_from_center = np.sqrt(dx**2 + dy**2)

        max_radius = 0.4 * min(nx, ny) / 2

        center_sources = self.sources[distance_from_center < max_radius]

        positions = np.vstack((center_sources['xcentroid'], center_sources['ycentroid'])).T

        tree = cKDTree(positions)

        distances, _ = tree.query(positions, k=2)

        nearest_distances = distances[:, 1]

        minimum_separation = 10

        isolated_mask = nearest_distances > minimum_separation

        isolated_sources = center_sources[isolated_mask]

        bright_center_sources = isolated_sources[np.argsort(isolated_sources['flux'][-5:])]

        test_aperture_sizes = np.arange(1, 15, 1)

        all_curves = {}

        for i, source in enumerate(bright_center_sources):
            position = (source['xcentroid'], source['ycentroid'])
            fluxes = []

            for test_aperture_size in test_aperture_sizes:
                aperture = CircularAperture([position], r=test_aperture_size)
                phot_table = aperture_photometry(self.background_subtracted_data, aperture)
                fluxes.append(phot_table['aperture_sum'][0])

            all_curves[i] = np.array(fluxes)

        all_suggested_radii = []
        for i, fluxes in all_curves.items():
            dflux = np.gradient(fluxes, test_aperture_sizes)

            flux_threshold = 0.1 * np.max(fluxes)

            plateau_idx = np.argmax(dflux < flux_threshold)

            suggested_radius = test_aperture_sizes[plateau_idx]
            all_suggested_radii.append(suggested_radius)
            log.debug(f"Suggested radius: {suggested_radius}")

        overall_suggested_radius = np.median(all_suggested_radii)
        log.info(f"Overall suggested radius: {overall_suggested_radius}")

        if True:

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
            # Selected Sources
            zcale = ZScaleInterval()

            z1, z2 = zcale.get_limits(self.background_subtracted_data)
            image = ax1.imshow(self.background_subtracted_data, origin='lower', cmap='gray', clim=(z1, z2))
            cbar = plt.colorbar(image, ax=ax1)
            cbar.set_label('Counts')


            # Highlight the 5 brightest near the center
            ax1.plot(bright_center_sources['xcentroid'], bright_center_sources['ycentroid'], 'o', markerfacecolor='none', markeredgecolor='r', markersize=5,
                     label='Bright center sources')


            ax1.set_title('Source Selection: Center & Bright & Isolated')
            ax1.set_xlabel('X Pixel')
            ax1.set_ylabel('Y Pixel')
            ax1.legend()
            ax1.grid(False)

            # Curve of Growth
            for i, fluxes in all_curves.items():
                ax2.plot(test_aperture_sizes, fluxes, label=f"Source {i + 1}")
            ax2.axvline(overall_suggested_radius, color='r', linestyle='--', label=f'Overall suggested radius: {overall_suggested_radius}')
            ax2.set_xlabel('Aperture Radius (pixels)')
            ax2.set_ylabel('Flux')
            ax2.set_title('Curve of Growth')
            ax2.legend()
            ax2.grid(True)
            plt.tight_layout()
            plt.show()


        sys.exit(1)

    def _aperture_photometry(self):
        positions = np.transpose((self.sources['xcentroid'], self.sources['ycentroid']))
        apertures = CircularAperture(positions=positions, r=self.aperture_radius)
        log.info(f"Running aperture photometry with fixed aperture of radius: {self.aperture_radius}")
        self.photometry_table = aperture_photometry(data=self.background_subtracted_data, apertures=apertures)

    def _psf_photometry(self):
        pass

    def _save_photometry_table(self):
        log.info("Preparing photometry table to save it.")
        if not self.photometry_table:
            log.error("No photometry table found")
            sys.exit(1)

        if not self.wcs:
            log.warning("No celestial WCS found. Will be unable to convert coordinates to sky coordinates.")
        else:
            ra, dec = self.wcs.all_pix2world(self.photometry_table['xcenter'], self.photometry_table['ycenter'], 0)
            sky_coords = SkyCoord(ra * u.deg, dec * u.deg)

            self.photometry_table['ra'] = sky_coords.ra.deg
            self.photometry_table['ra'].info.format = '%.6f'
            self.photometry_table['dec'] = sky_coords.dec.deg
            self.photometry_table['dec'].info.format = '%.6f'

        self.photometry_table['mag'] = -2.5 * np.log10(self.photometry_table['aperture_sum'])
        self.photometry_table['mag'].info.format = '%.2f'

        self.photometry_table_name = re.sub('.fits', '_phot.csv', self.filename)

        log.debug(f"New table name: {self.photometry_table_name}")

        try:
            self.photometry_table.write(self.photometry_table_name, format='csv', overwrite=self.overwrite)
            log.info(f"Photometry table written to {self.photometry_table_name}.")
        except OSError as e:
            log.debug(f"{e}")
            log.warning(f"Photometry Table  {self.photometry_table_name} already exists.")
            log.info(f"In order to overwrite files use the --overwrite flag. Or use -h or --help for more info.")
            sys.exit(0)

    def _get_photometric_zeropoint(self):
        sky_coords = self.wcs.pixel_to_world(
            self.sources['xcentroid'], self.sources['ycentroid']
        )

        photometry_coords = SkyCoord(sky_coords)

        matched_gaia_indices, separations, _ = match_coordinates_sky(photometry_coords, self.gaia_coords)

        maximum_separation = 1.0 * u.arcsec

        match_mask = separations < maximum_separation

        matched_photometry_sources = self.sources[match_mask]
        matched_gaia_sources = self.gaia_table[matched_gaia_indices[match_mask]]

        flux = matched_photometry_sources['flux']

        instrumental_magnitudes = -2.5 * np.log10(flux)
        if not self.gaia_photometry_column:
            self.gaia_photometry_column = self._get_gaia_band_for_filter()
        log.info(f"Retrieving matched gaia source's magnitudes from column {self.gaia_photometry_column}")
        gaia_magnitudes = matched_gaia_sources[self.gaia_photometry_column]

        corrected_gaia_magnitudes = self._apply_color_transformation_correction(standard_magnitudes=gaia_magnitudes)

        zero_points = corrected_gaia_magnitudes - instrumental_magnitudes

        valid_zero_points = zero_points[np.isfinite(zero_points)]
        self.zero_point_median = np.median(valid_zero_points)
        self.zero_point_std = np.std(valid_zero_points)

        log.info(f"Found zero point {self.zero_point_median:.3f} +/- {self.zero_point_std:.3f}")

    def _get_gaia_band_for_filter(self):
        """
        Maps your filter name to the closest usable Gaia band.

        Returns:
            str: Gaia band column name ('phot_g_mean_mag', 'phot_bp_mean_mag', or 'phot_rp_mean_mag')
        """

        filter_to_gaia_band = {
            # SDSS system
            'u-SDSS': None,
            'g-SDSS': 'phot_bp_mean_mag',
            'r-SDSS': 'phot_g_mean_mag',
            'i-SDSS': 'phot_rp_mean_mag',
            'z-SDSS': None,

            # Johnson/Cousins
            'U': None,
            'B': 'phot_bp_mean_mag',
            'V': 'phot_g_mean_mag',
            'Rc': 'phot_rp_mean_mag',
            'VR': 'phot_g_mean_mag',  # Approximate

            # Bessell
            'U-Bessel': None,
            'B-Bessel': 'phot_bp_mean_mag',
            'V-Bessel': 'phot_g_mean_mag',
            'R-Bessel': 'phot_rp_mean_mag',
            'I-Bessel': 'phot_rp_mean_mag',

            # Stromgren (no Gaia equivalents)
            'u-Stromgren': None,
            'v-Stromgren': None,
            'b-Stromgren': None,
            'y-Stromgren': None,
        }

        gaia_band = filter_to_gaia_band.get(self.filter_name, None)
        if gaia_band is None:
            log.warning(f"The filter {self.filter_name} does not have an equivalent band in gaia, will use 'phot_g_mean_mag' as default.")
        return gaia_band if gaia_band is not None else 'phot_g_mean_mag'

    def _apply_color_transformation_correction(self, standard_magnitudes):

        filter_corrections = {
            "g-SDSS": {"gaia_band": "phot_g_mean_mag", "a": 0.60, "b": -0.15},
            "r-SDSS": {"gaia_band": "phot_g_mean_mag", "a": 0.60, "b": -0.15},
            "i-SDSS": {"gaia_band": "phot_g_mean_mag", "a": 0.60, "b": -0.15},
            "z-SDSS": {"gaia_band": "phot_g_mean_mag", "a": -0.45, "b": 0.12},
            "V": {"gaia_band": "phot_g_mean_mag", "a": 0.017, "b": -0.02},
            "default": {"gaia_band": "phot_g_mean_mag", "a": 0.0, "b": 0.0},
        }

        # self.filter_name
        # self.gaia_table
        return standard_magnitudes
