import logging
import os
import re
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.table import Column
from astropy.table import hstack
from astropy.visualization import ZScaleInterval
from astropy.wcs import WCS
from astropy import units as u

from astroquery.gaia import Gaia

from photutils.aperture import aperture_photometry, CircularAperture
from photutils.background import Background2D, MedianBackground

from scipy.spatial import cKDTree

from typing import Union

from ..core import detect_point_sources, get_vigneting_mask, validate_fits_file_or_read

matplotlib._log.setLevel(logging.WARNING)
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
                 exposure_time_keyword: str = 'EXPTIME',
                 aperture_curve_of_growth: bool = False,
                 disable_mask_creation: bool = False,
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
        self.exposure_time_keyword = exposure_time_keyword
        self.sources = None
        self.gaia_table = None
        self.gaia_coords = None
        self.photometry_table = None
        self.photometry_table_name = ''
        self.filter_name: str = ''
        self.aperture_curve_of_growth = aperture_curve_of_growth
        self.disable_mask_creation = disable_mask_creation
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

        self._get_gaia_sources()

        self._get_photometric_zeropoint()

        self._save_photometry_table()

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

                        ax.set_title(f"Sources Downloaded from GAIA DR3")
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

    def _subtract_background(self):
        log.debug("Estimating background for image data")
        background = Background2D(
            data=self.image_data,
            box_size=(64, 64),
            filter_size=(3, 3),
            bkg_estimator=MedianBackground())
        log.info("Subtracting background to image data.")
        self.background_subtracted_data = self.image_data - background.background

    def _create_mask(self):
        if not self.disable_mask_creation:
            self.mask = get_vigneting_mask(self.background_subtracted_data, flat_data=self.flat_image_data)
        else:
            log.info("Mask creation disabled, remove --disable-mask-creation if this is not what you want.")
            self.mask = np.zeros_like(self.background_subtracted_data, dtype=np.bool)

        if self.plots:

            scale = ZScaleInterval()

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
            z1, z2 = scale.get_limits(self.background_subtracted_data)

            ax1.set_title("Original Image")
            ax1.imshow(self.background_subtracted_data, clim=(z1, z2), cmap='gray')

            masked_data = self.background_subtracted_data.copy()
            masked_data[self.mask] = np.nan

            # z1, z2 = scale.get_limits(masked_data)
            ax2.set_title(f"Masked Image {'(mask creation disabled)' if self.disable_mask_creation else ''}")
            ax2.imshow(masked_data, clim=(z1, z2), cmap='gray')
            plt.show()


    def _estimate_best_aperture_size(self):
        nx, ny = self.background_subtracted_data.shape
        center_x = nx / 2
        center_y = ny / 2

        dx = self.sources['xcentroid'] - center_x
        dy = self.sources['ycentroid'] - center_y
        distance_from_center = np.sqrt(dx**2 + dy**2)

        max_radius = 0.6 * min(nx, ny) / 2

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

        if self.plots:

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
        if not self.photometry_table:
            log.error("It appears that aperture photometry was unsuccessful")
            sys.exit(1)

        self.photometry_table = self.photometry_table[self.photometry_table['aperture_sum'] > 0]

        if not self.wcs:
            log.warning("No celestial WCS found. Will be unable to convert coordinates to sky coordinates.")
        else:
            ra, dec = self.wcs.all_pix2world(self.photometry_table['xcenter'], self.photometry_table['ycenter'], 0)
            sky_coords = SkyCoord(ra * u.deg, dec * u.deg)

            self.photometry_table['ra'] = sky_coords.ra.deg
            self.photometry_table['ra'].info.format = '%.6f'
            self.photometry_table['dec'] = sky_coords.dec.deg
            self.photometry_table['dec'].info.format = '%.6f'

        exposure_time = float(self.image_header[self.exposure_time_keyword])
        log.debug(f"Dividing aperture sum by exposure time: {exposure_time} seconds")

        self.photometry_table['mag_inst'] = -2.5 * np.log10(self.photometry_table['aperture_sum'] / exposure_time)
        self.photometry_table['mag_inst'].info.format = '%.2f'

        # self.photometry_table.pprint(max_lines=-1, max_width=-1)

    def _psf_photometry(self):
        pass

    def _save_photometry_table(self):
        log.info("Preparing photometry table to save it.")

        self.photometry_table_name = re.sub('.fits', '_phot.csv', self.filename)

        log.debug(f"New table name: {self.photometry_table_name}")

        try:
            # self.photometry_table.pprint(max_lines=-1, max_width=-1)
            self.photometry_table.write(self.photometry_table_name, format='csv', overwrite=self.overwrite)
            log.info(f"Photometry table written to {self.photometry_table_name}.")
        except OSError as e:
            log.debug(f"{e}")
            log.error(f"Photometry Table  {self.photometry_table_name} already exists.")
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

        matched_table = hstack([matched_photometry_sources, matched_gaia_sources])

        filtered_sources = matched_table[
            (matched_table['bp_rp_color'].value > 0.3) &
            (matched_table['bp_rp_color'].value < 2.5) &
            (matched_table['phot_g_mean_mag'].value < 18)
        ]


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

            ax.set_title(f"Matched and Selected Sources in {os.path.basename(self.filename)}")
            plt.colorbar(im, ax=ax, label='Pixel value')
            matched_coords = SkyCoord(ra=filtered_sources['ra'].to_value('deg').filled(np.nan), dec=filtered_sources['dec'].to_value('deg').filled(np.nan), unit='deg',
                                   frame='icrs')
            x_pix, y_pix = self.wcs.world_to_pixel(matched_coords, )
            #             print(x_pix, y_pix)
            ax.plot(x_pix, y_pix, 'o', markersize=5, markerfacecolor='none', markeredgecolor='cyan',
                    label='Gaia DR2')
            ax.legend(loc='upper right')
            # ax.set_xlim(0, nx)
            # ax.set_ylim(0, ny)

            plt.tight_layout()
            plt.show()

        if not self.gaia_photometry_column:
            self.gaia_photometry_column = self._get_gaia_band_for_filter()
        log.info(f"Retrieving matched gaia source's magnitudes from column {self.gaia_photometry_column}")

        log.info(f"Calculating instrumental zero point as Gaia's {self.gaia_photometry_column} minus instrumental magnitude in {self.filter_name}")
        zero_points = filtered_sources[self.gaia_photometry_column].value - filtered_sources['mag'].value

        valid_zero_points = zero_points[np.isfinite(zero_points)]
        self.zero_point_median = np.median(valid_zero_points)
        self.zero_point_std = np.std(valid_zero_points)

        self.photometry_table['zeropoint_inst'] = [self.zero_point_median] * len(self.photometry_table)
        self.photometry_table['zeropoint_std_inst'] = [self.zero_point_std] * len(self.photometry_table)

        log.info(f"Found instrumental zero point {self.zero_point_median:.3f} +/- {self.zero_point_std:.3f}")

        if self.filter_name in ['g-SDSS', 'r-SDSS', 'i-SDSS', 'z-SDSS', 'B', 'V']:
            log.info(f"Obtaining zero point converting GAIA's magnitude to {self.filter_name}'s system.")

            parameters = {
                'g-SDSS': [-0.2199, 0.6365, 0.1548, -0.0064],
                'r-SDSS': [0.09837, -0.08592, -0.1907, 0.1701],
                'i-SDSS': [0.293, -0.6404, 0.09609, 0.002104],
                'z-SDSS': [0.4619, -0.8992, 0.08271, -0.005029],
                'B': [-0.01448, 0.6874, 0.3604, 0.06718],
                'V': [0.02704, -0.01424, 0.2156, -0.01426]
            }

            params = parameters[self.filter_name]

            gaia_sloan_mag = self.gaia_to_filter_conversion(
                gaia_g=filtered_sources['phot_g_mean_mag'].value,
                gaia_bp_rp=filtered_sources['bp_rp_color'].value,
                param_0=params[0],
                param_1=params[1],
                param_2=params[2],
                param_3=params[3])
            calibrated_zero_points = gaia_sloan_mag - filtered_sources['mag'].value

            calibrated_zero_point_median = np.median(calibrated_zero_points)
            calibrated_zero_point_std = np.std(calibrated_zero_points)

            self.photometry_table['zeropoint'] = [calibrated_zero_point_median] * len(self.photometry_table)
            self.photometry_table['zeropoint_std'] = [calibrated_zero_point_std] * len(self.photometry_table)

            log.info(f"Calibrated zero point median: {calibrated_zero_point_median:.3f} +/- {calibrated_zero_point_std:.3f} for {self.filter_name}")

            calibrated_mag_column = Column(self.photometry_table['mag_inst'] + calibrated_zero_point_median, name='mag')

            self.photometry_table.add_column(calibrated_mag_column)

        else:
            log.error(f"Filter {self.filter_name} does not have system convertion for absolute magnitude.")


    @staticmethod
    def gaia_to_filter_conversion(gaia_g, gaia_bp_rp, param_0, param_1, param_2, param_3):
        return gaia_g + param_0 + param_1 * gaia_bp_rp + param_2 * gaia_bp_rp ** 2 + param_3 * gaia_bp_rp ** 3


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
