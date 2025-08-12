import glob
import os
import logging
import numpy as np
import re
import sys
import shutil
import subprocess

from astropy import units as u
from astropy.io import fits
from ccdproc import CCDData
from pathlib import Path
from typing import Union

from ..core import (create_xyls_table,
                    detect_point_sources,
                    get_vigneting_mask,
                    subtract_background_from_image_data,
                    validate_fits_file_or_read)

log = logging.getLogger(__name__)


class Astrometry(object):

    def __init__(self,
                 pixel_scale: float = 0.15,
                 pixel_scale_tolerance: float = 0.02,
                 scale_units: str = 'arcsecperpix',
                 downsample_factor: int = 2,
                 binning_keyword: str = 'CCDSUM',
                 ra_keyword: str = 'OBSRA',
                 dec_keyword: str = 'OBSDEC',
                 imaging_filter_keyword: str = 'FILTER',
                 index_directory: str= '',
                 ignore_goodman_vignetting: bool = False,
                 plots: bool = False,
                 overwrite: bool = False,
                 verbose: bool = False,
                 debug: bool = False):
        self.filename = None
        self.target_file = None
        self.flat_image_filename = None
        self.image_data = None
        self.image_header = None
        self.flat_image_data = None
        self.flat_image_header = None
        self.pixel_scale = pixel_scale * u.arcsec / u.pix
        self.serial_binning = 1
        self.parallel_binning = 1
        self.pixel_scale_tolerance = pixel_scale_tolerance * u.arcsec / u.pix
        self.scale_units = scale_units
        self.scale_low = pixel_scale
        self.scale_high = pixel_scale
        self.downsample_factor = downsample_factor
        self.ra_keyword = ra_keyword
        self.dec_keyword = dec_keyword
        self._ra = ''
        self._dec = ''
        self._radius = 0.0
        self.solve_field_executable = "solve-field"
        self.solve_field_full_path = None
        self.index_directory = index_directory
        self.xy_sources_position = ''
        self.ignore_goodman_vignetting = ignore_goodman_vignetting
        self.plots = plots
        self.overwrite = overwrite
        self.verbose = verbose
        self.debug = debug

        self.binning_keyword = binning_keyword
        self.imaging_filter_keyword = imaging_filter_keyword

        self._new_files = {}

        log.setLevel(logging.DEBUG if self.debug else logging.INFO)


    def __call__(self, filename: str, flat_image_filename: Union[str, None] = None):
        log.info(f"Processing file {filename}")
        self.filename = filename
        self.flat_image_filename = flat_image_filename

        self._initial_checks()

        self._set_parameters()

        if not self.ignore_goodman_vignetting:
            self.xy_sources_position = 'something-other-than-empty'

        return_code, solve_field_full_logs = self.astrometry_net__solve_field(
            filename=self.filename,
            ra=self._ra,
            dec=self._dec,
            radius=self._radius.to(u.deg).value,
            scale_low=self.scale_low.value,
            scale_high=self.scale_high.value,
            scale_units=self.scale_units,
            downsample=self.downsample_factor,
            solve_field_executable=self.solve_field_full_path,
            index_directory=self.index_directory,
            xy_sources_position=self.xy_sources_position,
            overwrite=self.overwrite,
            verbose=self.verbose)

        log.debug(f"Astrometry.net exited with code {return_code}")
        if return_code in [255]:
            sys.exit(return_code)

        self._create_file_with_wcs()

        self._detect_new_files()

        return {
            'full_logs': solve_field_full_logs,
            'new_files': self._new_files,
        }


    def _initial_checks(self):
        log.info("Running input checks")
        self.image_data, self.image_header = validate_fits_file_or_read(filename=self.filename)

        if self.flat_image_filename is not None:
            self.flat_image_data, self.flat_image_header = validate_fits_file_or_read(filename=self.flat_image_filename)
            log.debug(
                f"Image's filter {self.image_header[self.imaging_filter_keyword]} vs flat filter {self.flat_image_header[self.imaging_filter_keyword]}")

            if (self.image_header[self.imaging_filter_keyword] != self.flat_image_header[self.imaging_filter_keyword] or
                    self.image_data.shape != self.flat_image_data.shape):
                log.error(f"Flat image provided is not compatible with science data.")
                sys.exit(1)

        log.debug(f"Validating that executable {self.solve_field_executable} exists")
        self.solve_field_full_path = shutil.which(self.solve_field_executable)
        if self.solve_field_full_path is None or not os.path.exists(self.solve_field_full_path):
            log.error(f"Unable to locate executable {self.solve_field_executable}")
            sys.exit(1)
        else:
            log.debug(f"Executable  {self.solve_field_executable} found at {self.solve_field_full_path}")
        if self.index_directory:
            log.debug(f"Validate --index-directory {self.index_directory}")
            if not os.path.isdir(self.index_directory) or not os.path.exists(self.index_directory):
                log.error(f"Index directory {self.index_directory} does not exist")
                sys.exit(1)
            else:
                index_directory_length = len(os.listdir(self.index_directory))
                if index_directory_length > 0:
                    log.debug(f"Index directory {self.index_directory} exists and contains {index_directory_length} files.")
                else:
                    log.error(f"Index directory {self.index_directory} is empty")
                    sys.exit(1)
        else:
            log.debug(f"No custom --index-directory specified")

        log.info("All inputs are valid")

    def _set_parameters(self):
        log.info("Calculating parameters based on input.")
        self.serial_binning, self.parallel_binning = [int(x) for x in self.image_header[self.binning_keyword].split()]
        log.debug(f"Found serial and parallel binning {self.serial_binning} x {self.parallel_binning}")
        if self.serial_binning != self.parallel_binning:
            log.warning(f"Binning is asymmetric {self.serial_binning}x{self.parallel_binning}")
        self.scale_low = self.pixel_scale * self.serial_binning - self.pixel_scale_tolerance
        log.debug(f"Set scale low to {self.scale_low}")
        self.scale_high = self.pixel_scale * self.serial_binning + self.pixel_scale_tolerance
        log.debug(f"Set scale high to {self.scale_high}")

        log.debug(f"Updating RA and DEC from image's header.")
        self._ra = self.image_header[self.ra_keyword]
        self._dec = self.image_header[self.dec_keyword]

        log.debug(f"Finding image's radius.")
        image_larger_side = np.max(self.image_data.shape) * u.pix
        self._radius = self.pixel_scale * image_larger_side * self.serial_binning
        log.info(f"Setting RA {self._ra}, DEC {self._dec} and radius {self._radius.to(u.arcmin)} or {self._radius.to(u.deg)}")

    def _detect_new_files(self):
        all_matching_files = [_file  for _file in glob.glob(re.sub('.fits', '*', self.filename)) if _file != self.filename]

        extension_to_key = {
            ".axy": "augmented_xylist",
            ".corr": "matched_stars",
            ".match": "quad_match_info",
            ".rdls": "reference_catalog",
            ".solved": "solved_flag",
            ".wcs": "wcs_header",
            ".xyls": "extracted_sources",
            "_wcs.fits": "wcs_fits_image",
            "-indx.png": "index_overlay_image",
            "-ngc.png": "ngc_overlay_image",
            "-objs.png": "object_overlay_image"
        }
        self.new_files = {}

        for _file in all_matching_files:
            for suffix, key in extension_to_key.items():
                if _file.endswith(suffix):
                    self.new_files[key] = _file
                    break

    def _create_file_with_wcs(self):
        """Creates a new FITS file containing WCS information.

        Depending on the file type:
            - For `.fits` inputs: renames a corresponding `.new` file to `_wcs.fits`.
            - For `.xyls` inputs: merges the original image header with a `.wcs`
              header file and saves as `_wcs.fits`.

        Raises:
            FileNotFoundError: If the required `.new` or `.wcs` file does not exist.
            OSError: If file renaming or writing fails.
        """
        path = Path(self.target_file)

        if path.suffix == ".fits":
            self._rename_new_fits_file(path)
        elif path.suffix == ".xyls":
            self._merge_headers_and_save(path)
        else:
            log.warning(f"Unsupported file type: {path.suffix}")

    def _rename_new_fits_file(self, path: Path):
        """Renames an existing `.new` file to `_wcs.fits`."""
        new_file = path.with_suffix(".new")
        new_wcs_file = path.with_name(path.stem + "_wcs.fits")

        log.debug(f"Looking for new file {new_file}")

        if new_file.exists() and new_file.is_file():
            log.info("A '*.new' file exists. Renaming to '_wcs.fits'.")
            new_file.rename(new_wcs_file)
        else:
            raise FileNotFoundError(f"Unable to find expected file: {new_file}")

    def _merge_headers_and_save(self, path: Path):
        """Merges the image header with WCS header and saves as `_wcs.fits`."""
        wcs_header_file = path.with_suffix(".wcs")
        new_file_name = path.with_name(path.stem + "_wcs.fits")

        if not wcs_header_file.exists():
            raise FileNotFoundError(f"WCS header file not found: {wcs_header_file}")

        log.debug(f"Reading WCS header from {wcs_header_file}")
        with fits.open(wcs_header_file) as hdul:
            wcs_header = hdul[0].header

        new_header = fits.Header(self.image_header)
        new_header.extend(wcs_header)

        log.info(f"Saving new file with WCS solution to {new_file_name}")
        new_hdu = CCDData(data=self.image_data, meta=new_header, unit="adu")
        new_hdu.write(new_file_name, overwrite=self.overwrite)

    @staticmethod
    def astrometry_net__solve_field(
            filename: str,
            ra: Union[str, float],
            dec: Union[str, float],
            radius: float,
            scale_low: float,
            scale_high: float,
            scale_units: str='arsecperpix',
            downsample: int=2,
            solve_field_executable: str="solve-field",
            index_directory: str='',
            xy_sources_position: str='',
            overwrite: bool=False,
            verbose: bool=False):
        options = {
            "overwrite": overwrite,
            "scale-low": scale_low,
            "scale-high": scale_high,
            "scale-units": scale_units,
            "ra": ra,
            "dec": dec,
            "radius": radius,
            "xyls": xy_sources_position,
            "index-dir": index_directory,
            "downsample": downsample,
            "verbose": verbose
        }
        options_str = " ".join(
            f"--{key}" if isinstance(value, bool) and value
            else f"--{key} {value}"
            for key, value in options.items()
            if value is not None and value != "" and (not isinstance(value, bool) or value)
        )
        command = f"{solve_field_executable} {options_str} {filename}"

        log.debug(f"command: {command}")

        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)

        process_logs = []
        for line in process.stdout:
            log.info(f"{os.path.basename(solve_field_executable)}: {line.rstrip()}")
            process_logs.append(line)

        process.wait()

        return_code = process.returncode
        log.debug(f"Process 'solve-field' exit code: {return_code}")
        if return_code in [255]:
            log.error(f"Astrometry.net's solve-field failed to solve for {filename}")
            sys.exit(return_code)

        full_logs = "".join(process_logs)
        return return_code, full_logs
