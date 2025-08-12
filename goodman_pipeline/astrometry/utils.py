import argparse

from importlib.metadata import version


__version__ = version("goodman_pipeline")


def get_astrometry_config_args(arguments=None):
    """Parse command-line arguments for configuring the Astrometry class.

    This function defines and processes command-line arguments used to configure
    an Astrometry object, which handles astrometric solution processing for
    Goodman High Throughput Spectrograph data.

    Args:
        arguments (list, optional): A list of command-line arguments. If None,
            arguments are taken from `sys.argv`. Defaults to None.

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments.

    Parsed Arguments:
        pixel_scale (float): Expected pixel scale in arcsec/pixel (default: 0.15).
        pixel_scale_tolerance (float): Tolerance for scale matching (default: 0.02).
        scale_units (str): Units of the scale, typically 'arcsecperpix' (default: 'arcsecperpix').
        downsample_factor (int): Image downsampling factor for astrometry.net (default: 2).
        binning_keyword (str): FITS header keyword used to determine binning (default: 'CCDSUM').
        overwrite (bool): Overwrite existing astrometry results if present.
        debug (bool): Enable debug mode.
    """
    parser = argparse.ArgumentParser(
        description="Configure parameters for astrometric solution using the Astrometry.net.")

    parser.add_argument('filename', nargs='?', help="Path to file to process.")

    parser.add_argument(
        '--initial-fwhm',
        default=3.0,
        type=float,
        help="Initial FWHM in pixels for source detection. (default=3.0)")

    parser.add_argument(
        '--detection-threshold',
        default=5.0,
        type=float,
        help="Detection threshold above noise factor (default=5.0)")

    parser.add_argument(
        '--flat',
        '--flat-image',
        default=None,
        type=str,
        help="Path to flat image used to build mask. Must match filter and size of file to process. (default=None")

    parser.add_argument(
        '--pixel-scale',
        default=0.15,
        type=float,
        help="Expected pixel scale in arcsec/pixel (default: 0.15).")

    parser.add_argument(
        '--pixel-scale-tolerance',
        default=0.02,
        type=float,
        help="Tolerance for pixel scale matching (default: 0.02).")

    parser.add_argument(
        '--scale-units',
        default='arcsecperpix',
        type=str,
        help="Units for pixel scale, usually 'arcsecperpix' (default: 'arcsecperpix').")

    parser.add_argument(
        '--downsample-factor',
        default=2,
        type=int,
        help="Downsample factor for astrometry.net (default: 2).")

    parser.add_argument(
        '--binning-keyword',
        default='CCDSUM',
        type=str,
        help="FITS header keyword for CCD binning (default: 'CCDSUM').")

    parser.add_argument(
        '--ra-keyword',
        default='OBSRA',
        type=str,
        help="FITS header keyword for observed RA (default: 'OBSRA').")

    parser.add_argument(
        '--dec-keyword',
        default='OBSDEC',
        type=str,
        help="FITS header keyword for observed DEC (default: 'OBSDEC').")

    parser.add_argument(
        '--overwrite',
        action='store_true',
        help="Allow overwriting of existing astrometry results.")

    parser.add_argument(
        '--index-directory',
        default='',
        type=str,
        help="Index directory for astrometry.net (default: None)."
    )

    parser.add_argument(
        '--ignore-goodman-vignetting',
        action='store_true',
        help="Don't do pre source detection and masking, useful for processing images without vignetting. (default=False).")

    parser.add_argument(
        '--plots',
        action='store_true',
        help="Show plots.")

    parser.add_argument(
        '--debug',
        action='store_true',
        help="Enable debug mode.")

    parser.add_argument(
        '--verbose',
        action='store_true',
        help="Enable verbose mode on astrometry.net's logging.")
    parser.add_argument(
        '-v', '--version',
        action='store_true')

    args = parser.parse_args(args=arguments)

    if args.version:
        parser.exit(status=0, message=__version__)


    if not args.filename:
        parser.print_help()
        parser.exit(status=0, message="\nPlease specify a filename to process.\n")

    return args
