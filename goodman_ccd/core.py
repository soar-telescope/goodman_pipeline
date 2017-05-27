from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import sys
import time
import glob
import random
import logging
import ccdproc
from ccdproc import CCDData
from astropy.coordinates import EarthLocation
from astropy.time import Time, TimeDelta
from astroplan import Observer
from astropy import units as u
import numpy as np
from scipy import signal

log = logging.getLogger('goodmanccd.core')


def convert_time(in_time):
    """Converts time to seconds since epoch

    Args:
        in_time (str): time obtained from header's keyword DATE-OBS

    Returns:
        time in seconds since epoch

    """
    return time.mktime(time.strptime(in_time, "%Y-%m-%dT%H:%M:%S.%f"))


def fix_duplicated_keywords(night_dir):
    """Remove duplicated keywords

    There are some cases when the raw data comes with duplicated keywords.
    The origin has not been tracked down. The solution is to identify the
    duplicated keywords and the remove all but one from the end backwards.

    Args:
        night_dir (str): The full path for the raw data location

    """
    files = glob.glob(os.path.join(night_dir, '*.fits'))
    # Pick a random file to find duplicated keywords
    random_file = random.choice(files)
    ccd = CCDData.read(random_file)
    header = ccd.header
    # Put the duplicated keywords in a list
    multiple_keys = []
    for keyword in header.keys():
        if keyword != '':
            if header.count(keyword) > 1:
                if keyword not in multiple_keys:
                    multiple_keys.append(keyword)
    if multiple_keys != []:
        for image_file in files:
            try:
                ccd = CCDData.read(image_file)
                for keyword in multiple_keys:
                    while ccd.header.count(keyword) > 1:
                        ccd.header.remove(keyword,
                                          ccd.header.count(keyword) - 1)
                log.warning('Overwriting file with duplicated keywords removed')
                log.info('File %s overwritten', image_file)
                ccd.write(image_file, clobber=True)
            except IOError as error:
                log.error(error)


def ra_dec_to_deg(right_ascension, declination):
    """Converts right ascension and declination to degrees

    Args:
        right_ascension (str): Right ascension in the format hh:mm:ss.sss
        declination (str): Declination in the format dd:mm:ss.sss

    Returns:
        right_ascension_deg (float): Right ascension in degrees
        declination_deg (float): Declination in degrees

    """
    right_ascension = right_ascension.split(":")
    declination = declination.split(":")
    # RIGHT ASCENTION conversion
    right_ascension_deg = (float(right_ascension[0])
                           + (float(right_ascension[1])
                              + (float(right_ascension[2]) / 60.)) / 60.) * (360. / 24.)
    # DECLINATION conversion
    # sign = float(declination[0]) / abs(float(declination[0]))
    if float(declination[0]) == abs(float(declination[0])):
        sign = 1
    else:
        sign = -1
    declination_deg = sign * (abs(float(declination[0]))
                              + (float(declination[1])
                                 + (float(declination[2]) / 60.)) / 60.)
    return right_ascension_deg, declination_deg


def print_spacers(message):
    """Miscelaneous function to print uniform spacers

    Prints a spacer of 80 columns with  and 3 rows height. The first and last
    rows contains the symbol "=" repeated 80 times. The middle row contains the
    message centered and the extremes has one single "=" symbol.
    The only functionality of this is aesthetic.

    Args:
        message (str): A message to be printed

    Returns:
        True (bool): A True value

    """

    columns = 80
    if len(message) % 2 == 1 and int(columns) % 2 != 1:
        message += " "
    bar_length = int(columns)
    spacer_bar = "=" * bar_length
    blanks = bar_length - 2
    space_length = int((blanks - len(message)) / 2)
    message_bar = "=" + " " * space_length + message + " " * space_length + "="
    print(spacer_bar)

    print(message_bar)
    print(spacer_bar)
    return True


def print_progress(current, total):
    """Prints the percentage of a progress

    It works for FOR loops, requires to know the full length of the loop.
    Prints to the standard output.

    Args:
        current (int): Current value in the range of the loop.
        total (int): The length of the loop.

    Returns:
        Nothing

    """
    if current == total:
        sys.stdout.write("Progress {:.2%}\n".format(1.0 * current / total))
    else:
        sys.stdout.write("\rProgress {:.2%}".format(1.0 * current / total))
    sys.stdout.flush()
    return


def get_twilight_time(date_obs):
    """Get end/start time of evening/morning twilight

    Notes:
        Taken from David Sanmartim's development

    Args:
        date_obs (list): List of all the dates from data.

    Returns:
        twilight_evening (str): Evening twilight time in the format 'YYYY-MM-DDTHH:MM:SS.SS'
        twilight_morning (str): Morning twilight time in the format 'YYYY-MM-DDTHH:MM:SS.SS'
        sun_set_time (str): Sun set time in the format 'YYYY-MM-DDTHH:MM:SS.SS'
        sun_rise_time (str): Sun rise time in the format 'YYYY-MM-DDTHH:MM:SS.SS'

    """
    # observatory(str): Observatory name.
    observatory = 'SOAR Telescope'
    geodetic_location = ['-70d44m01.11s', '-30d14m16.41s', 2748]
    # longitude (str): Geographic longitude in string format
    longitude = geodetic_location[0]
    # latitude (str): Geographic latitude in string format.
    latitude = geodetic_location[1]
    # elevation (int): Geographic elevation in meters above sea level
    elevation = geodetic_location[2]
    # timezone (str): Time zone.
    timezone = 'UTC'
    # description(str): Observatory description
    description = 'Soar Telescope on Cerro Pachon, Chile'

    soar_loc = EarthLocation.from_geodetic(longitude,
                                           latitude,
                                           elevation * u.m,
                                           ellipsoid='WGS84')

    soar = Observer(name=observatory,
                    location=soar_loc,
                    timezone=timezone,
                    description=description)

    time_first_frame, time_last_frame = Time(min(date_obs)), Time(max(date_obs))

    twilight_evening = soar.twilight_evening_astronomical(Time(time_first_frame),
                                                          which='nearest').isot
    twilight_morning = soar.twilight_morning_astronomical(Time(time_last_frame),
                                                          which='nearest').isot
    sun_set_time = soar.sun_set_time(Time(time_first_frame),
                                     which='nearest').isot
    sun_rise_time = soar.sun_rise_time(Time(time_last_frame),
                                       which='nearest').isot
    log.debug('Sun Set ' + sun_set_time)
    log.debug('Sun Rise ' + sun_rise_time)
    return twilight_evening, twilight_morning, sun_set_time, sun_rise_time


def image_overscan(ccd, overscan_region, add_keyword=False):
    """Apply overscan to data

    Args:
        ccd (object): A ccdproc.CCDData instance
        overscan_region (str): The overscan region in the format '[x1:x2,y1:y2]'
        where x is the spectral axis and y is the spatial axis.
        add_keyword (bool): Tells ccdproc whether to add a keyword or not.
        Default False.

    Returns:
        ccd (object): Overscan corrected ccdproc.CCDData instance

    """
    ccd = ccdproc.subtract_overscan(ccd=ccd,
                                    median=True,
                                    fits_section=overscan_region,
                                    add_keyword=add_keyword)
    ccd.header.add_history('Applied overscan correction ' + overscan_region)
    return ccd


def image_trim(ccd, trim_section, add_keyword=False):
    """Trim image to a given section

    Args:
        ccd (object): A ccdproc.CCDData instance
        trim_section (str): The trimming section in the format '[x1:x2,y1:y2]'
        where x is the spectral axis and y is the spatial axis.
        add_keyword (bool): Tells ccdproc whether to add a keyword or not.
        Default False.

    Returns:
        ccd (object): Trimmed ccdproc.CCDData instance

    """
    ccd = ccdproc.trim_image(ccd=ccd,
                             fits_section=trim_section,
                             add_keyword=add_keyword)
    ccd.header.add_history('Trimmed image to ' + trim_section)
    return ccd


def get_slit_trim_section(master_flat):
    """Find the slit edges to trim all data

    Using a master flat, ideally good signal to noise ratio, this function will
    identify the edges of the slit projected into the detector. Having this data
    will allow to reduce the overall processing time and also reduce the
    introduction of artifacts due to non-illuminated regions in the detectors,
    such as NaNs -INF +INF, etc.

    Args:
        master_flat (object): A ccdproc.CCDData instance

    Returns:
        slit_trim_section (str): Trim section in spatial direction in the format
        [:,slit_lower_limit:slit_higher_limit]

    """
    x, y = master_flat.data.shape
    # Using the middle point to make calculations, usually flats have good
    # illumination already at the middle.
    middle = int(y / 2.)
    ccd_section = master_flat.data[:, middle:middle + 200]
    ccd_section_median = np.median(ccd_section, axis=1)
    # spatial_axis = range(len(ccd_section_median))

    # Spatial half will be used later to constrain the detection of the first
    # edge before the first half.
    spatial_half = len(ccd_section_median) / 2.
    # pseudo-derivative to help finding the edges.
    pseudo_derivative = np.array(
        [abs(ccd_section_median[i + 1] - ccd_section_median[i]) for i in range(0, len(ccd_section_median) - 1)])
    filtered_data = np.where(np.abs(pseudo_derivative > 0.5 * pseudo_derivative.max()),
                             pseudo_derivative,
                             None)
    peaks = signal.argrelmax(filtered_data, axis=0, order=3)[0]
    # print(peaks)

    slit_trim_section = None
    if len(peaks) > 2 or peaks == []:
        log.debug('No trim section')
    else:
        # print(peaks, flat_files.grating[flat_files.file == file], flat_files.slit[flat_files.file == file])
        if len(peaks) == 2:
            # This is the ideal case, when the two edges of the slit are found.
            low, high = peaks
            slit_trim_section = '[:,{:d}:{:d}]'.format(low, high)
        elif len(peaks) == 1:
            # when only one peak is found it will choose the largest region from the spatial axis center to one edge.
            if peaks[0] <= spatial_half:
                slit_trim_section = '[:,{:d}:{:d}]'.format(peaks[0], len(ccd_section_median))
            else:
                slit_trim_section = '[:,{:d}:{:d}]'.format(0, peaks[0])
    return slit_trim_section


def cosmicray_rejection(ccd, mask_only=False):
    """Do cosmic ray rejection

      Notes:
          OBS: cosmic ray rejection is working pretty well by defining gain = 1. It's not working when we use the real
          gain of the image. In this case the sky level changes by a factor equal the gain.
          Function to determine the sigfrac and objlim: y = 0.16 * exptime + 1.2

      Args:
          ccd (object): CCDData Object
          mask_only (bool): In some cases you may want to obtain the cosmic
          rays mask only

      Returns:
          ccd (object): The modified CCDData object

      """
    # TODO (simon): Validate this method
    if ccd.header['OBSTYPE'] == 'OBJECT':
        value = 0.16 * float(ccd.header['EXPTIME']) + 1.2
        log.info('Cleaning cosmic rays... ')
        # ccd.data= ccdproc.cosmicray_lacosmic(ccd, sigclip=2.5, sigfrac=value, objlim=value,
        ccd.data, mask = ccdproc.cosmicray_lacosmic(ccd.data, sigclip=2.5, sigfrac=value, objlim=value,
                                                     gain=float(ccd.header['GAIN']),
                                                     readnoise=float(ccd.header['RDNOISE']),
                                                     satlevel=np.inf, sepmed=True, fsmode='median',
                                                     psfmodel='gaussy', verbose=False)
        ccd.header.add_history("Cosmic rays rejected with LACosmic")
        log.info("Cosmic rays rejected with LACosmic")
        if mask_only:
            return mask
        else:
            return ccd
    else:
        log.info('Skipping cosmic ray rejection for image of datatype: {:s}'.format(ccd.header['OBSTYPE']))
        return ccd



def get_best_flat(flat_name):
    """Look for matching master flat

    Given a basename for masterflats this function will find the name of the files that matches the basename and then
    will choose the first. Ideally this should go further as to check signal, time gap, etc.
    After it identifies the file it will load it using ccdproc.CCDData and return it along the filename.
    In case if fails it will return None instead of master_flat and another None instead of master_flat_name.

    Args:
        flat_name (str): Full path of masterflat basename. Ends in '*.fits'.

    Returns:
        master_flat (object): A ccdproc.CCDData instance
        master_flat_name (str): Full path to the chosen masterflat.

    """
    flat_list = glob.glob(flat_name)
    print(flat_name)
    print(flat_list)
    if len(flat_list) > 0:
        if len(flat_list) == 1:
            master_flat_name = flat_list[0]
        else:
            master_flat_name = flat_list[0]
        # elif any('dome' in flat for flat in flat_list):
        #     master_flat_name =

        master_flat = CCDData.read(master_flat_name, unit=u.adu)
        log.info('Found suitable master flat: ' + master_flat_name)
        return master_flat, master_flat_name
    else:
        log.error('There is no flat available')
        return None, None


def print_default_args(args):
    """Print default values of arguments.

    This is mostly helpful for debug but people not familiar with the software
    might find it useful as well

    Args:
        args (object): An argparse instance

    """
    arg_name = {'auto_clean': '--auto-clean',
                'clean_cosmic': '-c, --cosmic',
                'debug_mode': '--debug',
                'flat_normalize': '--flat-normalize',
                'ignore_bias': '--ignore-bias',
                'log_to_file': '--log-to-file',
                'norm_order': '--flat-norm-order',
                'raw_path': '--raw-path',
                'red_path': '--red-path',
                'saturation_limit': '--saturation',
                'destiny': '-d --proc-path',
                'interactive_ws': '-i --non-interactive',
                'lamp_all_night': '-r --reference-lamp',
                'lamp_file': '-l --lamp-file',
                'output_prefix': '-o --output-prefix',
                'pattern': '-s --search-pattern',
                'procmode': '-m --proc-mode',
                'reference_dir': '-R --reference-files',
                'source': '-p --data-path'}
    for key in args.__dict__:
        log.info('Default value for {:s} is {:s}'.format(arg_name[key],
                                       str(args.__getattribute__(key))))