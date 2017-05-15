from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import sys
import time
import glob
import random
import logging
from ccdproc import CCDData
from astropy.coordinates import EarthLocation
from astropy.time import Time, TimeDelta
from astroplan import Observer
from astropy import units as u

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

    There are some cases when the raw data comes with duplicated keywords. The origin has not been tracked down.
    The solution is to identify the duplicated keywords and the remove all but one from the end backwards.

    Args:
        night_dir (str): The full path for the raw data location

    """

    files = glob.glob(os.path.join(night_dir, '*.fits'))
    random_file = random.choice(files)
    ccd = CCDData.read(random_file)
    random_header = ccd.header
    multiple_keys = []
    for keyword in random_header.keys():
        if keyword != '':
            if random_header.count(keyword) > 1:
                if keyword not in multiple_keys:
                    multiple_keys.append(keyword)
    if multiple_keys != []:
        for image_file in files:
            try:
                print(image_file)
                ccd = CCDData.read(image_file)
                for keyword in multiple_keys:
                    while ccd.header.count(keyword) > 1:
                        ccd.header.remove(keyword, ccd.header.count(keyword) - 1)
                log.warning('Overwriting file with duplicated keywords removed')
                log.info('File %s overwritten', image_file)
                ccd.write(image_file, clobber=True)
            except IOError as error:
                log.error(error)


def ra_dec_to_deg(right_ascension, declination):
    """Converts right ascension and declination to degrees

    Args:
        right_ascension (str): right ascension
        declination (str): declination

    Returns:
        right ascension and declination in degrees

    """
    right_ascension = right_ascension.split(":")
    declination = declination.split(":")
    # RIGHT ASCENTION conversion
    right_ascension_deg = (float(right_ascension[0])
                           + (float(right_ascension[1])
                              + (float(right_ascension[2]) / 60.)) / 60.) * (360. / 24.)
    # DECLINATION conversion
    # print(declination)
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

    Prints a spacer of 80 columns with  and 3 rows height. The first and last rows contains the symbol "="
    repeated 80 times. The middle row contains the message centered and the extremes has one single "=" symbol.
    The only functionality of this is aesthetic.

    Args:
        message (str): a message to be printed

    Returns:
        a True boolean

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

    It works for FOR loops, requires to know the full length of the loop. Prints to the standard output.

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

    soar_loc = EarthLocation.from_geodetic(longitude, latitude, elevation * u.m, ellipsoid='WGS84')

    soar = Observer(name=observatory, location=soar_loc, timezone=timezone, description=description)

    time_first_frame, time_last_frame = Time(min(date_obs)), Time(max(date_obs))

    twilight_evening = soar.twilight_evening_astronomical(Time(time_first_frame), which='nearest').isot
    twilight_morning = soar.twilight_morning_astronomical(Time(time_last_frame), which='nearest').isot
    sun_set_time = soar.sun_set_time(Time(time_first_frame), which='nearest').isot
    sun_rise_time = soar.sun_rise_time(Time(time_last_frame), which='nearest').isot
    log.debug('Sun Set ' + sun_set_time)
    log.debug('Sun Rise ' + sun_rise_time)
    return twilight_evening, twilight_morning, sun_set_time, sun_rise_time
