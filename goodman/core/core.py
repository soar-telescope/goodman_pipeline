from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import re
import sys
import time
import glob
import random
import logging
import ccdproc
import numpy as np
import numpy.ma as ma
import matplotlib
import pandas
import shutil
import subprocess
import goodman

from threading import Timer

# matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
from ccdproc import CCDData, ImageFileCollection
from astropy.coordinates import EarthLocation
from astropy.time import Time, TimeDelta
from astropy.stats import sigma_clip
from astroplan import Observer
from astropy import units as u
from astropy.io import fits
from astropy.modeling import (models, fitting, Model)
from scipy import signal

log_ccd = logging.getLogger('goodmanccd.core')
log_spec = logging.getLogger('redspec.core')


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
    log_ccd.debug('Finding duplicated keywords')
    log_ccd.warning('Files headers will be overwritten')
    files = glob.glob(os.path.join(night_dir, '*.fits'))
    # Pick a random file to find duplicated keywords
    random_file = random.choice(files)
    ccd = CCDData.read(random_file, unit=u.adu)
    header = ccd.header
    # Put the duplicated keywords in a list
    multiple_keys = []
    for keyword in header.keys():
        if keyword != '':
            if header.count(keyword) > 1:
                if keyword not in multiple_keys:
                    multiple_keys.append(keyword)
    if multiple_keys != []:
        log_ccd.debug('Found {:d} duplicated keyword '
                 '{:s}'.format(len(multiple_keys),
                               's' if len(multiple_keys) > 1 else ''))

        for image_file in files:
            log_ccd.debug('Processing Image File: {:s}'.format(image_file))
            try:
                ccd = CCDData.read(image_file, unit=u.adu)
                for keyword in multiple_keys:
                    while ccd.header.count(keyword) > 1:
                        ccd.header.remove(keyword,
                                          ccd.header.count(keyword) - 1)
                log_ccd.warning('Overwriting file with duplicated keywords removed')
                log_ccd.debug('File %s overwritten', image_file)
                ccd.write(image_file, clobber=True)
            except IOError as error:
                log_ccd.error(error)


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
                              + (float(right_ascension[2]) / 60.)) / 60.) * \
                          (360. / 24.)

    # DECLINATION conversion
    if float(declination[0]) == abs(float(declination[0])):
        sign = 1
    else:
        sign = -1
    declination_deg = sign * (abs(float(declination[0]))
                              + (float(declination[1])
                                 + (float(declination[2]) / 60.)) / 60.)
    return right_ascension_deg, declination_deg


def print_spacers(message):
    """Miscellaneous function to print uniform spacers

    Prints a spacer of 80 columns with  and 3 rows height. The first and last
    rows contains the symbol "=" repeated 80 times. The middle row contains the
    message centered and the extremes has one single "=" symbol.
    The only functionality of this is aesthetic.

    Args:
        message (str): A message to be printed

    Returns:
        True (bool): A True value

    """
    # define the width of the message
    columns = 80
    if len(message) % 2 == 1 and int(columns) % 2 != 1:
        message += " "

    bar_length = int(columns)

    # compose bars top and bottom
    spacer_bar = "=" * bar_length

    blanks = bar_length - 2
    space_length = int((blanks - len(message)) / 2)

    # compose the message
    message_bar = "=" + " " * space_length + message + " " * space_length + "="

    print(spacer_bar)
    print(message_bar)
    print(spacer_bar)
    return True


def print_progress(current, total):
    """Prints the percentage of a progress

    It works for FOR loops, requires to know the full length of the loop.
    Prints to the standard output.

    Notes:
        A possible improvement for this is to run it using multithreading

    Args:
        current (int): Current value in the range of the loop.
        total (int): The length of the loop.

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
        twilight_evening (str): Evening twilight time in the format
            'YYYY-MM-DDTHH:MM:SS.SS'
        twilight_morning (str): Morning twilight time in the format
            'YYYY-MM-DDTHH:MM:SS.SS'
        sun_set_time (str): Sun set time in the format 'YYYY-MM-DDTHH:MM:SS.SS'
        sun_rise_time (str): Sun rise time in the format
            'YYYY-MM-DDTHH:MM:SS.SS'

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

    twilight_evening = soar.twilight_evening_astronomical(
        Time(time_first_frame), which='nearest').isot

    twilight_morning = soar.twilight_morning_astronomical(
        Time(time_last_frame), which='nearest').isot

    sun_set_time = soar.sun_set_time(
        Time(time_first_frame), which='nearest').isot

    sun_rise_time = soar.sun_rise_time(
        Time(time_last_frame), which='nearest').isot

    log_ccd.debug('Sun Set ' + sun_set_time)
    log_ccd.debug('Sun Rise ' + sun_rise_time)

    return twilight_evening, twilight_morning, sun_set_time, sun_rise_time


def read_fits(full_path, technique='Unknown'):
    assert os.path.isfile(full_path)
    ccd = CCDData.read(full_path, unit=u.adu)

    ccd.header.set('GSP_VERS',
                   value=goodman.__version__,
                   comment='Goodman Spectroscopic Pipeline Version')

    ccd.header.set('GSP_FNAM',
                   value=os.path.basename(full_path),
                   comment='Original file name')

    ccd.header.set('GSP_PATH',
                   value=os.path.dirname(full_path),
                   comment='Location at moment of reduce')

    ccd.header.set('GSP_TECH',
                   value=technique,
                   comment='Observing technique')

    ccd.header.set('GSP_DATE',
                   value=time.strftime("%Y-%m-%d"),
                   comment='Processing date')

    ccd.header.set('GSP_OVER',
                   value='none',
                   comment='Overscan Region')

    ccd.header.set('GSP_TRIM',
                   value='none',
                   comment='Trim section')

    ccd.header.set('GSP_SLIT',
                   value='none',
                   comment='Slit trim section, slit illuminated area only')

    ccd.header.set('GSP_BIAS',
                   value='none',
                   comment='Master bias image')

    ccd.header.set('GSP_FLAT',
                   value='none',
                   comment='Master flat image')

    ccd.header.set('GSP_NORM',
                   value='none',
                   comment='Flat normalization method')

    ccd.header.set('GSP_COSM',
                   value='none',
                   comment='Cosmic ray rejection method')

    ccd.header.set('GSP_WRMS',
                   value='none',
                   comment='Wavelength solution RMS Error')

    ccd.header.set('GSP_WPOI',
                   value='none',
                   comment='Number of points used to '
                           'calculate wavelength solution')

    ccd.header.set('GSP_WREJ',
                   value='none',
                   comment='Number of points rejected')

    ccd.header.add_blank('-- Goodman Spectroscopic Pipeline --', before='GSP_VERS')
    # ccd.header.set('', value='', comment='')
    # ccd.header.set('', value='', comment='')
    # ccd.header.set('', value='', comment='')
    # ccd.header.set('', value='', comment='')
    # ccd.header.set('', value='', comment='')
    # ccd.header.set('', value='', comment='')
    return ccd

def image_overscan(ccd, overscan_region, add_keyword=False):
    """Apply overscan to data

    Uses ccdproc.subtract_overscan to perform the task.

    Notes:
        The overscan_region argument uses FITS convention, just like IRAF,
        therefore is 1 based. i.e. it starts in 1 not 0.

    Args:
        ccd (object): A ccdproc.CCDData instance
        overscan_region (str): The overscan region in the format '[x1:x2,y1:y2]'
            where x is the spectral axis and y is the spatial axis.
        add_keyword (bool): Tells ccdproc whether to add a keyword or not.
            Default False.

    Returns:
        ccd (object): Overscan corrected ccdproc.CCDData instance

    """
    log_ccd.debug('Applying overscan Correction: {:s}'.format(overscan_region))
    ccd = ccdproc.subtract_overscan(ccd=ccd,
                                    median=True,
                                    fits_section=overscan_region,
                                    add_keyword=add_keyword)

    ccd.header['GSP_OVER'] = (overscan_region, 'Overscan region')
    return ccd


def image_trim(ccd, trim_section, trim_type, add_keyword=False):
    """Trim image to a given section

    Notes:
        The overscan_region argument uses FITS convention, just like IRAF,
        therefore is 1 based. i.e. it starts in 1 not 0.

    Args:
        ccd (object): A ccdproc.CCDData instance
        trim_section (str): The trimming section in the format '[x1:x2,y1:y2]'
            where x is the spectral axis and y is the spatial axis.
        trim_type (str): default or slit trim
        add_keyword (bool): Tells ccdproc whether to add a keyword or not.
            Default False.

    Returns:
        ccd (object): Trimmed ccdproc.CCDData instance

    """
    ccd = ccdproc.trim_image(ccd=ccd,
                             fits_section=trim_section,
                             add_keyword=add_keyword)
    if trim_type == 'trimsec':
        ccd.header['GSP_TRIM'] = (trim_section, 'Trimsection from TRIMSEC')
    elif trim_type == 'slit':
        ccd.header['GSP_SLIT'] = (trim_section,
                                  'Slit trim section, slit illuminated '
                                  'area only.')
    else:
        log_ccd.warning('Unrecognized trim type')

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
    spatial_axis = range(len(ccd_section_median))

    # set values for initial box model definition
    box_max = np.max(ccd_section_median)
    box_center = len(ccd_section_median) / 2.
    box_width = .75 * len(ccd_section_median)

    # box model definition
    box_model = models.Box1D(amplitude=box_max, x_0=box_center, width=box_width)

    box_fitter = fitting.SimplexLSQFitter()

    fitted_box = box_fitter(box_model, spatial_axis, ccd_section_median)

    # the number of pixels that will be removed from the detected edge of the
    # image on each side
    offset = 10

    # this defines a preliminary set of slit limit
    l_lim = fitted_box.x_0.value - fitted_box.width.value / 2. + offset

    h_lim = fitted_box.x_0.value + fitted_box.width.value / 2. - offset

    # Here we force the slit limits within the boundaries of the data (image)
    low_lim = int(np.max([0 + offset, l_lim]))

    high_lim = int(np.min([h_lim, len(ccd_section_median) - offset]))

    # define the slit trim section as (IRA
    slit_trim_section = '[:,{:d}:{:d}]'.format(low_lim, high_lim)

    # debugging plots that have to be manually turned on
    if False:
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
        plt.title('Slit Edge Detection')
        plt.plot(box_model(spatial_axis), color='c', label='Initial Box1D')
        plt.plot(fitted_box(spatial_axis), color='k', label='Fitted Box1D')
        plt.plot(ccd_section_median, label='Median Along Disp.')
        # plt.plot(pseudo_derivative, color='g', label='Pseudo Derivative')
        plt.axvline(None, color='r', label='Detected Edges')
        plt.axvline(low_lim, color='r')
        plt.axvline(high_lim, color='r')
        # for peak in peaks:
        #     plt.axvline(peak, color='r')
        plt.legend(loc='best')
        plt.show()

    # plt.imshow(master_flat.data[low_lim:high_lim, :])
    # plt.axvline(low_lim, color='r')
    # plt.axvline(high_lim, color='r')
    # plt.show()

    return slit_trim_section


def dcr_cosmicray_rejection(data_path, in_file, prefix, dcr_par_dir,
                            delete=False):
    """Runs an external code for cosmic ray rejection

    DCR was created by Wojtek Pych and the code can be obtained from
    http://users.camk.edu.pl/pych/DCR/ and is written in C. Contrary to
    ccdproc's LACosmic it actually applies the correction, and also doesn't
    update the mask attribute since it doesn't work with CCDData instances.

    The binary takes three positional arguments, they are: 1. input image,
    2. output image and 3. cosmic rays images. Also it needs that a dcr.par file
    is located in the directory. All this is implemented in this function if
    delete is True it will remove the original image and the cosmic rays image.
    The removal of the original image is absolutely safe when used in the
    context of the goodman pipeline, however if you want to implement it
    somewhere else, be careful.

    Notes:
        This function operates an external code therefore it doesn't return
        anything, instead it creates a new image.

    Args:
        data_path (str): Data location
        in_file (str): Name of the file to have its cosmic rays removed
        prefix (str): Prefix to add to the file with the cosmic rays removed
        dcr_par_dir (str): Directory of default dcr.par file
        delete (bool): True for deleting the input and cosmic ray file.

    """

    log_ccd.info('Removing cosmic rays using DCR by Wojtek Pych')
    log_ccd.debug('See http://users.camk.edu.pl/pych/DCR/')

    # add the prefix for the output file
    out_file = prefix + in_file

    # define the name for the cosmic rays file
    cosmic_file = 'cosmic_' + '_'.join(in_file.split('_')[1:])

    # define full path for all the files involved
    full_path_in = os.path.join(data_path, in_file)
    full_path_out = os.path.join(data_path, out_file)
    full_path_cosmic = os.path.join(data_path, cosmic_file)

    # this is the command for running dcr, all arguments are required
    command = 'dcr {:s} {:s} {:s}'.format(full_path_in,
                                          full_path_out,
                                          full_path_cosmic)

    log_ccd.debug('DCR command:')
    log_ccd.debug(command)
    # print(command.split(' '))

    # get the current working directory to go back to it later in case the
    # the pipeline has not been called from the same data directory.
    cwd = os.getcwd()

    # move to the directory were the data is, dcr is expecting a file dcr.par
    os.chdir(data_path)

    # check if file dcr.par exists
    while not os.path.isfile('dcr.par'):

        log_spec.debug('File dcr.par does not exist. Copying default one.')
        dcr_par_path = os.path.join(dcr_par_dir, 'dcr.par')
        log_ccd.debug('dcr.par full path: {:s}'.format(dcr_par_path))
        if os.path.isfile(dcr_par_path):
            shutil.copy2(dcr_par_path, data_path)
        else:
            log_ccd.error('Could not find dcr.par file')
    else:
        log_ccd.debug('File dcr.par exists.')

    # call dcr
    try:

        dcr = subprocess.Popen(command.split(),
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)

    except OSError as error:
        log_ccd.error(error)
        sys.exit('Your system can not locate the executable file dcr, try '
                 'moving it to /bin or create a symbolic link\n\n\tcd /bin\n\t'
                 'sudo ln -s /full/path/to/dcr')

        # return False

    # if the process is taking too long to respond, kill it
    # kill_process = lambda process: process.kill()
    def kill_process(process): process.kill()

    dcr_timer = Timer(5, kill_process, [dcr])
    try:
        dcr_timer.start()
        stdout, stderr = dcr.communicate()
    finally:
        dcr_timer.cancel()

    # wait for dcr to terminate
    # dcr.wait()

    # go back to the original directory. Could be the same.
    os.chdir(cwd)

    # If no error stderr is an empty string
    if stderr != b'':
        log_ccd.error(stderr)
        if b'dcr: not found' in stderr:
            sys.exit('Your system can not locate the executable file dcr, try '
                     'moving it to /bin or create a symbolic link\n\n\tcd '
                     '/bin\n\tsudo ln -s /full/path/to/dcr')
    else:
        for output_line in stdout.split(b'\n'):
            log_ccd.debug(output_line)

    # delete extra files only if the execution ended without error
    if delete and stderr == b'' and b'USAGE:' not in stdout:
        try:
            log_ccd.warning('Removing input file: {:s}'.format(full_path_in))
            os.unlink(full_path_in)
        except OSError as error:
            log_ccd.error(error)

        try:
            log_ccd.warning(
                'Removing cosmic rays file: {:s}'.format(full_path_cosmic))
            os.unlink(full_path_cosmic)
        except OSError as error:
            log_ccd.error(error)


def lacosmic_cosmicray_rejection(ccd, mask_only=False):
    """Do cosmic ray rejection using ccdproc.LACosmic

    This function in fact does not apply any correction, it detects the cosmic
    rays and updates the attribute mask of the ccd object (CCDData instance).
    The attribute mask is used later as a mask for the pixels hit by cosmic rays

      Notes:
          OBS: cosmic ray rejection is working pretty well by defining gain = 1.
          It's not working when we use the real gain of the image. In this case
          the sky level changes by a factor equal the gain.
          Function to determine the sigfrac and objlim: y = 0.16 * exptime + 1.2

      Args:
          ccd (object): A ccdproc.CCDData instance.
          mask_only (bool): In some cases you may want to obtain the cosmic
              rays mask only.

      Returns:
          ccd (object): A CCDData instance with the mask attribute updated or an
          array which corresponds to the mask.

    """
    if ccd.header['OBSTYPE'] == 'OBJECT':
        value = 0.16 * float(ccd.header['EXPTIME']) + 1.2
        log_ccd.info('Cleaning cosmic rays... ')

        ccd = ccdproc.cosmicray_lacosmic(
            ccd,
            sigclip=2.5,
            sigfrac=value,
            objlim=value,
            gain=float(ccd.header['GAIN']),
            readnoise=float(ccd.header['RDNOISE']),
            satlevel=np.inf,
            sepmed=True,
            fsmode='median',
            psfmodel='gaussy',
            verbose=False)

        ccd.header['GSP_COSM'] = ("LACosmic", "Cosmic ray rejection method")
        log_ccd.info("Cosmic rays rejected with LACosmic")
        if mask_only:
            return ccd.mask
        else:
            return ccd
    else:
        log_ccd.debug('Skipping cosmic ray rejection for image of datatype: '
                 '{:s}'.format(ccd.header['OBSTYPE']))
        return ccd


def call_cosmic_rejection(ccd, image_name, out_prefix, red_path,
                          dcr_par, keep_files=False, prefix='c', method='dcr'):
    """Call for the appropriate cosmic ray rejection method

    There are three options when dealing with cosmic ray rejection in this
    pipeline, the first is ``dcr`` which is a program written in C by Wojtek
    Pych (http://users.camk.edu.pl/pych/DCR/) and works very well for
    spectroscopy the only negative aspect is that integration with python was
    difficult and not native.

    Args:
        ccd (object): a ccdproc.CCDData instance.
        image_name (str): Science image name.
        out_prefix (str): Partial prefix to be added to the image name. Related
            to previous processes and not cosmic ray rejection.
        red_path (str): Path to reduced data directory.
        dcr_par (str): Path to dcr.par file.
        keep_files (bool): If True, the original file and the cosmic ray mask
            will not be deleted. Default is False.
        prefix (str): Cosmic ray rejection related prefix to be added to image
            name.
        method (str): Method to use for cosmic ray rejection. There are three
            options: dcr, lacosmic and none.

    """

    if method == 'dcr':
        log_ccd.warning('DCR does apply the correction to images if you want '
                        'the mask use --keep-cosmic-files')
        full_path = os.path.join(red_path, out_prefix + image_name)

        ccd.header['GSP_COSM'] = ("DCR", "Cosmic ray rejection method")

        ccd.write(full_path, clobber=True)
        log_ccd.info('Saving image: {:s}'.format(full_path))

        in_file = out_prefix + image_name

        dcr_cosmicray_rejection(data_path=red_path,
                                in_file=in_file,
                                prefix=prefix,
                                dcr_par_dir=dcr_par,
                                delete=keep_files)

    elif method == 'lacosmic':
        log_ccd.warning('LACosmic does not apply the correction to images '
                        'instead it updates the mask attribute for CCDData '
                        'objects. For saved files the mask is a fits extension')

        ccd = lacosmic_cosmicray_rejection(ccd=ccd)

        out_prefix = prefix + out_prefix
        full_path = os.path.join(red_path, out_prefix + image_name)

        ccd.write(full_path, clobber=True)
        log_ccd.info('Saving image: {:s}'.format(full_path))

    elif method == 'none':
        full_path = os.path.join(red_path, out_prefix + image_name)
        log_ccd.warning("--cosmic set to 'none'")
        ccd.write(full_path, clobber=True)
        log_ccd.info('Saving image: {:s}'.format(full_path))

    else:
        log_ccd.error('Unrecognized Cosmic Method {:s}'.format(method))


def get_best_flat(flat_name):
    """Look for matching master flat

    Given a basename for masterflats defined as a combination of key parameters
    extracted from the header of the image that we want to flatfield, this
    function will find the name of the files that matches the base name and then
    will choose the first. Ideally this should go further as to check signal,
    time gap, etc.
    After it identifies the file it will load it using ccdproc.CCDData and
    return it along the filename.
    In case if fails it will return None instead of master_flat and another
    None instead of master_flat_name.

    Args:
        flat_name (str): Full path of masterflat basename. Ends in '*.fits' for
            globbing.

    Returns:
        master_flat (object): A ccdproc.CCDData instance
        master_flat_name (str): Full path to the chosen masterflat.

    """
    flat_list = glob.glob(flat_name)
    log_ccd.debug('Flat base name {:s}'.format(flat_name))
    log_ccd.debug('Matching master flats found: {:d}'.format(len(flat_list)))
    if len(flat_list) > 0:
        if len(flat_list) == 1:
            master_flat_name = flat_list[0]
        else:
            master_flat_name = flat_list[0]
        # elif any('dome' in flat for flat in flat_list):
        #     master_flat_name =

        master_flat = CCDData.read(master_flat_name, unit=u.adu)
        log_ccd.debug('Found suitable master flat: {:s}'.format(master_flat_name))
        return master_flat, master_flat_name
    else:
        log_ccd.error('There is no flat available')
        return None, None


def print_default_args(args):
    """Print default values of arguments.

    This is mostly helpful for debug but people not familiar with the software
    might find it useful as well

    Notes:
        This function is deprecated.

    Notes:
        This is not dynamically updated so use with caution

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
                'interactive_ws': '-i --interactive',
                'lamp_all_night': '-r --reference-lamp',
                'lamp_file': '-l --lamp-file',
                'output_prefix': '-o --output-prefix',
                'pattern': '-s --search-pattern',
                'procmode': '-m --proc-mode',
                'reference_dir': '-R --reference-files',
                'source': '-p --data-path',
                'save_plots': '--save-plots',
                'dcr_par_dir': '--dcr-par-dir'}
    for key in args.__dict__:
        log_ccd.debug('Default value for {:s} is {:s}'.format(
            arg_name[key],
            str(args.__getattribute__(key))))


def normalize_master_flat(master, name, method='simple', order=15):
    """ Master flat normalization method

    This function normalize a master flat in three possible ways:
     *mean*: simply divide the data by its mean

     *simple*: Calculates the median along the spatial axis in order to obtain
     the dispersion profile. Then fits a Chebyshev1D model and apply this to all
     the data.

     *full*: This is for experimental purposes only because it takes a lot of
     time to process. It will fit a model to each line along the dispersion axis
     and then divide it by the fitted model. I do not recommend this method
     unless you have a good reason as well as a powerful computer..

    Args:
        master (object): Master flat. Has to be a ccdproc.CCDData instance
        name (str): Full path of master flat prior to normalization
        method (str): Normalization method, 'mean', 'simple' or 'full'
        order (int): Order of the polinomial to be fitted.

    Returns:
        master (object):  The normalized master flat. ccdproc.CCDData instance

    """
    assert isinstance(master, CCDData)

    # define new name, base path and full new name
    new_name = 'norm_' + name.split('/')[-1]
    path = '/'.join(name.split('/')[0:-1])
    norm_name = os.path.join(path, new_name)

    if method == 'mean':
        log_ccd.debug('Normalizing by mean')
        master.data /= master.data.mean()

        master.header['GSP_NORM'] = ('mean', 'Flat normalization method')

    elif method == 'simple' or method == 'full':
        log_ccd.debug('Normalizing flat by {:s} model'.format(method))

        # Initialize Fitting models and fitter
        model_init = models.Chebyshev1D(degree=order)
        model_fitter = fitting.LevMarLSQFitter()

        # get data shape
        x_size, y_size = master.data.shape
        x_axis = range(y_size)

        if method == 'simple':
            # get profile along dispersion axis to fit a model to use for
            # normalization
            profile = np.median(master.data, axis=0)

            # do the actual fit
            fit = model_fitter(model_init, x_axis, profile)

            # convert fit into an array
            fit_array = fit(x_axis)

            # pythonic way to divide an array by a vector
            master.data = master.data / fit_array[None, :]

            # master.header.add_history('Flat Normalized by simple model')
            master.header['GSP_NORM'] = ('simple', 'Flat normalization method')

        elif method == 'full':
            log_ccd.warning('This part of the code was left here for '
                            'experimental purposes only')
            log_ccd.warning('This procedure takes a lot to process, you might '
                            'want to see other method such as "simple" or '
                            '"mean".')
            for i in range(x_size):
                fit = model_fitter(model_init, x_axis, master.data[i])
                master.data[i] = master.data[i] / fit(x_axis)
            # master.header.add_history('Flat Normalized by full model')
            master.header['GSP_NORM'] = ('full', 'Flat normalization method')

    # write normalized flat to a file
    master.write(norm_name, clobber=True)

    return master, norm_name


def get_central_wavelength(grating, grt_ang, cam_ang):
    """Calculates the central wavelength for a given spectroscopic mode

    The equation used to calculate the central wavelength is the following


    .. math::
        \\lambda_{central} = \\frac{1e6}{GRAT}
        \\sin\\left(\\frac{\\alpha \\pi}{180}\\right) +
        \\sin\\left(\\frac{\\beta \\pi}{180}\\right)


    Args:
        grating (str): Grating frequency as a string. Example '400'.
        grt_ang (str): Grating Angle as a string. Example '12.0'.
        cam_ang (str): Camera Angle as a string. Example '20.0'

    Returns:
        central_wavelength (float): Central wavelength as a float value.

    """

    grating_frequency = float(grating)
    alpha = float(grt_ang)
    beta = float(cam_ang) - float(grt_ang)

    central_wavelength = (1e6 / grating_frequency) * \
                         (np.sin(alpha * np.pi / 180.) +
                          np.sin(beta * np.pi / 180.))

    log_spec.debug('Found {:.3f} as central wavelength'.format(central_wavelength))

    return central_wavelength


def remove_conflictive_keywords(path, file_list):
    """Removes problematic keywords

    The blue camera has a set of keywords whose comments contain non-ascii
    characters, in particular the degree symbol. Those keywords are not
    needed in any stage of the data reduction therefore they are removed.
    The data will be overwritten with the keywords removed. The user will
    need to have backups of raw data.

    Notes:
        This function solves a problem with old data, new headers are compliant
        with the headers.

    Args:
        path (str): Path to the folder containing the files
        file_list (list): List of files to remove keywords

    """
    log_ccd.debug('Removing conflictive keywords in Blue Camera Headers')
    log_ccd.warning('Files will be overwritten')
    for blue_file in file_list:
        full_path = os.path.join(path, blue_file)
        log_ccd.debug('Processing file {:s}'.format(blue_file))
        try:
            data, header = fits.getdata(full_path,
                                        header=True,
                                        ignore_missing_end=True)

            keys_to_remove = ['PARAM0',
                              'PARAM61',
                              'PARAM62',
                              'PARAM63',
                              'NAXIS3']

            if data.ndim == 3:
                header['NAXIS'] = 2
                data = data[0]

                log_ccd.debug('Modified file to be 2D instead of 3D '
                              '(problematic)')

            for keyword in keys_to_remove:
                header.remove(keyword)

                log_ccd.debug('Removed conflictive keyword '
                              '{:s}'.format(keyword))

            log_ccd.debug('Updated headers')

            fits.writeto(full_path,
                         data,
                         header,
                         clobber=True)

        except KeyError as error:
            log_ccd.debug(error)


# spectroscopy specific functions

def classify_spectroscopic_data(path, search_pattern):
    """Classify data by grouping them as pandas.DataFrame instances

    This functions uses ImageFileCollection from ccdproc. First it creates a
    collection of information regarding the images located in *path* that match
    the pattern *search_pattern*
    The information obtained are all keywords listed in the list *keywords*
    The ImageFileCollection is translated into pandas.DataFrame and then is used
    much like an SQL database to select and filter values and in that way put
    them in groups that are pandas.DataFrame instances.


    Args:
        path (str): Path to data location
        search_pattern (str): Prefix to match files.

    Returns:
        data_container (object): Instance of NightDataContainer

    """

    search_path = os.path.join(path, search_pattern + '*.fits')

    file_list = glob.glob(search_path)

    if file_list == []:
        log_spec.error('No file found using search pattern '
                       '"{:s}"'.format(search_pattern))

        sys.exit('Please use the argument --search-pattern to define the '
                 'common prefix for the files to be processed.')

    data_container = NightDataContainer(path=path,
                                        instrument=str('Red'),
                                        technique=str('Spectroscopy'))

    keywords = ['date',
                'slit',
                'date-obs',
                'obstype',
                'object',
                'exptime',
                'obsra',
                'obsdec',
                'grating',
                'cam_targ',
                'grt_targ',
                'filter',
                'filter2',
                'gain',
                'rdnoise']

    ifc = ImageFileCollection(path, keywords=keywords, filenames=file_list)

    pifc = ifc.summary.to_pandas()

    pifc['radeg'] = ''
    pifc['decdeg'] = ''
    for i in pifc.index.tolist():
        radeg, decdeg = ra_dec_to_deg(pifc.obsra.iloc[i], pifc.obsdec.iloc[i])

        pifc.iloc[i, pifc.columns.get_loc('radeg')] = '{:.2f}'.format(radeg)

        pifc.iloc[i, pifc.columns.get_loc('decdeg')] = '{:.2f}'.format(decdeg)
        # now we can compare using degrees

    confs = pifc.groupby(['slit',
                          'radeg',
                          'decdeg',
                          'grating',
                          'cam_targ',
                          'grt_targ',
                          'filter',
                          'filter2',
                          'gain',
                          'rdnoise']).size().reset_index().rename(
        columns={0: 'count'})

    for i in confs.index:
        spec_group = pifc[((pifc['slit'] == confs.iloc[i]['slit']) &
                           (pifc['radeg'] == confs.iloc[i]['radeg']) &
                           (pifc['decdeg'] == confs.iloc[i]['decdeg']) &
                           (pifc['grating'] == confs.iloc[i]['grating']) &
                           (pifc['cam_targ'] == confs.iloc[i]['cam_targ']) &
                           (pifc['grt_targ'] == confs.iloc[i]['grt_targ']) &
                           (pifc['filter'] == confs.iloc[i]['filter']) &
                           (pifc['filter2'] == confs.iloc[i]['filter2']) &
                           (pifc['gain'] == confs.iloc[i]['gain']) &
                           (pifc['rdnoise'] == confs.iloc[i]['rdnoise']))]

        group_obstype = spec_group.obstype.unique()

        if 'COMP' in group_obstype and len(group_obstype) == 1:
            log_spec.debug('Adding COMP group')
            data_container.add_comp_group(comp_group=spec_group)
        elif 'OBJECT' in group_obstype and len(group_obstype) == 1:
            log_spec.debug('Adding OBJECT group')
            data_container.add_object_group(object_group=spec_group)
        else:
            log_spec.debug('Adding OBJECT-COMP group')
            data_container.add_spec_group(spec_group=spec_group)

    return data_container


def search_comp_group(object_group, comp_groups):
    """Search for a suitable comparison lamp group

    In case a science target was observed without comparison lamps, usually
    right before or right after, this function will look for a compatible set
    obtained at a different time or pointing.

    Notes:
        This methodology is not recommended for radial velocity studies.

    Args:
        object_group (object): A pandas.DataFrame instances containing a group
            of images for a given scientific target.
        comp_groups (list): A list in which every element is a pandas.DataFrame
            that contains information regarding groups of comparison lamps.

    Returns:

    """
    log_spec.debug('Finding a suitable comparison lamp group')

    object_confs = object_group.groupby(['slit',
                                         'grating',
                                         'cam_targ',
                                         'grt_targ',
                                         'filter',
                                         'filter2']
                                        ).size().reset_index()
    # .rename(columns={0: 'count'})

    for comp_group in comp_groups:

        if ((comp_group['slit'] == object_confs.iloc[0]['slit']) &
                (comp_group['grating'] == object_confs.iloc[0]['grating']) &
                (comp_group['cam_targ'] == object_confs.iloc[0]['cam_targ']) &
                (comp_group['grt_targ'] == object_confs.iloc[0]['grt_targ']) &
                (comp_group['filter'] == object_confs.iloc[0]['filter']) &
                (comp_group['filter2'] == object_confs.iloc[0]['filter2'])).all():
            log_spec.debug('Found a matching comparison lamp group')

            return comp_group

    raise NoMatchFound


def spectroscopic_extraction(ccd, extraction,
                             comp_list=None,
                             nfind=3,
                             n_sigma_extract=10,
                             plots=False):
    """This function does not do the actual extraction but prepares the data

    There are several steps involved in a spectroscopic extraction, this
    function manages them.

    Args:
        ccd (object): A ccdproc.CCDData Instance
        extraction (str): Extraction type name. _simple_ or _optimal_
        comp_list (list): List of ccdproc.CCDData instances of COMP lamps data
        nfind (int): Maximum number of targets to be returned
        n_sigma_extract (int): Number of sigmas to be used for extraction
        plots (bool): If plots will be shown or not.

    Returns:
        extracted (list): List of ccdproc.CCDData instances
        comp_zones (list): List of ccdproc.CCDData instances

    Raises:
        NoTargetException (Exception): A NoTargetException if there is no target
           found.

    """

    assert isinstance(ccd, CCDData)

    comp_zones = []
    extracted = []

    if comp_list is None:
        comp_list = []
    # print(comp_list)

    iccd = remove_background_by_median(ccd=ccd)

    profile_model = identify_targets(ccd=iccd, nfind=nfind, plots=plots)
    del (iccd)

    if profile_model is None:
        log_spec.critical('Target identification FAILED!')
        raise NoTargetException
    else:
        background_image = create_background_image(ccd=ccd,
                                                   profile_model=profile_model,
                                                   nsigma=n_sigma_extract,
                                                   separation=5)

    if isinstance(profile_model, Model):
        traces = trace_targets(ccd=ccd, profile=profile_model, plots=plots)
        # extract(ccd=ccd,
        #         spatial_profile=profile_model,
        #         n_sigma_extract=10,
        #         sampling_step=5)
        if 'CompoundModel' in profile_model.__class__.name:
            log_spec.debug(profile_model.submodel_names)
            for m in range(len(profile_model.submodel_names)):
                submodel_name = profile_model.submodel_names[m]

                ntrace = traces[m]
                model = profile_model[submodel_name]

                zone = get_extraction_zone(
                    ccd=ccd,
                    extraction=extraction,
                    trace=ntrace,
                    trace_index=m,
                    model=profile_model[submodel_name],
                    n_sigma_extract=n_sigma_extract,
                    plots=plots)

                for comp in comp_list:
                    comp_zone = get_extraction_zone(ccd=comp,
                                                    extraction=extraction,
                                                    trace=ntrace,
                                                    trace_index=m,
                                                    zone=zone,
                                                    plots=plots)
                    # since a comparison lamp only needs only the relative line
                    # center in the dispersion direction, therefore the flux is
                    # not important we are only calculating the median along the
                    # spatial direction
                    comp_zone.data = np.median(comp_zone.data, axis=0)
                    comp_zones.append(comp_zone)

                background_level = get_background_value(
                    background_image=background_image,
                    zone=zone)

                extracted_ccd = extract(ccd=ccd,
                                        trace=ntrace,
                                        spatial_profile=model,
                                        extraction=extraction,
                                        zone=zone,
                                        background_level=background_level,
                                        sampling_step=10,
                                        plots=plots)
                extracted.append(extracted_ccd)

                # if plots:
                #     plt.imshow(nccd.data)
                #     plt.show()

        else:
            ntrace = traces[0]

            zone = get_extraction_zone(
                ccd=ccd,
                extraction=extraction,
                trace=traces[0],
                trace_index=0,
                model=profile_model,
                n_sigma_extract=n_sigma_extract,
                plots=plots)

            for comp in comp_list:
                comp_zone = get_extraction_zone(ccd=comp,
                                                extraction=extraction,
                                                trace=ntrace,
                                                trace_index=0,
                                                zone=zone,
                                                plots=plots)

                # since a comparison lamp only needs the relative line
                # center in the dispersion direction, therefore the flux is not
                # important we are only calculating the median along the spatial
                # direction
                comp_zone.data = np.median(comp_zone.data, axis=0)
                comp_zones.append(comp_zone)

            background_level = get_background_value(
                background_image=background_image,
                zone=zone)

            extracted_ccd = extract(ccd=ccd,
                                    trace=ntrace,
                                    spatial_profile=profile_model,
                                    extraction=extraction,
                                    zone=zone,
                                    background_level=background_level,
                                    sampling_step=10,
                                    plots=plots)

            extracted.append(extracted_ccd)

            # if plots:
            #     plt.imshow(nccd.data)
            #     plt.show()

        # print(extracted)
        # print(comp_zones)
        return extracted, comp_zones

    elif profile_model is None:
        log_spec.warning("Didn't receive identified targets "
                         "from {:s}".format(ccd.header['GSP_FNAM']))
        raise NoTargetException
    else:
        log_spec.error('Got wrong input')


def identify_targets(ccd, nfind=3, plots=False):
    """Identify spectroscopic targets in an image

    This function collapses the image along the dispersion direction using a
    median, This highlights the spatial features present in a 2D spectrum
    (image), Then does a sigma clip to remove any features in order to fit the
    background level and shape, the fit is a linear function. Once the
    background has been removed it will equal to zero all negative values. It
    will perform a new sigma clipping but this time to determinate the
    background amplitude. Finally it finds all the peaks above the background
    level and pick the n largets ones. n is defined by nfind.

    Args:
        ccd (object): a ccdproc.CCDData instance
        nfind (int): Maximum number of targets to be returned
        plots (bool): to show debugging plots

    Returns:
        profile_model (object): an astropy.modeling.Model instance, it could be
            a Gaussian1D or CompoundModel (several Gaussian1D). Each of them
            represent a point source spectrum found.

    """
    if isinstance(ccd, CCDData):
        slit_size = re.sub('[a-zA-Z"]', '', ccd.header['SLIT'])
        serial_binning = int(ccd.header['CCDSUM'].split()[0])
        # order will be used for finding the peaks later but also as an initial
        # estimate for stddev of gaussian
        order = int(round(float(slit_size) / (0.15 * serial_binning)))

        if plots:
            plt.title(ccd.header['GSP_FNAM'])
            plt.imshow(ccd.data, clim=(30, 250))
            plt.xlabel('Dispersion Axis (x)')
            plt.ylabel('Spatial Axis (y)')
            plt.show()

        median_profile = np.median(ccd.data, axis=1)

        # Fitting Background

        # Before doing the fitting we will do a sigma clipping in order to
        # remove any feature.
        clipped_profile = sigma_clip(median_profile, sigma=2, iters=5)

        linear_model = models.Linear1D(slope=0,
                                       intercept=np.median(median_profile))

        linear_fitter = fitting.LinearLSQFitter()

        # the fitters do not support masked arrays so we need to have a new
        # array without the masked (clipped) elements.
        new_profile = clipped_profile[~clipped_profile.mask]

        # also the indexes are different
        new_x_axis = [i for i in range(len(clipped_profile)) if not clipped_profile.mask[i]]

        fitted_background = linear_fitter(linear_model, new_x_axis, new_profile)

        if plots:
            plt.title('Background Fitting Model Defined')
            plt.plot(median_profile, color='k')
            plt.plot(linear_model(range(ccd.data.shape[0])), color='r')
            plt.show()

            plt.title('Background Fitted Model')
            plt.plot(median_profile, color='k')
            plt.plot(fitted_background(range(ccd.data.shape[0])), color='r')
            plt.show()

        # Removing Background
        # Remove the background and set negative values to zero

        # build an array of the same dimensions of the profile
        background_array = fitted_background(range(len(median_profile)))

        background_subtracted = median_profile - background_array

        # set negative values to zero
        background_subtracted[background_subtracted < 0] = 0

        final_profile = background_subtracted.copy()

        # sigma clip and then get some features of the noise.
        clipped_final_profile = sigma_clip(final_profile, sigma=3, iters=3)

        # define the propper x-axis
        new_x_axis = [i for i in range(len(clipped_final_profile)) if
                      not clipped_final_profile.mask[i]]

        clipped_final_profile = clipped_final_profile[~clipped_final_profile.mask]



        background_level = np.abs(np.max(clipped_final_profile) - np.min(clipped_final_profile))
        # print('MEAN: ', np.mean(clipped_final_profile))
        # print('MEDIAN: ', np.median(clipped_final_profile))
        # print('STDDEV: ', np.std(clipped_final_profile))
        # print('RANGE: ', background_level)

        # TODO (simon): Add information to plots
        if plots:
            plt.ioff()
            plt.close()
            # if plt.isinteractive():
            #     plt.ioff()
            plt.title('Median Along Dispersion Axis (spatial)')
            plt.plot(background_subtracted, label='Background Subtracted Data')
            plt.plot(new_x_axis, clipped_final_profile, color='r', label='Sigma Clip Data')

            plt.axhline(background_level, color='m', label='Min-Max Difference')
            # plt.plot(final_profile, color='r')
            # plt.plot(median_profile)
            # plt.plot(background_array)
            plt.legend(loc='best')
            if plt.isinteractive():
                plt.draw()
                plt.pause(5)
            else:
                plt.show()

        # Identify targets
        # Now that the profile is flat it should be easier to identify targets.

        filtered_profile = np.where(np.abs(
            final_profile > final_profile.min() + 0.03 * final_profile.max()),
            final_profile,
            None)

        _upper_limit = final_profile.min() + 0.03 * final_profile.max()
        # print(_upper_limit)
        # print(np.median(final_profile))

        # replace the None elements to zero.
        none_to_zero_prof = [0 if it is None else it for it in filtered_profile]

        # convert the list to array
        filtered_profile = np.array(none_to_zero_prof)

        # find the peaks
        peaks = signal.argrelmax(filtered_profile, axis=0, order=order)[0]

        # find profile values for peaks found
        values = [final_profile[i] for i in peaks]

        # sort values and reverse the order so that larger values are first
        sorted_values = np.sort(values)[::-1]

        # pick nfind top values
        n_top_values = sorted_values[:nfind]
        # print(n_top_values)

        # retrieve the original index (real location) of the peaks
        selected_peaks = []
        for val in n_top_values:
            # TODO (simon): replace the 3 below by a parameter in a conf file.
            # discard peaks smaller than twice the level of background
            if val > 3 * background_level:
                index = np.where(values == val)[0]
                #     print(index[0])
                selected_peaks.append(peaks[index[0]])
            else:
                log_spec.debug('Discarding peak: {:.3f}'.format(val))

        if plots:
            plt.ioff()
            plt.plot(final_profile)
            plt.axhline(_upper_limit, color='g')
            for peak in selected_peaks:
                plt.axvline(peak, color='r')
            plt.show()

        # build the model to return
        fitter = fitting.LevMarLSQFitter()
        best_stddev = None

        profile_model = None
        for peak in selected_peaks:
            peak_value = median_profile[peak]
            gaussian = models.Gaussian1D(amplitude=peak_value,
                                         mean=peak,
                                         stddev=order).rename(
                'Gaussian_{:d}'.format(peak))

            # fixes mean and amplitude already found, just finding stddev
            gaussian.mean.fixed = True
            gaussian.amplitude.fixed = True
            fitted_gaussian = fitter(gaussian,
                                     range(len(median_profile)),
                                     median_profile)

            # after being fitted, unfix the parameters and now fix stddev
            fitted_gaussian.mean.fixed = False
            fitted_gaussian.amplitude.fixed = False
            fitted_gaussian.stddev.fixed = True

            # manually forcing the use of the best stddev if possitive
            # this disables the pipeline for extended sources
            if best_stddev is None:
                best_stddev = fitted_gaussian.stddev.value
            elif best_stddev < fitted_gaussian.stddev.value:
                fitted_gaussian.stddev.value = best_stddev
            else:
                best_stddev = fitted_gaussian.stddev.value
            if best_stddev < 0:
                best_stddev = None

            # print(fitted_gaussian.stddev.value)
            # plt.plot(median_profile, color='b')
            # plt.plot(fitted_gaussian(range(len(median_profile))), color='r')
            # plt.show()

            # this ensures the profile returned are valid
            if fitted_gaussian.stddev.value > 0:
                if profile_model is None:
                    profile_model = fitted_gaussian
                else:
                    profile_model += fitted_gaussian
        if plots:
            plt.plot(median_profile, color='b')
            plt.plot(profile_model(range(len(median_profile))), color='r')
            plt.show()

        # plt.imshow(ccd.data, clim=(50, 200), cmap='gray')
        # for peak in selected_peaks:
        #     plt.axhline(peak, color='r')
        # plt.show()

        if profile_model is None:
            return None
        else:
            return profile_model


def trace(ccd, model, trace_model, fitter, sampling_step, nsigmas=2):
    """Find the trace of a spectrum

    This function is called by the `trace_targets` targets, the difference is
    that it only takes single models only not CompoundModels so this function
    is called for every single target.

    Notes:
        This method forces the trace to go withing a rectangular region of
        center `model.mean.value` and width `2 * nsigmas`, this is for allowing
        the trace of low SNR targets. The assumption is valid since the spectra
        are always well aligned to the detectors's pixel columns. (dispersion
        axis)

    Args:
        ccd (object): A ccdproc.CCDData instance, 2D image.
        model (object): An astropy.modeling.Model instance that contains
            information regarding the target to be traced.
        trace_model (object): An astropy.modeling.Model instance, usually a low
            order polynomial.
        fitter (object): An astropy.modeling.fitting.Fitter instance. Will fit
            the sampled points to construct the trace model
        sampling_step (int): Step for sampling the spectrum.
        nsigmas (int): Number of stddev to each side of the mean to be used for
            searching the trace.

    Returns:
        An astropy.modeling.Model instance, that defines the trace of the
            spectrum.

    """
    spatial_length, dispersion_length = ccd.data.shape

    sampling_axis = range(0, dispersion_length, sampling_step)
    sample_values = []

    model_stddev = model.stddev.value
    model_mean = model.mean.value

    sample_center = float(model_mean)

    for point in sampling_axis:

        lower_limit = int(sample_center - nsigmas * model_stddev)
        upper_limit = int(sample_center + nsigmas * model_stddev)

        # print(sample_center, nsigmas, model_stddev, lower_limit, upper_limit)

        sample = ccd.data[lower_limit:upper_limit, point:point + sampling_step]
        sample_median = np.median(sample, axis=1)

        try:
            sample_peak = np.argmax(sample_median)
            # print(sample_peak + lower_limit)
        except ValueError:
            # plt.plot(model(range(spatial_length)))
            # plt.plot(ccd.data[:,point])
            # plt.show()
            print('Nsigmas ', nsigmas)
            print('Model Stddev ', model_stddev)
            print('sample_center ', sample_center)
            print('sample ', sample)
            print('sample_median ', sample_median)
            print('lower_limit ', lower_limit)
            print('upper_limit ', upper_limit)
            print('point ', point)
            print('point + sampling_step ', point + sampling_step)
            print(spatial_length, dispersion_length)
            sys.exit()

        sample_values.append(sample_peak + lower_limit)

        if np.abs(sample_peak + lower_limit - model_mean) < nsigmas * model_stddev:
            sample_center = int(sample_peak + lower_limit)
        else:
            # print(np.abs(sample_peak + lower_limit - model_mean), nsigmas * model_stddev)
            sample_center = float(model_mean)

    fitted_trace = fitter(trace_model, sampling_axis, sample_values)

    if False:
        plt.title(ccd.header['GSP_FNAM'])
        plt.imshow(ccd.data, clim=(30, 200))
        plt.plot(sampling_axis, sample_values, color='y', marker='o')
        plt.axhspan(lower_limit,
                    upper_limit,
                    alpha=0.4,
                    color='g')
        plt.plot(fitted_trace(range(dispersion_length)), color='c')
        # plt.plot(model(range(spatial_length)))
        if plt.isinteractive():
            plt.draw()
            plt.pause(2)
        else:
            plt.show()
        # print(dispersion_length)
        # print(sampling_axis)

    # fitted_trace = None

    return fitted_trace


def trace_targets(ccd, profile, sampling_step=5, pol_deg=2, plots=True):
    """Find the trace of the target's spectrum on the image

    This function defines a low order polynomial that trace the location of the
    spectrum. The attributes pol_deg and sampling_step define the polynomial
    degree and the spacing in pixels for the samples. For every sample a
    gaussian model is fitted and the center (mean) is recorded and since
    spectrum traces vary smoothly this value is used as a new center for the
    base model used to fit the spectrum profile.

    Notes:
        This doesn't work for extended sources. Also this calls for the function
        trace for doing the actual trace, the difference is that this method is
        at a higher level.

    Args:
        ccd (object): Instance of ccdproc.CCDData
        profile (object): Instance of astropy.modeling.Model, contains the
            spatial profile of the 2D spectrum.
        sampling_step (int): Frequency of sampling in pixels
        pol_deg (int): Polynomial degree for fitting the trace
        plots (bool): If True will show plots (debugging)

    Returns:
        all_traces (list): List that contains traces that are
            astropy.modeling.Model instance

    """

    # added two assert for debugging purposes
    assert isinstance(ccd, CCDData)
    assert isinstance(profile, Model)

    # Initialize model fitter
    model_fitter = fitting.LevMarLSQFitter()

    # Initialize the model to fit the traces
    trace_model = models.Polynomial1D(degree=pol_deg)

    # List that will contain all the Model instances corresponding to traced
    # targets
    all_traces = None

    if 'CompoundModel' in profile.__class__.name:
        log_spec.debug(profile.__class__.name)
        # TODO (simon): evaluate if targets are too close together.

        stddev_keys = [key for key in profile._param_names if 'stddev' in key]

        mean_keys = [key for key in profile._param_names if 'mean' in key]

        stddev_values = [
            profile.__getattribute__(key).value for key in stddev_keys]

        mean_values = [
            profile.__getattribute__(key).value for key in mean_keys]

        # if len(mean_values) == 1:
        #     nsigmas = 20
        # else:
        #     # get maximum width
        #     for i in range(len(mean_values)-1):


        for m in range(len(profile.submodel_names)):
            submodel_name = profile.submodel_names[m]

            model = profile[submodel_name]

            single_trace = trace(ccd=ccd,
                                 model=model,
                                 trace_model=trace_model,
                                 fitter=model_fitter,
                                 sampling_step=sampling_step)

            if all_traces is None:
                all_traces = [single_trace]
            else:
                all_traces.append(single_trace)

        return all_traces
    else:
        single_trace = trace(ccd=ccd,
                             model=profile,
                             trace_model=trace_model,
                             fitter=model_fitter,
                             sampling_step=sampling_step,
                             nsigmas=10)
        return [single_trace]


def get_extraction_zone(ccd,
                        extraction=None,
                        trace=None,
                        trace_index=None,
                        model=None,
                        n_sigma_extract=None,
                        zone=None,
                        plots=False):
    """Get the minimum rectangular CCD zone that fully contains the spectrum

    In this context, *fully contained* means the spectrum plus some region for
    background subtraction.

    Notes:
        For Goodman HTS the alignment of the spectrum with the detector lines
        is quite good, that's why this function does not consider the trace.
        Also because the `model` argument is based on the median throughout all
        the detector along the dispersion axis, so if there is a strong
        misalignment it will result in a wider Gaussian Profile

    Args:
        ccd (object): A ccdproc.CCDData instance, the image from which the zone
            will be extracted
        extraction (str): Extraction type, `simple` or `optimal`
        trace (object): An astropy.modeling.Model instance that correspond to
            the trace of the spectrum
        trace_index (int): The index number of the spectrum. 0 based.
        model (object): An astropy.modeling.Model instance that was previously
            fitted to the spatial profile.
        n_sigma_extract (int): Total number of sigmas to be extracted.
        plots (bool): If True will show plots, similar to a debugging mode.
        zone (list): Low and high limits to extract

    Returns:
        nccd (object): Instance of ccdproc.CCDData that contains only the region
            extracted from the full image. The header is updated with a new HISTORY
            keyword that contain the region of the original image extracted.
        model (object): Instance of astropy.modeling.Model with an updated mean
            to match the new center in pixel units.
        zone (list): Low and high limits of extraction zone

    """
    # make a copy of the ccd image
    nccd = ccd.copy()

    if zone is None and extraction is not None:
        assert (model is not None) and (n_sigma_extract is not None)
        assert isinstance(trace, Model)
        log_spec.debug('Extracting zone centered at: {:.3f}'.format(model.mean.value))

        spatial_length, dispersion_length = nccd.data.shape

        # get maximum variation in spatial direction
        trace_array = trace(range(dispersion_length))
        trace_inclination = trace_array.max() - trace_array.min()
        log_spec.debug('Trace Min-Max difference: {:.3f}'.format(trace_inclination))

        # m_mean = model.mean.value
        m_stddev = model.stddev.value
        extract_width = n_sigma_extract // 2 * m_stddev

        low_lim = np.max([0, int(trace_array.min() - extract_width)])

        hig_lim = np.min([int(trace_array.max() + extract_width),
                          spatial_length])

        # in some rare cases the low_lim turns out larger than hig_lim which
        # creates a series of problems regarding extraction, here I just reverse
        # them
        if low_lim > hig_lim:
            low_lim, hig_lim = hig_lim, low_lim

        log_spec.debug('Zone: low {:d}, high {:d}'.format(low_lim, hig_lim))
        zone = [low_lim, hig_lim]

        # This is to define the APNUM1 Keyword for the header.

        apnum_1 = '{:d} {:d} {:d} {:d}'.format(trace_index + 1,
                                              1,
                                              low_lim,
                                              hig_lim)

        nccd.header['APNUM1'] = apnum_1

        # this is necessary since we are cutting a piece of the full ccd.
        trace.c0.value -= low_lim
        log_spec.debug('Changing attribute c0 from trace, this is to adjust it to '
                 'the new extraction zone which is smaller that the full CCD.')

        log_spec.debug('Changing attribute mean of profile model')
        model.mean.value = extract_width



        nccd.data = nccd.data[low_lim:hig_lim, :]
        if nccd.mask is not None:
            log_spec.debug('Trimming mask')
            nccd.mask = nccd.mask[low_lim:hig_lim, :]
        nccd.header['HISTORY'] = 'Subsection of CCD ' \
                                 '[{:d}:{:d}, :]'.format(low_lim, hig_lim)

        if plots:
            plt.imshow(ccd.data, clim=(0, 60))
            plt.axhspan(low_lim, hig_lim, color='r', alpha=0.2)
            plt.show()

        # return nccd, trace, model, zone
        return zone

    else:

        low_lim, hig_lim = zone

        apnum_1 = '{:d} {:d} {:d} {:d}'.format(trace_index + 1,
                                               1,
                                               low_lim,
                                               hig_lim)

        nccd = ccd.copy()
        nccd.header['APNUM1'] = apnum_1

        nccd.data = nccd.data[low_lim:hig_lim, :]
        if nccd.mask is not None:
            log_spec.debug('Trimming mask')
            nccd.mask = nccd.mask[low_lim:hig_lim, :]
        nccd.header['HISTORY'] = 'Subsection of CCD ' \
                                 '[{:d}:{:d}, :]'.format(low_lim, hig_lim)
        return nccd


def add_wcs_keys(header):
    """Adds generic keyword to the header

    Linear wavelength solutions require a set of standard fits keywords. Later
    on they will be updated accordingly
    The main goal of putting them here is to have consistent and nicely ordered
    headers

    Notes:
        This does NOT add a WCS solution, just the keywords

    Args:
        header (object): New header without WCS entries

    Returns:
        header (object): Modified header with added WCS keywords

    """
    try:
        header['BANDID1'] = 'spectrum - background none, weights none, clean no'
        header['APNUM1'] = '1 1 0 0'
        header['WCSDIM'] = 1
        header['CTYPE1'] = 'LINEAR'
        header['CRVAL1'] = 1
        header['CRPIX1'] = 1
        header['CDELT1'] = 1
        header['CD1_1'] = 1
        header['LTM1_1'] = 1
        header['WAT0_001'] = 'system=equispec'
        header['WAT1_001'] = 'wtype=linear label=Wavelength units=angstroms'
        header['DC-FLAG'] = 0
        header['DCLOG1'] = 'REFSPEC1 = non set'
        return header
    except TypeError as err:
        log_spec.error("Can't add wcs keywords to header")
        log_spec.debug(err)


def remove_background_by_median(ccd, plots=False):
    """Remove Background of a ccd spectrum image

    Notes:
        This function works well for images without strong sky lines. Or for
        targets embedded in extended sources.

    Args:
        ccd (object): A ccdproc.CCDData instance.

    Returns:
        ccd (object): The modified ccdproc.CCDData instance.

    """
    new_ccd = ccd.copy()
    data = ma.masked_invalid(new_ccd.data)
    # x, y = ccd.data.shape
    median = ma.median(data, axis=0)

    data -= median
    data.set_fill_value(-np.inf)
    new_ccd.data = data.filled()

    # ccd.write('/user/simon/dummy_{:d}.fits'.format(g), clobber=True)
    return new_ccd


def get_background_value(background_image, zone_sep=0.5, zone=None):
    """finds background value for each dispersion line

    A background image is an image whose spectrum has been masked with zeros,
    the spectrum zone is retrieved by analyzing the pixels values and then two
    background zones are defined at a distance from the edges of the target zone
    defined by zone_sep, the width of the background zone is the same as the
    target's. After validating the background zone they averaged by median along
    the spatial direction. If there are two zones they are averaged as well.

    Args:
        background_image (object): ccdproc.CCDData instance. Spectra are masked.
        zone_sep (float): How far the background zone should be from the edges
            of the target zone. Default 0.5.
        zone (list): Alternative you could parse the zone as a list with each
            element a limit. Therefore there should be an even number of
            elements.

    Returns:
        A 1D array the same length of the dispersion length.

    """
    spatial_length, dispersion_length = background_image.data.shape
    if zone is None:
        background_profile = np.median(background_image.data, axis=1)

        target_zones_points = np.where(background_profile == 0)[0]

        target_zones = [i + 1 for i in range(len(target_zones_points) - 1) \
                        if
                        target_zones_points[i + 1] - target_zones_points[i] > 1]

        target_zones.append(np.argmin(target_zones_points))
        target_zones.append(len(target_zones_points))
        target_zones.sort()

        for i in range(len(target_zones) - 1):

            if target_zones[i + 1] - target_zones[i] > 2:

                log_spec.debug('Found valid target zone:'
                               ' {:d}:{:d}'.format(target_zones[i],
                                                   target_zones[i + 1]))

                full_zone = target_zones_points[
                            target_zones[i]:target_zones[i + 1]]
                zone = [full_zone[0], full_zone[-1]]
                # zone_width = len(zone)


            else:
                log_spec.warning('Target Zone is too small')

    # here do the evaluation if background zones are valid.
    first_background = None
    second_background = None
    background_median = None

    zone_width = zone[1] - zone[0]

    # first background zone
    back_first_low = int(zone[0] - (1 + zone_sep) * zone_width)
    back_first_high = int(zone[0] - zone_sep * zone_width)
    if 0 < back_first_low < back_first_high:

        first_background = np.median(background_image.data[
                                     back_first_low:back_first_high,
                                     :],
                                     axis=0)
    else:
        log_spec.debug('Zone [{:d}:{:d}] is forbidden'.format(
            back_first_low,
            back_first_high))

    # second background zone
    back_second_low = int(zone[-1] + zone_sep * zone_width)
    back_second_high = int(zone[-1] + (1 + zone_sep) * zone_width)
    if back_second_low < back_second_high < spatial_length:

        second_background = np.mean(background_image.data[
                                    back_second_low:back_second_high,
                                    :],
                                    axis=0)
    else:
        log_spec.debug('Zone [{:d}:{:d}] is forbidden'.format(
            back_second_low,
            back_second_high))

    if first_background is not None and second_background is not None:
        background_mean = np.mean([first_background, second_background], axis=0)

        return background_mean
    elif first_background is not None and second_background is None:
        return first_background
    elif first_background is None and second_background is not None:
        return second_background
    else:
        log_spec.error(
            'Not possible to get a background extraction zone')
        return 0


def create_background_image(ccd, profile_model, nsigma, separation):
    """Creates a background-only image

    Using a profile model and assuming the spectrum is misaligned only a little
    bit (i.e. a couple of pixels from end to end) with respect to the lines of
    the detector. The number of sigmas determines the width of the zone to be
    masked and the separation is an offset that is added.

    Args:
        ccd (object): A ccdproc.CCDData instance.
        profile_model (object): An astropy.modeling.Model instance. Describes
            the spatial profile of the target.
        nsigma (float): Number of sigmas. Used to calculate the width of the
            target zone to be masked.
        separation (float): Additional offset that adds to the width of the
            target zone.

    Returns:
        A ccdproc.CCDData instance with the spectrum masked.

    """
    background_ccd = ccd.copy()
    spatial_length, dispersion_length = background_ccd.data.shape
    target_profiles = []
    if 'CompoundModel' in profile_model.__class__.name:
        log_spec.debug(profile_model.submodel_names)
        for m in range(len(profile_model.submodel_names)):
            submodel_name = profile_model.submodel_names[m]

            target_profiles.append(profile_model[submodel_name])
    else:
        target_profiles.append(profile_model)

    for target in target_profiles:
        target_mean = target.mean.value
        target_stddev = target.stddev.value

        data_low_lim = np.int(np.max(
            [0, target_mean - (nsigma / 2. + separation) * target_stddev]))

        data_high_lim = np.int(np.min([spatial_length, int(
            target_mean + (nsigma / 2. + separation) * target_stddev)]))

        background_ccd.data[data_low_lim:data_high_lim, :] = 0

    if False:
        plt.title('Background Image')
        plt.imshow(background_ccd.data, clim=(0, 50))
        plt.show()

    return background_ccd


def extract(ccd,
            trace,
            spatial_profile,
            extraction,
            zone,
            background_level=0,
            sampling_step=1,
            plots=False):
    """Performs spectrum extraction

    This function is designed to perform two types of spectrum extraction, a
    simple sum in the spatial direction and an optimal extraction.

    Notes:
        For the beta release the optimal extraction is not implemented.

    Args:
        ccd (object): Instance of ccdproc.CCDData containing a 2D spectrum
        trace (object): Instance of astropy.modeling.Model, a low order
            polynomial that defines the trace of the spectrum in the ccd object.
        spatial_profile (object): Instance of astropy.modeling.Model, a Gaussian
            model previously fitted to the spatial profile of the 2D spectrum
            contained in the ccd object.
        extraction (str): Extraction type, can be `simple` or `optimal` for the
            beta release the optimal extraction is not implemented yet.
        zone (list):
        background_level (array):
        sampling_step (int): The optimal extraction needs to sample the spatial
            profile, this value defines the intervals at which get such
            sampling.
        plots (bool): Determines whether display plots or not.

    Returns:
        ccd (object): Instance of ccdproc.CCDData containing a 1D spectrum. The
            attribute 'data' is replaced by the 1D array resulted from the
            extraction process.

    Raises:
        NotImplementedError: When `extraction` is optimal, this is valid for the
            beta release

    """
    assert isinstance(ccd, CCDData)
    assert isinstance(trace, Model)

    nccd = ccd.copy()

    spatial_length, dispersion_length = nccd.data.shape

    apnum1 = '{:d} {:d} {:d} {:d}'.format(1, 1, zone[0], zone[1])
    log_spec.debug('APNUM1 Keyword: {:s}'.format(apnum1))

    # create variance model
    rdnoise = float(nccd.header['RDNOISE'])
    gain = float(nccd.header['GAIN'])
    log_spec.debug('Original Name {:s}'.format(nccd.header['GSP_FNAM']))

    variance_2d = (rdnoise + np.absolute(nccd.data) * gain) / gain
    cr_mask = np.ones(nccd.data.shape, dtype=int)
    # if nccd.mask is None and nccd.header['OBSTYPE'] == 'OBJECT':
    #     log_spec.debug('Finding cosmic rays to create mask')
    #     cr_mask = cosmicray_rejection(ccd=ccd, mask_only=True)
    #     cr_mask = np.log_spec.cal_not(cr_mask).astype(int)
    # elif nccd.mask is None and nccd.header['OBSTYPE'] != 'OBJECT':
    #     log_spec.debug('Only OBSTYPE == OBJECT get cosmic ray rejection.')
    #     cr_mask = np.ones(nccd.data.shape, dtype=int)
    # else:
    #     log_spec.debug('Cosmic ray mask already exists.')
    #     cr_mask = np.logical_not(nccd.mask).astype(int)

    model_fitter = fitting.LevMarLSQFitter()

    # print(spatial_profile.mean.value)
    # print(trace.c0.value)

    if isinstance(spatial_profile, models.Gaussian1D):
        amplitude = spatial_profile.amplitude.value
        mean = spatial_profile.mean.value
        stddev = spatial_profile.stddev.value
        new_model = models.Gaussian1D(amplitude=amplitude,
                                      mean=mean,
                                      stddev=stddev)
        # print('Fixed ', new_model.mean.fixed)
    else:
        raise NotImplementedError
    log_spec.debug('nccd.data is a masked array: '
                   '{:s}'.format(str(np.ma.isMaskedArray(nccd.data))))

    nccd.data = np.ma.masked_invalid(nccd.data)
    # print(np.ma.isMaskedArray(nccd.data))
    np.ma.set_fill_value(nccd.data, 0)

    if extraction == 'simple':


        # print(indexes)
        if plots:
            indexes = np.argwhere(cr_mask == 0)
            fig = plt.figure(1)
            ax1 = fig.add_subplot(111)
            for index in indexes:
                x, y = index
                ax1.plot(y, x, marker='o', color='r')

            ax1.imshow(nccd.data, interpolation='none')
            if plt.isinteractive():
                plt.draw()
                plt.pause(1)
            else:
                plt.show()

        # print(np.ma.isMaskedArray(nccd.data))
        # spectrum zone limit
        low_lim, high_lim = zone
        spectrum_masked = nccd.data * cr_mask
        # plt.imshow(spectrum_masked, clim=(10, 70))
        # plt.show()
        # TODO (simon): Add fractional pixel
        spectrum_sum = np.ma.sum(spectrum_masked[low_lim:high_lim, :], axis=0)

        background_sum = np.abs(high_lim - low_lim) * background_level

        nccd.data = spectrum_sum - background_sum

        nccd.header['APNUM1'] = apnum1

        if plots:
            fig = plt.figure()
            fig.canvas.set_window_title('Simple Extraction')
            # ax = fig.add_subplot(111)
            manager = plt.get_current_fig_manager()
            if plt.get_backend() == u'GTK3Agg':
                manager.window.maximize()
            elif plt.get_backend() == u'Qt5Agg':
                manager.window.showMaximized()

            plt.title(nccd.header['OBJECT'])
            plt.xlabel('Dispersion Axis (Pixels)')
            plt.ylabel('Intensity (Counts)')
            # plt.plot(simple_sum, label='Simple Sum', color='k', alpha=0.5)
            plt.plot(nccd.data, color='k',
                     label='Simple Extracted')
            plt.plot(background_sum, color='r', label='Background')
            plt.plot(spectrum_sum, color='b', label='Spectrum Raw Sum')
            plt.xlim((0, len(nccd.data)))
            plt.legend(loc='best')
            if plt.isinteractive():
                plt.draw()
                plt.pause(1)
            else:
                plt.show()

    elif extraction == 'optimal':
        raise NotImplementedError
        # out_spectrum = np.empty(dispersion_length)
        # for i in range(0, dispersion_length, sampling_step):
        #     # force the model to follow the trace
        #     new_model.mean.value = trace(i)
        #
        #     # warn if the difference of the spectrum position in the trace at the
        #     # extremes of the sampling range is larger than 1 pixel.
        #     if np.abs(trace(i) - trace(i + sampling_step)) > 1:
        #         log_spec.warning('Sampling step might be too large')
        #
        #     sample = np.median(nccd.data[:, i:i + sampling_step], axis=1)
        #     fitted_profile = model_fitter(model=new_model,
        #                                   x=range(len(sample)),
        #                                   y=sample)
        #
        #     profile = fitted_profile(range(sample.size))
        #
        #     # enforce positivity
        #     pos_profile = np.array([np.max([0, x]) for x in profile])
        #
        #     # enforce normalization
        #     nor_profile = np.array([x / pos_profile.sum() for x in pos_profile])
        #
        #     if sampling_step > 1:
        #         # TODO (simon): Simplify to Pythonic way
        #         right = min((i + sampling_step), dispersion_length)
        #
        #         for e in range(i, right, 1):
        #             mask = cr_mask[:, e]
        #             data = ma.masked_invalid(nccd.data[:, e])
        #             # print(ma.isMaskedArray(data))
        #             V = variance_2d[:, e]
        #             P = nor_profile
        #             a = [(P[z] / V[z]) / np.sum(P ** 2 / V) for z in
        #                  range(P.size)]
        #             weights = (nor_profile / variance_2d[:, e]) / np.sum(
        #                 nor_profile ** 2 / variance_2d[:, e])
        #             # print('SUMN ', np.sum(a), np.sum(weights), np.sum(nor_profile), np.sum(P * weights))
        #
        #
        #
        #             # if e in range(5, 4001, 500):
        #             #     plt.plot(nor_profile * data.max()/ nor_profile.max(), label=str(e))
        #             #     plt.plot(data, label='Data')
        #             #     plt.legend(loc='best')
        #             #     plt.show()
        #
        #             out_spectrum[e] = np.sum(data * mask * nor_profile)
        # nccd.data = out_spectrum
        # if plots:
        #     fig = plt.figure()
        #     fig.canvas.set_window_title('Optimal Extraction')
        #     # ax = fig.add_subplot(111)
        #     manager = plt.get_current_fig_manager()
        #
        #     if plt.get_backend() == u'GTK3Agg':
        #         manager.window.maximize()
        #     elif plt.get_backend() == u'Qt5Agg':
        #         manager.window.showMaximized()
        #
        #     plt.title(nccd.header['OBJECT'])
        #     plt.xlabel('Dispersion Axis (Pixels)')
        #     plt.ylabel('Intensity (Counts)')
        #     # plt.plot(simple_sum, label='Simple Sum', color='k', alpha=0.5)
        #     plt.plot(nccd.data, color='k',
        #              label='Optimal Extracted')
        #     plt.xlim((0, len(nccd.data)))
        #     plt.legend(loc='best')
        #     if plt.isinteractive():
        #         plt.draw()
        #         plt.pause(1)
        #     else:
        #         plt.show()

    # nccd.data = out_spectrum
    return nccd


# classes definition

class NightDataContainer(object):
    """This class is designed to be the organized data container. It doesn't
    store image data but list of pandas.DataFrame objects. Also it stores
    critical variables such as sunrise and sunset times.

    """

    def __init__(self, path, instrument, technique):
        """Initializes all the variables for the class

        Args:
            path (str): Full path to the directory where raw data is located
            instrument (str): `Red` or `Blue` stating whether the data was taken
                using the Red or Blue Goodman Camera.
            technique (str): `Spectroscopy` or `Imaging` stating what kind of
                data was taken.
        """

        self.full_path = path
        self.instrument = instrument
        self.technique = technique
        self.is_empty = True

        """For imaging use"""
        self.bias = None
        self.day_flats = None
        self.dome_flats = None
        self.sky_flats = None
        self.data_groups = None

        """For spectroscopy use"""

        # comp_groups will store pandas.DataFrame (groups) that contain only
        # OBSTYPE == COMP, they should be requested only when needed, for the
        # science case when for every science target is observed with comparison
        # lamps and quartz (if)
        self.comp_groups = None

        # object_groups will store pandas.DataFrame (groups) with only
        # OBSTYPE == OBJECT this is the case when the observer takes comparison
        # lamps only at the beginning or end of the night.
        self.object_groups = None

        # spec_groups will store pandas.DataFrame (groups) with a set of OBJECT
        # and COMP, this is usually the case for radial velocity studies.
        self.spec_groups = None

        """Time reference points"""
        self.sun_set_time = None
        self.sun_rise_time = None
        self.evening_twilight = None
        self.morning_twilight = None

    def add_bias(self, bias_group):
        """Adds a bias group

        Args:
            bias_group (pandas.DataFrame): Contains a set of keyword values of
                grouped image metadata

        """

        if len(bias_group) < 2:
            if self.technique == 'Imaging':

                log_ccd.error('Imaging mode needs BIAS to work properly. '
                          'Go find some.')

            else:
                log_ccd.warning('BIAS are needed for optimal results.')
        else:
            if self.bias is None:
                self.bias = [bias_group]
            else:
                self.bias.append(bias_group)
        if self.bias is not None:
            self.is_empty = False

    def add_day_flats(self, day_flats):
        """"Adds a daytime flat group

        Args:
            day_flats (pandas.DataFrame): Contains a set of keyword values of
                grouped image metadata

        """

        if self.day_flats is None:
            self.day_flats = [day_flats]
        else:
            self.day_flats.append(day_flats)
        if self.day_flats is not None:
            self.is_empty = False

    def add_data_group(self, data_group):
        """Adds a data group

        Args:
            data_group (pandas.DataFrame): Contains a set of keyword values of
                grouped image metadata

        """

        if self.data_groups is None:
            self.data_groups = [data_group]
        else:
            self.data_groups.append(data_group)
        if self.data_groups is not None:
            self.is_empty = False

    def add_comp_group(self, comp_group):
        """Adds a comp-only group

        Args:
            comp_group (pandas.DataFrame): Contains a set of keyword values of
                grouped image metadata

        """

        if self.comp_groups is None:
            self.comp_groups = [comp_group]
        else:
            self.comp_groups.append(comp_group)
        if self.comp_groups is not None:
            self.is_empty = False

    def add_object_group(self, object_group):
        """Adds a object-only group

        Args:
            object_group (pandas.DataFrame): Contains a set of keyword values of
                grouped image metadata

        """

        if self.object_groups is None:
            self.object_groups = [object_group]
        else:
            self.object_groups.append(object_group)
        if self.object_groups is not None:
            self.is_empty = False

    def add_spec_group(self, spec_group):
        """Adds a data group containing object and comp

        Args:
            spec_group (pandas.DataFrame): Contains a set of keyword values of
                grouped image metadata

        """

        if self.spec_groups is None:
            self.spec_groups = [spec_group]
        else:
            self.spec_groups.append(spec_group)
        if self.spec_groups is not None:
            self.is_empty = False

    def set_sun_times(self, sun_set, sun_rise):
        """Sets values for sunset and sunrise

        Args:
            sun_set (str): Sun set time in the format 'YYYY-MM-DDTHH:MM:SS.SS'
            sun_rise (str):Sun rise time in the format 'YYYY-MM-DDTHH:MM:SS.SS'

        """

        self.sun_set_time = sun_set
        self.sun_rise_time = sun_rise

    def set_twilight_times(self, evening, morning):
        """Sets values for evening and morning twilight

        Args:
            evening (str): Evening twilight time in the format
                'YYYY-MM-DDTHH:MM:SS.SS'
            morning (str): Morning twilight time in the format
                'YYYY-MM-DDTHH:MM:SS.SS'

        """

        self.evening_twilight = evening
        self.morning_twilight = morning


class SpectroscopicMode(object):

    def __init__(self):
        """Init method for the Spectroscopic Mode

        This method defines a pandas.DataFrame instance that contains all the
        current standard wavelength modes for Goodman HTS.

        """
        columns = ['grating_freq', 'wavmode', 'camtarg', 'grttarg', 'ob_filter']
        spec_mode = [['400', 'm1', '11.6', '5.8', 'None'],
                     ['400', 'm2', '16.1', '7.5', 'GG455'],
                     ['600', 'UV', '15.25', '7.0', 'None'],
                     ['600', 'Blue', '17.0', '7.0', 'None'],
                     ['600', 'Mid', '20.0', '10.0', 'GG385'],
                     ['600', 'Red', '27.0', '12.0', 'GG495'],
                     ['930', 'm1', '20.6', '10.3', 'None'],
                     ['930', 'm2', '25.2', '12.6', 'None'],
                     ['930', 'm3', '29.9', '15.0', 'GG385'],
                     ['930', 'm4', '34.6', '18.3', 'GG495'],
                     ['930', 'm5', '39.4', '19.7', 'GG495'],
                     ['930', 'm6', '44.2', '22.1', 'OG570'],
                     ['1200', 'm0', '26.0', '16.3', 'None'],
                     ['1200', 'm1', '29.5', '16.3', 'None'],
                     ['1200', 'm2', '34.4', '18.7', 'None'],
                     ['1200', 'm3', '39.4', '20.2', 'None'],
                     ['1200', 'm4', '44.4', '22.2', 'GG455'],
                     ['1200', 'm5', '49.6', '24.8', 'GG455'],
                     ['1200', 'm6', '54.8', '27.4', 'GG495'],
                     ['1200', 'm7', '60.2', '30.1', 'OG570'],
                     ['1800', 'Custom', 'None', 'None', 'None'],
                     ['2100', 'Custom', 'None', 'None', 'None'],
                     ['2400', 'Custom', 'None', 'None', 'None']
                     ]
        self.modes_data_frame = pandas.DataFrame(spec_mode, columns=columns)
        # print(self.modes_data_frame)

    def __call__(self,
                 header=None,
                 grating=None,
                 camera_targ=None,
                 grating_targ=None,
                 blocking_filter=None):
        """Get spectroscopic mode

        This method can be called either parsing a header alone or the rest of
        values separated.
        Args:
            header (object): FITS header.
            grating (str): Grating as in the FITS header.
            camera_targ (str): Camera target angle as in the FITS header.
            grating_targ (str): Grating target angle as in the FITS header.
            blocking_filter (str): Order blocking filter as in the FITS header.

        Returns:
            string that defines the instrument wavelength mode.

        """

        if all(x is None for x in (
                grating, camera_targ, grating_targ, blocking_filter)) and \
                        header is not None:

            grating = str(re.sub('[A-Za-z_-]', '', header['grating']))
            camera_targ = str(header['cam_targ'])
            grating_targ = str(header['grt_targ'])
            blocking_filter = str(header['filter2'])
            # print(grating, camera_targ, grating_targ, blocking_filter)

            return self.get_mode(grating=grating,
                                 camera_targ=camera_targ,
                                 grating_targ=grating_targ,
                                 blocking_filter=blocking_filter)
        else:
            grating = re.sub('[A-Za-z_-]', '', grating)

            return self.get_mode(grating=grating,
                                 camera_targ=camera_targ,
                                 grating_targ=grating_targ,
                                 blocking_filter=blocking_filter)

    def get_mode(self, grating, camera_targ, grating_targ, blocking_filter):
        """Get the camera's optical configuration mode.

        This method is useful for data that does not have the WAVMODE keyword

        Args:
            grating (str): Grating frequency as string
            camera_targ (str): Camera target angle as in the header.
            grating_targ (str): Grating target angle as in the header.
            blocking_filter (str): Order blocking filter listed on the header.

        Returns:
            string that defines the wavelength mode used

        """

        # print(grating, camera_targ, grating_targ)
        if any(grat == grating for grat in ('1800', '2100', '2400')):
            central_wavelength = get_central_wavelength(grating=grating,
                                                        grt_ang=grating_targ,
                                                        cam_ang=camera_targ)
            return 'Custom_{:d}nm'.format(int(round(central_wavelength)))

        else:
            _mode = self.modes_data_frame[
                ((self.modes_data_frame['grating_freq'] == grating) &
                 (self.modes_data_frame['camtarg'] == camera_targ) &
                 (self.modes_data_frame['grttarg'] == grating_targ))]
            if _mode.empty:
                central_wavelength = get_central_wavelength(grating=grating,
                                                            grt_ang=grating_targ,
                                                            cam_ang=camera_targ)
                return 'Custom_{:d}nm'.format(int(round(central_wavelength)))
            else:
                # print('%s %s' %(grating, _mode.wavmode))
                # print(_mode['wavmode'].to_string)
                return _mode['wavmode'].to_string(index=False)

class NoTargetException(Exception):
    """Exception to be raised when no target is identified"""
    def __init__(self):
        Exception.__init__(self, 'No targets identified.')


class NoMatchFound(Exception):
    def __init__(self):
        Exception.__init__(self, 'Did not find a match')


class NotEnoughLinesDetected(Exception):
    def __init__(self):
        Exception.__init__(self, 'Not enough lines detected.')


class CriticalError(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)
