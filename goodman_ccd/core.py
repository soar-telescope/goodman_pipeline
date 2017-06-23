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
import shutil
import subprocess
from threading import Timer

matplotlib.use('GTK3Agg')
from matplotlib import pyplot as plt
from ccdproc import CCDData, ImageFileCollection
from astropy.coordinates import EarthLocation
from astropy.time import Time, TimeDelta
from astropy.stats import sigma_clip
from astroplan import Observer
from astropy import units as u
from astropy.modeling import (models, fitting, Model)

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
    log.info('Finding duplicated keywords')
    log.warning('Files will be overwritten')
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
        log.info('Found {:d} duplicated keyword '
                 '{:s}'.format(len(multiple_keys),
                               's' if len(multiple_keys) > 1 else ''))

        for image_file in files:
            log.debug('Processing Image File: {:s}'.format(image_file))
            try:
                ccd = CCDData.read(image_file)
                for keyword in multiple_keys:
                    while ccd.header.count(keyword) > 1:
                        ccd.header.remove(keyword,
                                          ccd.header.count(keyword) - 1)
                log.warning('Overwriting file with duplicated keywords removed')
                log.debug('File %s overwritten', image_file)
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
        [abs(ccd_section_median[i + 1] - ccd_section_median[i])
         for i in range(0, len(ccd_section_median) - 1)])

    filtered_data = np.where(
        np.abs(pseudo_derivative > 0.5 * pseudo_derivative.max()),
        pseudo_derivative,
        None)

    peaks = signal.argrelmax(filtered_data, axis=0, order=3)[0]

    slit_trim_section = None
    if len(peaks) > 2 or peaks == []:
        log.debug('No trim section')
    else:
        if len(peaks) == 2:
            # This is the ideal case, when the two edges of the slit are found.
            low, high = peaks
            slit_trim_section = '[:,{:d}:{:d}]'.format(low, high)
        elif len(peaks) == 1:
            # when only one peak is found it will choose the largest region from
            # the spatial axis center to one edge.
            if peaks[0] <= spatial_half:
                slit_trim_section = '[:,{:d}:{:d}]' \
                                    ''.format(peaks[0], len(ccd_section_median))
            else:
                slit_trim_section = '[:,{:d}:{:d}]'.format(0, peaks[0])
    return slit_trim_section


def dcr_cosmicray_rejection(data_path, in_file, prefix, dcr_par_dir,
                            delete=False):
    """Runs an external code for cosmic ray rejection

    DCR was created by Wojtek Pych and the code can be obtained from
    http://users.camk.edu.pl/pych/DCR/ and is written in C

    Args:
        data_path (str): Data location
        in_file (str): Name of the file to have its cosmic rays removed
        prefix (str): Prefix to add to the file with the cosmic rays removed
        dcr_par_dir (str): Directory of default dcr.par file
        delete (bool): True for deleting the input and cosmic ray file.

    Returns:

    """

    log.info('Removing cosmic rays using DCR by Wojtek Pych')
    log.info('See http://users.camk.edu.pl/pych/DCR/')

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

    log.debug('DCR command:')
    log.debug(command)
    # print(command.split(' '))

    # get the current working directory to go back to it later in case the
    # the pipeline has not been called from the same data directory.
    cwd = os.getcwd()

    # move to the directory were the data is, dcr is expecting a file dcr.par
    os.chdir(data_path)

    # check if file dcr.par exists
    while not os.path.isfile('dcr.par'):

        log.error('File dcr.par does not exist. Copying default one.')
        dcr_par_path = os.path.join(dcr_par_dir, 'dcr.par')
        log.debug('dcr.par full path: {:s}'.format(dcr_par_path))
        if os.path.isfile(dcr_par_path):
            shutil.copy2(dcr_par_path, data_path)
        else:
            log.error('Could not find dcr.par file')
    else:
        log.info('File dcr.par exists.')

    # call dcr
    try:

        dcr = subprocess.Popen(command.split(),
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)

    except OSError as error:
        log.error(error)
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

    # read stdout and stderr
    # stdout = dcr.stdout.read()
    # stderr = dcr.stderr.read()

    # If no error stderr is an empty string
    if stderr != '':
        log.error(stderr)
        if 'dcr: not found' in stderr:
            sys.exit('Your system can not locate the executable file dcr, try '
                 'moving it to /bin or create a symbolic link\n\n\tcd /bin\n\t'
                     'sudo ln -s /full/path/to/dcr')
    else:
        for output_line in stdout.split('\n'):
            log.info(output_line)

    # delete extra files only if the execution ended without error
    if delete and stderr == '' and 'USAGE:' not in stdout:
        try:
            log.warning('Removing input file: {:s}'.format(full_path_in))
            os.unlink(full_path_in)
        except OSError as error:
            log.error(error)

        try:
            log.warning(
                'Removing cosmic rays file: {:s}'.format(full_path_cosmic))
            os.unlink(full_path_cosmic)
        except OSError as error:
            log.error(error)



def cosmicray_rejection(ccd, mask_only=False):
    """Do cosmic ray rejection

    This function in fact does not apply any correction, it detects the cosmic
    rays and updates the attribute mask of the ccd object (CCDData instance).
    The attribute mask is used later as a mask for the pixels hit by cosmic rays

      Notes:
          OBS: cosmic ray rejection is working pretty well by defining gain = 1.
          It's not working when we use the real gain of the image. In this case
          the sky level changes by a factor equal the gain.
          Function to determine the sigfrac and objlim: y = 0.16 * exptime + 1.2

      Args:
          ccd (object): CCDData Object
          mask_only (bool): In some cases you may want to obtain the cosmic
          rays mask only

      Returns:
          ccd (object): A CCDData instance with the mask attribute updated.

      """
    # TODO (simon): Validate this method
    if ccd.header['OBSTYPE'] == 'OBJECT':
        value = 0.16 * float(ccd.header['EXPTIME']) + 1.2
        log.info('Cleaning cosmic rays... ')

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

        ccd.header.add_history("Cosmic rays rejected with LACosmic")
        log.info("Cosmic rays rejected with LACosmic")
        if mask_only:
            return ccd.mask
        else:
            return ccd
    else:
        log.info('Skipping cosmic ray rejection for image of datatype: '
                 '{:s}'.format(ccd.header['OBSTYPE']))
        return ccd


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
        flat_name (str): Full path of masterflat basename. Ends in '*.fits'.

    Returns:
        master_flat (object): A ccdproc.CCDData instance
        master_flat_name (str): Full path to the chosen masterflat.

    """
    flat_list = glob.glob(flat_name)
    log.debug('Flat base name {:s}'.format(flat_name))
    log.debug('Matching master flats found: {:d}'.format(len(flat_list)))
    if len(flat_list) > 0:
        if len(flat_list) == 1:
            master_flat_name = flat_list[0]
        else:
            master_flat_name = flat_list[0]
        # elif any('dome' in flat for flat in flat_list):
        #     master_flat_name =

        master_flat = CCDData.read(master_flat_name, unit=u.adu)
        log.info('Found suitable master flat: {:s}'.format(master_flat_name))
        return master_flat, master_flat_name
    else:
        log.error('There is no flat available')
        return None, None


def print_default_args(args):
    """Print default values of arguments.

    This is mostly helpful for debug but people not familiar with the software
    might find it useful as well

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
        log.info('Default value for {:s} is {:s}'.format(arg_name[key],
                                                         str(
                                                             args.__getattribute__(
                                                                 key))))


def normalize_master_flat(master, name, method='simple', order=15):
    """ Master flat normalization method

    This function normalize a master flat in three possible ways:
     _mean_: simply divide the data by its mean

     _simple_: Calculates the median along the spatial axis in order to obtain
     the dispersion profile. Then fits a Chebyshev1D model and apply this to all
     the data.

     _full_: This is for experimental purposes only because it takes a lot of
     time to process. It will fit a model to each line along the dispersion axis
     and then divide it by the fitted model. I do not recommend this method
     unless you have a good reason as well as a powerful computer.

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
        log.info('Normalizing by mean')
        master.data /= master.data.mean()

        master.header.add_history('Flat Normalized by Mean')

    elif method == 'simple' or method == 'full':
        log.info('Normalizing flat by {:s} model'.format(method))

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

            master.header.add_history('Flat Normalized by simple model')

        elif method == 'full':
            log.warning('This part of the code was left here for experimental '
                        'purposes only')
            log.info('This procedure takes a lot to process, you might want to'
                     'see other method such as simple or mean.')
            for i in range(x_size):
                fit = model_fitter(model_init, x_axis, master.data[i])
                master.data[i] = master.data[i] / fit(x_axis)
            master.header.add_history('Flat Normalized by full model')

    # write normalized flat to a file
    master.write(norm_name, clobber=True)

    return master


def get_central_wavelength(grating, grt_ang, cam_ang):
    """Calculates the central wavelength for a given spectroscopic mode

    The equation used to calculate the central wavelength is the following


    $$\lambda_{central} = \frac{1e6}{GRAT}
    \sin\left(\frac{\alpha \pi}{180}\right) +
    \sin\left(\frac{\beta \pi}{180}\right)$$

    Args:
        grating:
        grt_ang:
        cam_ang:

    Returns:

    """

    grating_frequency = float(grating)
    alpha = float(grt_ang)
    beta = float(cam_ang) - float(grt_ang)

    central_wavelength = (1e6 / grating_frequency) * \
                         (np.sin(alpha * np.pi / 180.) +
                          np.sin(beta * np.pi / 180.))

    log.debug('Found {:.3f} as central wavelength'.format(central_wavelength))
    return central_wavelength


# spectroscopy specific functions

def classify_spectroscopic_data(path, search_pattern):
    """Classify data by grouping them as pandas.DataFrame instances

    This functions uses ImageFileCollection from ccdproc. First it creates a
    collection of information regarding the images located in _path_ that match
    the pattern _search_pattern_
    The information obtained are all keywords listed in the list _keywords_
    The ImageFileCollection is translated into pandas.DataFrame and then is used
    much like SQL database to select and filter values and in that way put them
    in groups that are pandas.DataFrame instances.


    Args:
        path (str): Path to data location
        search_pattern (str): Prefix to match files.

    Returns:
        data_container (object): Instance of NightDataContainer

    """

    search_path = os.path.join(path, search_pattern + '*.fits')

    data_container = NightDataContainer(path=path,
                                        instrument=str('Red'),
                                        technique=str('Spectroscopy'))

    file_list = glob.glob(search_path)

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

        data_container.add_spec_group(spec_group=spec_group)

    return data_container


def spectroscopic_extraction(ccd, extraction,
                             comp_list=None,
                             n_sigma_extract=10,
                             plots=False):
    """This function does not do the actual extraction but prepares the data



    Args:
        ccd (object): A ccdproc.CCDData Instance
        extraction (str): Extraction type name. _simple_ or _optimal_
        comp_list (list): List of ccdproc.CCDData instances of COMP lamps data
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

    iccd = remove_background(ccd=ccd)

    profile_model = identify_targets(ccd=iccd, plots=plots)
    del (iccd)

    if isinstance(profile_model, Model):
        traces = trace_targets(ccd=ccd, profile=profile_model, plots=plots)
        # extract(ccd=ccd,
        #         spatial_profile=profile_model,
        #         n_sigma_extract=10,
        #         sampling_step=5)
        if 'CompoundModel' in profile_model.__class__.name:
            log.debug(profile_model.submodel_names)
            for m in range(len(profile_model.submodel_names)):

                submodel_name = profile_model.submodel_names[m]

                nccd, trace, model, zone = get_extraction_zone(
                    ccd=ccd,
                    extraction=extraction,
                    trace=traces[m],
                    model=profile_model[submodel_name],
                    n_sigma_extract=n_sigma_extract,
                    plots=plots)

                for comp in comp_list:
                    comp_zone = get_extraction_zone(ccd=comp,
                                                    extraction=extraction,
                                                    trace=trace,
                                                    zone=zone,
                                                    plots=plots)
                    # since a comparison lamp only needs only the relative line
                    # center in the dispersion direction, therefore the flux is
                    # not important we are only calculating the median along the
                    # spatial direction
                    comp_zone.data = np.median(comp_zone.data, axis=0)
                    comp_zones.append(comp_zone)

                nccd = remove_background(ccd=nccd,
                                         plots=plots)

                extracted_ccd = extract(ccd=nccd,
                                        trace=trace,
                                        spatial_profile=model,
                                        extraction=extraction,
                                        sampling_step=10,
                                        plots=plots)
                extracted.append(extracted_ccd)

                if plots:
                    plt.imshow(nccd)
                    plt.show()

        else:
            nccd, trace, model, zone = get_extraction_zone(
                ccd=ccd,
                extraction=extraction,
                trace=traces[0],
                model=profile_model,
                n_sigma_extract=n_sigma_extract,
                plots=plots)

            for comp in comp_list:
                comp_zone = get_extraction_zone(ccd=comp,
                                                extraction=extraction,
                                                trace=trace,
                                                zone=zone,
                                                plots=plots)

                # since a comparison lamp only needs the relative line
                # center in the dispersion direction, therefore the flux is not
                # important we are only calculating the median along the spatial
                # direction
                comp_zone.data = np.median(comp_zone.data, axis=0)
                comp_zones.append(comp_zone)

            nccd = remove_background(ccd=nccd,
                                     plots=plots)

            extracted_ccd = extract(ccd=nccd,
                                    trace=trace,
                                    spatial_profile=model,
                                    extraction=extraction,
                                    sampling_step=10,
                                    plots=plots)

            extracted.append(extracted_ccd)

            if plots:
                plt.imshow(nccd)
                plt.show()

        # print(extracted)
        # print(comp_zones)
        return extracted, comp_zones

    elif profile_model is None:
        log.warning("Didn't receive identified targets "
                    "from {:s}".format(ccd.header['OFNAME']))
        raise NoTargetException
    else:
        log.error('Got wrong input')


def identify_targets(ccd, plots=False):
    """Identify targets cross correlating spatial profile with a gaussian model

    The method of cross-correlating a gaussian model to the spatial profile was
    mentioned in Marsh 1989, then I created my own implementation. The spatial
    profile is obtained by finding the median across the full dispersion axis.
    For goodman the spectrum and ccd are very well aligned, there is a slight
    variation from one configuration to another but in general is acceptable.

    Args:
        ccd (object): a ccdproc.CCDData instance

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
        # averaging overall spectral dimension because in goodman the spectra is
        # deviated very little
        profile_median = np.median(ccd.data, axis=1)

        # Gaussian has to be defined at the middle for it to work
        gaussian = models.Gaussian1D(amplitude=profile_median.max(),
                                     mean=profile_median.size // 2,
                                     stddev=order).rename('Gaussian')

        # do cross correlation of data with gaussian
        # this is useful for cases with low signal to noise ratio
        cross_corr = signal.correlate(in1=profile_median,
                                      in2=gaussian(range(profile_median.size)),
                                      mode='same')
        # filter cross correlation data
        filt_cross_corr = np.where(np.abs(cross_corr > cross_corr.min()
                                          + 0.03 * cross_corr.max()),
                                   cross_corr,
                                   None)
        peaks = signal.argrelmax(filt_cross_corr, axis=0, order=order)[0]

        profile_model = None
        for i in range(len(peaks)):
            low_lim = np.max([0, peaks[i] - 5])
            hig_lim = np.min([peaks[i] + 5, profile_median.size])
            amplitude = profile_median[low_lim: hig_lim].max()
            gaussian = models.Gaussian1D(amplitude=amplitude,
                                         mean=peaks[i],
                                         stddev=order).rename('Gaussian_'
                                                              '{:d}'.format(i))
            if profile_model is not None:
                profile_model += gaussian
            else:
                profile_model = gaussian
        # plt.axvline(peaks[i])
        # print(profile_model)

        # fit model to profile
        fit_gaussian = fitting.LevMarLSQFitter()
        profile_model = fit_gaussian(model=profile_model,
                                     x=range(profile_median.size),
                                     y=profile_median)
        # print(fit_gaussian.fit_info)
        if plots:
            # plt.plot(cross_corr, label='Cross Correlation', linestyle='--')
            plt.plot(profile_model(range(profile_median.size)),
                     label='Fitted Model')
            plt.plot(profile_median, label='Profile (median)', linestyle='--')
            plt.legend(loc='best')
            plt.show()

        # print(profile_model)
        if fit_gaussian.fit_info['ierr'] not in [1, 2, 3, 4]:
            log.warning('There is some problem with the fitting.'
                        'Returning None.')
            return None
        else:
            return profile_model

    else:
        log.error('Not a ccdproc.CCDData instance')
        return None


def trace_targets(ccd, profile, sampling_step=5, pol_deg=2, plots=True):
    """Find the trace of the target's spectrum on the image

    This function defines a low order polynomial that trace the location of the
    spectrum. The attributes pol_deg and sampling_step define the polynomial
    degree and the spacing in pixels for the samples. For every sample a
    gaussian model is fitted and the center (mean) is recorded and since
    spectrum traces vary smoothly this value is used as a new center for the
    base model used to fit the spectrum profile.

    Notes:
        This doesn't work for extended sources.

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

    # Get image dimensions
    spatial_length, dispersion_length = ccd.data.shape

    # Initialize model fitter
    model_fitter = fitting.LevMarLSQFitter()

    # Initialize the model to fit the traces
    trace_model = models.Polynomial1D(degree=pol_deg)

    # Will store the arrays for the fitted location of the target obtained
    # in the fitting
    trace_points = None

    # List that will contain all the Model instances corresponding to traced
    # targets
    all_traces = None

    # Array defining the points to be sampled
    sampling_axis = range(0,
                          dispersion_length // sampling_step * sampling_step,
                          sampling_step)

    # Loop to go through all the sampling points and gather the points
    for i in sampling_axis:
        # Fit the inital model to the data
        fitted_profile = model_fitter(model=profile,
                                      x=range(ccd.data[:, i].size),
                                      y=ccd.data[:, i])
        log.debug('Profile fitted')

        # if plots or log.isEnabledFor(logging.DEBUG):
        #     # print(fitted_profile)
        #     x_axis = range(ccd.data[:, i].size)
        #     plt.ion()
        #     plt.title('Column {:d}'.format(i))
        #     plt.plot(ccd.data[:, i], label='Data')
        #     plt.plot(x_axis, fitted_profile(x_axis), label='Model')
        #     plt.legend(loc='best')
        #     # if model_fitter.fit_info['ierr'] == 1:
        #     #     # print(model_fitter.fit_info)
        #     #     plt.ioff()
        #     #     plt.show()
        #     #     plt.ion()
        #     plt.draw()
        #     plt.pause(.5)
        #     plt.clf()
        #     plt.ioff()

        if model_fitter.fit_info['ierr'] not in [1, 2, 3, 4]:
            log.error(
                "Fitting did not work, fit_info['ierr'] = \
                {:d}".format(model_fitter.fit_info['ierr']))

        # alternatively could use fitted_profile.param_names
        # the mean_keys is another way to tell how many submodels are in the
        # model that was parsed.
        mean_keys = [key for key in dir(fitted_profile) if 'mean' in key]

        # Since traces are assumed to vary smoothly with wavelength I change the
        # previously found mean to the profile model used to fit the profile
        # data
        for mean_key in mean_keys:
            profile.__getattribute__(mean_key).value = \
                fitted_profile.__getattribute__(mean_key).value

        if trace_points is None:
            trace_points = np.ndarray((len(mean_keys),
                                       dispersion_length // sampling_step))

        # store the corresponding value in the proper array for later fitting
        # a low order polynomial
        for e in range(trace_points.shape[0]):
            log.debug('Median Trace'
                      ' {:d}: {:.3f}'.format(e, np.median(trace_points[e])))

            trace_points[e][i // sampling_step] = \
                fitted_profile.__getattribute__(mean_keys[e]).value

    # fit a low order polynomial for the trace
    for trace_num in range(trace_points.shape[0]):

        fitted_trace = model_fitter(model=trace_model,
                                    x=sampling_axis,
                                    y=trace_points[trace_num])

        fitted_trace.rename('Trace_{:d}'.format(trace_num))

        # TODO (simon): Validate which kind of errors are represented here
        if model_fitter.fit_info['ierr'] not in [1, 2, 3, 4]:
            log.error(model_fitter.fit_info['ierr'])
        else:
            # RMS Error calculation
            rms_error = np.sqrt(
                np.sum(np.array([fitted_trace(sampling_axis) -
                                 trace_points[trace_num]]) ** 2))
            log.info('Trace Fit RMSE: {:.3f}'.format(rms_error))

            # append for returning
            if all_traces is None:
                all_traces = [fitted_trace]
            else:
                all_traces.append(fitted_trace)

        if plots:
            plt.plot(sampling_axis,
                     trace_points[trace_num],
                     marker='o',
                     label='Data Points')
            plt.plot(fitted_trace(range(dispersion_length)),
                     label='Model RMSE: {:.2f}'.format(rms_error))

    if plots:
        plt.title('Targets Trace')
        plt.xlabel('Dispersion Axis')
        plt.ylabel('Spatial Axis')
        plt.imshow(ccd.data, cmap='YlGnBu', clim=(0, 300))

        for trace in trace_points:
            plt.plot(sampling_axis, trace, marker='.', linestyle='--')
            # print(trace)
        plt.legend(loc='best')
        plt.show()
    return all_traces


def get_extraction_zone(ccd,
                        extraction=None,
                        trace=None,
                        model=None,
                        n_sigma_extract=None,
                        zone=None,
                        plots=False):
    """Get the minimum rectangular CCD zone that fully contains the spectrum

    In this context, _fully contained_ means the spectrum plus some region for
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
        extraction (str): Extraction type, 'simple' or 'optimal'
        trace (object): An astropy.modeling.Model instance that correspond to
        the trace of the spectrum
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

    if zone is None and extraction is not None:
        assert (model is not None) and (n_sigma_extract is not None)
        assert isinstance(trace, Model)
        log.info('Extracting zone centered at: {:.3f}'.format(model.mean.value))
        spatial_length, dispersion_length = ccd.data.shape

        # get maximum variation in spatial direction
        trace_array = trace(range(dispersion_length))
        trace_inclination = trace_array.max() - trace_array.min()
        log.info('Trace Min-Max difference: {:.3f}'.format(trace_inclination))

        # m_mean = model.mean.value
        m_stddev = model.stddev.value
        extract_width = n_sigma_extract // 2 * m_stddev

        low_lim = np.max([0, int(trace_array.min() - extract_width)])

        hig_lim = np.min([int(trace_array.max() + extract_width),
                          spatial_length])

        zone = [low_lim, hig_lim]

        # this is necessary since we are cutting a piece of the full ccd.
        trace.c0.value -= low_lim
        log.info('Changing attribute c0 from trace, this is to adjust it to '
                 'the new extraction zone which is smaller that the full CCD.')

        log.info('Changing attribute mean of profile model')
        model.mean.value = extract_width

        nccd = ccd.copy()
        nccd.data = ccd.data[low_lim:hig_lim, :]
        if nccd.mask is not None:
            log.debug('Trimming mask')
            nccd.mask = ccd.mask[low_lim:hig_lim, :]
        nccd.header['HISTORY'] = 'Subsection of CCD ' \
                                 '[{:d}:{:d}, :]'.format(low_lim, hig_lim)

        if plots:
            plt.imshow(ccd.data)
            plt.axhspan(low_lim, hig_lim, color='r', alpha=0.2)
            plt.show()

        return nccd, trace, model, zone

    else:

        low_lim, hig_lim = zone

        nccd = ccd.copy()
        nccd.data = ccd.data[low_lim:hig_lim, :]
        if nccd.mask is not None:
            log.debug('Trimming mask')
            nccd.mask = ccd.mask[low_lim:hig_lim, :]
        nccd.header['HISTORY'] = 'Subsection of CCD ' \
                                 '[{:d}:{:d}, :]'.format(low_lim, hig_lim)
        return nccd


def add_wcs_keys(header):
    """Adds generic keyword to the header
    Linear wavelength solutions require a set of standard fits keywords. Later
    on they will be updated accordingly
    The main goal of putting them here is to have consistent and nicely ordered
    headers

    Args:
        header (object): New header without WCS entries

    Returns:
        header (object): Modified header

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
        log.error("Can't add wcs keywords to header")
        log.debug(err)


def remove_background(ccd, plots=False):
    """Remove Background of a ccd spectrum image

    Args:
        ccd (object): A ccdproc.CCDData instance.

    Returns:
        ccd (object): The modified

    """
    # ccd_copia = ccd.copy()
    data = ma.masked_invalid(ccd.data)
    # x, y = ccd.data.shape
    median = ma.median(data, axis=0)

    data -= median
    data.set_fill_value(-np.inf)
    ccd.data = data.filled()

    # ccd.write('/user/simon/dummy_{:d}.fits'.format(g), clobber=True)
    return ccd


def extract(ccd, trace, spatial_profile, extraction, sampling_step=1,
            plots=False):
    assert isinstance(ccd, CCDData)
    assert isinstance(trace, Model)

    spatial_length, dispersion_length = ccd.data.shape

    # create variance model
    rdnoise = float(ccd.header['RDNOISE'])
    gain = float(ccd.header['GAIN'])
    log.debug('Original Name {:s}'.format(ccd.header['OFNAME']))

    variance_2d = (rdnoise + np.absolute(ccd.data) * gain) / gain
    if ccd.mask is None and ccd.header['OBSTYPE'] == 'OBJECT':
        log.info('Finding cosmic rays to create mask')
        cr_mask = cosmicray_rejection(ccd=ccd, mask_only=True)
        cr_mask = np.logical_not(cr_mask).astype(int)
    elif ccd.mask is None and ccd.header['OBSTYPE'] != 'OBJECT':
        log.info('Only OBSTYPE == OBJECT get cosmic ray rejection.')
        cr_mask = np.ones(ccd.data.shape, dtype=int)
    else:
        log.debug('Cosmic ray mask already exists.')
        cr_mask = np.logical_not(ccd.mask).astype(int)

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
    log.debug('ccd.data is a masked array: '
              '{:s}'.format(str(np.ma.isMaskedArray(ccd.data))))

    ccd.data = np.ma.masked_invalid(ccd.data)
    # print(np.ma.isMaskedArray(ccd.data))
    np.ma.set_fill_value(ccd.data, 0)

    if extraction == 'simple':

        indexes = np.argwhere(cr_mask == 0)
        # print(indexes)
        if plots:
            fig = plt.figure(1)
            ax1 = fig.add_subplot(111)
            for index in indexes:
                x, y = index
                ax1.plot(y, x, marker='o', color='r')

            ax1.imshow(ccd.data, interpolation='none')
            if plt.isinteractive():
                plt.draw()
                plt.pause(1)
            else:
                plt.show()
        # print(np.ma.isMaskedArray(ccd.data))
        to_sum = ccd.data * cr_mask
        simple_sum = np.ma.sum(to_sum, axis=0)

        ccd.data = simple_sum

        if plots:
            fig = plt.figure()
            fig.canvas.set_window_title('Simple Extraction')
            # ax = fig.add_subplot(111)
            manager = plt.get_current_fig_manager()
            manager.window.maximize()

            plt.title(ccd.header['OBJECT'])
            plt.xlabel('Dispersion Axis (Pixels)')
            plt.ylabel('Intensity (Counts)')
            # plt.plot(simple_sum, label='Simple Sum', color='k', alpha=0.5)
            plt.plot(ccd.data, color='k',
                     label='Simple Extracted')
            plt.xlim((0, len(ccd.data)))
            plt.legend(loc='best')
            if plt.isinteractive():
                plt.draw()
                plt.pause(1)
            else:
                plt.show()

    elif extraction == 'optimal':
        out_spectrum = np.empty(dispersion_length)
        for i in range(0, dispersion_length, sampling_step):
            # force the model to follow the trace
            new_model.mean.value = trace(i)

            # warn if the difference of the spectrum position in the trace at the
            # extremes of the sampling range is larger than 1 pixel.
            if np.abs(trace(i) - trace(i + sampling_step)) > 1:
                log.warning('Sampling step might be too large')

            sample = np.median(ccd.data[:, i:i + sampling_step], axis=1)
            fitted_profile = model_fitter(model=new_model,
                                          x=range(len(sample)),
                                          y=sample)

            profile = fitted_profile(range(sample.size))

            # enforce positivity
            pos_profile = np.array([np.max([0, x]) for x in profile])

            # enforce normalization
            nor_profile = np.array([x / pos_profile.sum() for x in pos_profile])

            if sampling_step > 1:
                # TODO (simon): Simplify to Pythonic way
                right = min((i + sampling_step), dispersion_length)

                for e in range(i, right, 1):
                    mask = cr_mask[:, e]
                    data = ma.masked_invalid(ccd.data[:, e])
                    # print(ma.isMaskedArray(data))
                    V = variance_2d[:, e]
                    P = nor_profile
                    a = [(P[z] / V[z]) / np.sum(P ** 2 / V) for z in
                         range(P.size)]
                    weights = (nor_profile / variance_2d[:, e]) / np.sum(
                        nor_profile ** 2 / variance_2d[:, e])
                    # print('SUMN ', np.sum(a), np.sum(weights), np.sum(nor_profile), np.sum(P * weights))



                    # if e in range(5, 4001, 500):
                    #     plt.plot(nor_profile * data.max()/ nor_profile.max(), label=str(e))
                    #     plt.plot(data, label='Data')
                    #     plt.legend(loc='best')
                    #     plt.show()

                    out_spectrum[e] = np.sum(data * mask * nor_profile)
        ccd.data = out_spectrum
        if plots:
            fig = plt.figure()
            fig.canvas.set_window_title('Optimal Extraction')
            # ax = fig.add_subplot(111)
            manager = plt.get_current_fig_manager()
            manager.window.maximize()

            plt.title(ccd.header['OBJECT'])
            plt.xlabel('Dispersion Axis (Pixels)')
            plt.ylabel('Intensity (Counts)')
            # plt.plot(simple_sum, label='Simple Sum', color='k', alpha=0.5)
            plt.plot(ccd.data, color='k',
                     label='Optimal Extracted')
            plt.xlim((0, len(ccd.data)))
            plt.legend(loc='best')
            if plt.isinteractive():
                plt.draw()
                plt.pause(1)
            else:
                plt.show()

    # ccd.data = out_spectrum
    return ccd


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
            instrument (str): 'Red' or 'Blue' stating whether the data was taken
            using the Red or Blue Goodman Camera.
            technique (str): 'Spectroscopy' or 'Imaging' stating what kind of
            data was taken.
        """

        self.full_path = path
        self.instrument = instrument
        self.technique = technique
        self.is_empty = True
        self.bias = None
        self.day_flats = None
        self.dome_flats = None
        self.sky_flats = None
        self.data_groups = None
        self.spec_groups = None
        # time reference points
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

                log.error('Imaging mode needs BIAS to work properly. '
                          'Go find some.')

            else:
                log.warning('BIAS are needed for optimal results.')
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

    def add_spec_group(self, spec_group):
        """Adds a data group

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


class NoTargetException(Exception):
    def __init__(self):
        Exception.__init__(self, 'No targets identified.')
