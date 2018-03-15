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
import math
import pandas
import scipy
import shutil
import subprocess

from threading import Timer
# matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
from ccdproc import CCDData, ImageFileCollection
from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy.stats import sigma_clip
from astroplan import Observer
from astropy import units as u
from astropy.modeling import (models, fitting, Model)
from scipy import signal
from ..info import __version__

log = logging.getLogger(__name__)


def add_wcs_keys(ccd):
    """Adds generic keyword to the header

    Linear wavelength solutions require a set of standard fits keywords. Later
    on they will be updated accordingly
    The main goal of putting them here is to have consistent and nicely ordered
    headers

    Notes:
        This does NOT add a WCS solution, just the keywords

    Args:
        ccd (object): ccdproc.CCDData instance with no wcs keywords

    Returns:
        ccd (object): ccdproc.CCDData instance with modified header with added
          WCS keywords

    """
    try:
        log.debug("Adding FITS LINEAR wcs keywords to header.")
        ccd.header.set('BANDID1',
                       value='spectrum - background none, weights none, '
                             'clean no',
                       comment='')

        ccd.header.set('APNUM1',
                       value='1 1 0 0',
                       comment='')

        ccd.header.set('WCSDIM',
                       value=1,
                       comment='')

        ccd.header.set('CTYPE1',
                       value='LINEAR',
                       comment='')

        ccd.header.set('CRVAL1',
                       value=1,
                       comment='')

        ccd.header.set('CRPIX1',
                       value=1,
                       comment='')

        ccd.header.set('CDELT1',
                       value=1,
                       comment='')

        ccd.header.set('CD1_1',
                       value=1,
                       comment='')

        ccd.header.set('LTM1_1',
                       value=1,
                       comment='')

        ccd.header.set('WAT0_001',
                       value='system=equispec',
                       comment='')

        ccd.header.set('WAT1_001',
                       value='wtype=linear label=Wavelength units=angstroms',
                       comment='')

        ccd.header.set('DC-FLAG',
                       value=0,
                       comment='')

        ccd.header.set('DCLOG1',
                       value='REFSPEC1 = non set',
                       comment='')

        return ccd

    except TypeError as err:
        log.error("Can't add wcs keywords to header")
        log.debug(err)


def call_cosmic_rejection(ccd, image_name, out_prefix, red_path,
                          dcr_par, keep_files=False, prefix='c', method='dcr',
                          save=False):
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
        save (bool): Disables by default saving the images

    """
    if ccd.header['OBSTYPE'] == 'COMP':
        log.info("Setting cosmic ray rejection method to 'none' for "
                 "comparison lamp. 'c' prefix will be added anyways.")

        method = 'none'
        out_prefix = prefix + out_prefix

    if method == 'dcr':
        log.warning('DCR does apply the correction to images if you want '
                    'the mask use --keep-cosmic-files')

        full_path = os.path.join(red_path, out_prefix + image_name)

        ccd.header.set('GSP_COSM',
                       value="DCR",
                       comment="Cosmic ray rejection method")

        # ccd.write(full_path, clobber=True)
        write_fits(ccd=ccd, full_path=full_path)
        log.info('Saving image: {:s}'.format(full_path))

        in_file = out_prefix + image_name

        # This is to return the prefix that will be used by dcr
        # Not to be used by dcr_cosmicray_rejection
        out_prefix = prefix + out_prefix

        ccd = dcr_cosmicray_rejection(data_path=red_path,
                                      in_file=in_file,
                                      prefix=prefix,
                                      dcr_par_dir=dcr_par,
                                      delete=keep_files,
                                      save=save)
        return ccd, out_prefix

    elif method == 'lacosmic':
        log.warning('LACosmic does not apply the correction to images '
                    'instead it updates the mask attribute for CCDData '
                    'objects. For saved files the mask is a fits extension')

        ccd = lacosmic_cosmicray_rejection(ccd=ccd)

        out_prefix = prefix + out_prefix
        full_path = os.path.join(red_path, out_prefix + image_name)

        # ccd.write(full_path, clobber=True)
        log.info('Saving image: {:s}'.format(full_path))
        if save:
            write_fits(ccd=ccd, full_path=full_path)
        return ccd, out_prefix

    elif method == 'none':
        full_path = os.path.join(red_path, out_prefix + image_name)
        log.warning("--cosmic set to 'none'")
        # ccd.write(full_path, clobber=True)
        if save:
            write_fits(ccd=ccd, full_path=full_path)
        log.info('Saving image: {:s}'.format(full_path))

        return ccd, out_prefix

    else:
        log.error('Unrecognized Cosmic Method {:s}'.format(method))


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
    log.debug("Spectroscopic Data Classification")

    search_path = os.path.join(path, search_pattern + '*.fits')

    file_list = glob.glob(search_path)

    if file_list == []:
        log.error('No file found using search pattern '
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
            log.debug('Adding COMP group')
            data_container.add_comp_group(comp_group=spec_group)
        elif 'OBJECT' in group_obstype and len(group_obstype) == 1:
            log.debug('Adding OBJECT group')
            data_container.add_object_group(object_group=spec_group)
        else:
            log.debug('Adding OBJECT-COMP group')
            data_container.add_spec_group(spec_group=spec_group)

    # print(data_container)
    return data_container


def combine_data(image_list, dest_path, prefix=None, output_name=None,
                 method="median",
                 save=False):
    """Combine a list of CCDData instances.

    Args:
        image_list (list): Each element should an instance of CCDData
        dest_path (str): Path to where the new image should saved
        prefix (str): Prefix to add to the image file name
        output_name (str): Alternatively a file name can be parsed, this will
          ignore `prefix`.
        method (str): Method for doing the combination, this goes straight to
          the call of `ccdproc.combine` function
        save (bool): If True will save the combined images. If False it will
        ignore `prefix` or `output_name`.

    Returns:
        A combined image as a CCDData object.

    """
    assert len(image_list) > 1

    # This defines a default filename that should be deleted below
    combined_full_path = os.path.join(dest_path, "combined.fits")
    if output_name is not None:
        combined_full_path = os.path.join(dest_path, output_name)
    elif prefix is not None:
        combined_base_name = ''
        target_name = image_list[0].header["OBJECT"]

        grating_name = re.sub('[A-Za-z_-]',
                              '',
                              image_list[0].header["GRATING"])

        slit_size = re.sub('[A-Za-z" ]',
                           '',
                           image_list[0].header["SLIT"])

        for field in [prefix,
                      'combined',
                      target_name,
                      grating_name,
                      slit_size]:

            value = re.sub('[_ /]',
                           '',
                           field)

            combined_base_name += "{:s}_".format(value)

        # print(combined_base_name)

        number = len(glob.glob(
            os.path.join(dest_path,
                         combined_base_name + "*.fits")))

        combined_full_path = os.path.join(
            dest_path,
            combined_base_name + "{:02d}.fits".format(number + 1))

    # combine image
    combined_image = ccdproc.combine(image_list,
                                     method=method,
                                     sigma_clip=True,
                                     sigma_clip_low_thresh=1.0,
                                     sigma_clip_high_thresh=1.0,
                                     add_keyword=False)

    # add name of files used in the combination process
    for i in range(len(image_list)):
        image_name = image_list[i].header['GSP_FNAM']
        combined_image.header.set("GSP_IC{:02d}".format(i + 1),
                                  value=image_name,
                                  comment='Image used to create combined')

    if save:
        write_fits(combined_image,
                   full_path=combined_full_path,
                   combined=True)

    return combined_image


def convert_time(in_time):
    """Converts time to seconds since epoch

    Args:
        in_time (str): time obtained from header's keyword DATE-OBS

    Returns:
        time in seconds since epoch

    """
    return time.mktime(time.strptime(in_time, "%Y-%m-%dT%H:%M:%S.%f"))


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
        log.debug(profile_model.submodel_names)
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


def dcr_cosmicray_rejection(data_path, in_file, prefix, dcr_par_dir,
                            delete=False, save=True):
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
        save (bool):

    """

    log.info('Removing cosmic rays using DCR by Wojtek Pych')
    log.debug('See http://users.camk.edu.pl/pych/DCR/')

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
    # command = 'dcr {:s} {:s} {:s}'.format(in_file,
    #                                       out_file,
    #                                       cosmic_file)
    # print(command)

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

        log.debug('File dcr.par does not exist. Copying default one.')
        dcr_par_path = os.path.join(dcr_par_dir, 'dcr.par')
        log.debug('dcr.par full path: {:s}'.format(dcr_par_path))
        if os.path.isfile(dcr_par_path):
            shutil.copy2(dcr_par_path, data_path)
        else:
            log.error('Could not find dcr.par file')
    else:
        log.debug('File dcr.par exists.')

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
    def kill_process(process):
        log.error("DCR Timed out")
        process.kill()

    dcr_timer = Timer(10, kill_process, [dcr])
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
        log.error(stderr)
        if b'dcr: not found' in stderr:
            sys.exit('Your system can not locate the executable file dcr, try '
                     'moving it to /bin or create a symbolic link\n\n\tcd '
                     '/bin\n\tsudo ln -s /full/path/to/dcr')
    elif b'ERROR' in stdout:
        for output_line in stdout.split(b'\n'):
            log.error(output_line.decode("utf-8"))
    else:
        for output_line in stdout.split(b'\n'):
            log.debug(output_line)

    # delete extra files only if the execution ended without error
    if delete and stderr == b'' and b'USAGE:' not in stdout \
            and b'ERROR! calc_submean() failed' not in stdout:
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

    # recovers the saved file and returns the CCDData instance
    if os.path.isfile(full_path_out):
        ccd = CCDData.read(full_path_out, unit=u.adu)
        if not save:
            log.warning("Removing file because the attribute 'save' "
                        "is set to False")
            os.unlink(full_path_out)
        return ccd


def extraction(ccd,
               target_trace,
               spatial_profile,
               extraction_name):
    """Calls appropriate spectrum extraction routine

    This function calls the appropriate extraction function based on
    `extraction_name`

    Notes:
        Optimal extraction is not implemented.

    Args:
        ccd (object): Instance of ccdproc.CCDData containing a 2D spectrum
        target_trace (object): Instance of astropy.modeling.Model, a low order
            polynomial that defines the trace of the spectrum in the ccd object.
        spatial_profile (object): Instance of astropy.modeling.Model, a Gaussian
            model previously fitted to the spatial profile of the 2D spectrum
            contained in the ccd object.
        extraction_name (str): Extraction type, can be `fractional` or
            `optimal` though the optimal extraction is not implemented yet.

    Returns:
        ccd (object): Instance of ccdproc.CCDData containing a 1D spectrum. The
            attribute 'data' is replaced by the 1D array resulted from the
            extraction process.

    Raises:
        NotImplementedError: When `extraction_name` is `optimal`.

    """
    assert isinstance(ccd, CCDData)
    assert isinstance(target_trace, Model)

    if extraction_name == 'fractional':
        extracted, background = extract_fractional_pixel(
            ccd=ccd,
            target_trace=target_trace,
            target_stddev=spatial_profile.stddev.value,
            extraction_width=2)
        return extracted

    elif extraction_name == 'optimal':
        raise NotImplementedError


def extract_fractional_pixel(ccd, target_trace, target_stddev, extraction_width,
                             background_spacing=3):
    """Performs an spectrum extraction using fractional pixels.

    Args:
        ccd (object): Instance of ccdproc.CCDData that contains a 2D spectrum.
        target_trace (object):  Instance of astropy.modeling.models.Model that
          defines the trace of the target on the image (ccd).
        target_stddev (float): Standard deviation value for the spatial profile
          fitted to the target.
        extraction_width (int): Width of the extraction area as a function of
          `target_stddev`. For instance if `extraction_with` is set to 1 the
          function extract 0.5 to each side from the center of the traced
          target.
        background_spacing (float): Number of `target_stddev` to separate the
          target extraction to the background. This is from the edge of the
          extraction zone to the edge of the background region.
    """
    assert isinstance(ccd, CCDData)
    assert isinstance(target_trace, Model)
    log.info("Fractional Pixel Extraction for "
             "{:s}".format(ccd.header['GSP_FNAM']))

    spat_length, disp_length = ccd.data.shape

    disp_axis = range(disp_length)
    trace_points = target_trace(disp_axis)

    apnum1 = None  # '{:d} {:d} {:d} {:d}'.format(1, 1, 1, 1)

    non_background_sub = []
    extracted_spectrum = []
    background_list = []

    if ccd.header['OBSTYPE'] != 'OBJECT':
        log.debug("No background subtraction for OBSTYPE = "
                  "{:s}".format(ccd.header['OBSTYPE']))

    for i in disp_axis:

        # this defines the extraction limit for every column
        low_limit = trace_points[i] - 0.5 * extraction_width * target_stddev
        high_limit = trace_points[i] + 0.5 * extraction_width * target_stddev

        # low_limits_list.append(low_limit)
        # high_limits_list.append(high_limit)

        if apnum1 is None:
            # TODO (simon): add secondary targets
            apnum1 = '{:d} {:d} {:.2f} {:.2f}'.format(1,
                                                      1,
                                                      low_limit,
                                                      high_limit)
            ccd.header.set('APNUM1',
                           value=apnum1,
                           comment="Aperture in first column")

        column_sum = fractional_sum(data=ccd.data,
                                    index=i,
                                    low_limit=low_limit,
                                    high_limit=high_limit)

        non_background_sub.append(column_sum)

        if ccd.header['OBSTYPE'] == 'OBJECT':
            # background limits

            # background_spacing is the distance from the edge of the target's
            # limits defined by `int_low_limit` and
            # `int_high_limit` in stddev units

            background_width = high_limit - low_limit

            # define pixel values for background subtraction
            # low_background_zone
            high_1 = low_limit - background_spacing * target_stddev
            low_1 = high_1 - background_width
            # print(low_1,high_1)

            # High background zone
            low_2 = high_limit + background_spacing * target_stddev
            high_2 = low_2 + background_width
            # print(low_1,'-',high_1,':',low_2,'-',high_2,)

            # validate background subtraction zones
            background_1 = None
            background_2 = None

            # this has to be implemented, leaving it True assumes there is no
            # restriction for background selection.
            # TODO (simon): Implement background subtraction zones validation
            neighbouring_target_condition = True

            if low_1 > 0 and neighbouring_target_condition:
                # integer limits
                background_1 = fractional_sum(data=ccd.data,
                                              index=i,
                                              low_limit=low_1,
                                              high_limit=high_1)
            else:
                print("Invalid Zone 1")

            if high_2 < spat_length and neighbouring_target_condition:
                background_2 = fractional_sum(data=ccd.data,
                                              index=i,
                                              low_limit=low_2,
                                              high_limit=high_2)
            else:
                print("Invalid Zone 2")

            # background = 0
            if background_1 is not None and background_2 is None:
                background = background_1
            elif background_1 is None and background_2 is not None:
                background = background_2
            else:
                background = np.mean([background_1, background_2])

            # actual background subtraction
            background_subtracted_column_sum = column_sum - background

            # append column value to list
            extracted_spectrum.append(background_subtracted_column_sum)
            background_list.append(background)
        else:
            extracted_spectrum.append(column_sum)
    # plt.plot(low_limits_list, color='r')
    # plt.plot(high_limits_list, color='r')
    # plt.show()
    # plt.plot(non_background_sub, color='b', label='RAW')
    # plt.plot(extracted_spectrum, color='k', label='Extracted')
    # plt.plot(background_list, color='r', label='Background')
    # plt.legend(loc='best')
    # plt.show()

    new_ccd = ccd.copy()
    new_ccd.data = np.asarray(extracted_spectrum)
    return new_ccd, np.asarray(background_list)


def extract_optimal():
    raise NotImplementedError


def fix_duplicated_keywords(night_dir):
    """Remove duplicated keywords

    There are some cases when the raw data comes with duplicated keywords.
    The origin has not been tracked down. The solution is to identify the
    duplicated keywords and the remove all but one from the end backwards.

    Args:
        night_dir (str): The full path for the raw data location

    """
    log.debug('Finding duplicated keywords')
    log.warning('Files headers will be overwritten')
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
        log.debug('Found {:d} duplicated keyword '
                  '{:s}'.format(len(multiple_keys),
                                's' if len(multiple_keys) > 1 else ''))

        for image_file in files:
            log.debug('Processing Image File: {:s}'.format(image_file))
            try:
                ccd = CCDData.read(image_file, unit=u.adu)
                for keyword in multiple_keys:
                    while ccd.header.count(keyword) > 1:
                        ccd.header.remove(keyword,
                                          ccd.header.count(keyword) - 1)
                log.warning('Overwriting file with duplicated keywords removed')
                log.debug('File %s overwritten', image_file)
                ccd.write(image_file, clobber=True)
            except IOError as error:
                log.error(error)


def fractional_sum(data, index, low_limit, high_limit):
    """Performs a fractional pixels sum

    A fractional pixels sum is required several times while
    extracting a 1D spectrum from a 2D spectrum. The method
    is actually very simple.

    It requires the full data, the column and the range to sum, this
    range is given as real numbers. First it separates the limits values as an
    integer and fractional parts. Then it will sum the integer's interval and
    subtract the `low_limit`'s fractional part and sum the `high_limit`'s
    fractional part.

    The sum is performed in one operation. It does not do
    background subtraction, for which this very same method is used to
    get the background sum to be subtracted later.

    Args:
        data (numpy.ndarray): 2D array that contains the 2D spectrum/image
        index (int): Index of the column to be summed.
        low_limit (float): Lower limit for the range to be summed.
        high_limit (float): Higher limit for the range to be summed.

    Returns:
        Sum in ADU of all pixels and fractions between `low_limit` and
        `high_limit`.
    """
    # these are the limits within the full amount of flux on each pixel is
    # summed
    low_fraction, low_integer = math.modf(low_limit)
    high_fraction, high_integer = math.modf(high_limit)

    column_sum = np.sum(data[int(low_integer):int(high_integer), index]) - \
        data[int(low_integer), index] * low_fraction + \
        data[int(high_integer), index] * high_fraction

    # print(low_limit,
    #       high_limit,
    #       low_integer,
    #       high_integer,
    #       high_integer + 1,
    #       high_fraction,
    #       column_sum)

    return column_sum


def get_best_flat(flat_name):
    """Look for matching master flat

    Given a basename for master flats defined as a combination of key parameters
    extracted from the header of the image that we want to flat field, this
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
        log.debug('Found suitable master flat: {:s}'.format(master_flat_name))
        return master_flat, master_flat_name
    else:
        log.error('There is no flat available')
        return None, None


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

    log.debug('Found {:.3f} as central wavelength'.format(central_wavelength))

    return central_wavelength


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

    if fitted_box.width.value < x:
        log.debug("Slit detected. Adding a 10 pixels offset")
    else:
        log.debug("Slit limits not detected. Setting additional "
                  "offset to 0")
        offset = 0
    # Here we force the slit limits within the boundaries of the data (image)
    # this defines a preliminary set of slit limit
    l_lim = 1 + fitted_box.x_0.value - 0.5 * fitted_box.width.value + offset

    h_lim = 1 + fitted_box.x_0.value + 0.5 * fitted_box.width.value - offset

    low_lim = int(np.max([1 + offset, l_lim]))

    high_lim = int(np.min([h_lim, len(ccd_section_median) - offset]))

    # define the slit trim section as (IRAF)
    # convert o 1-based
    slit_trim_section = '[1:{:d},{:d}:{:d}]'.format(y,
                                                    low_lim,
                                                    high_lim)
    # print(slit_trim_section)

    # debugging plots that have to be manually turned on
    if False:
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
        plt.title('Slit Edge Detection')
        plt.plot(box_model(spatial_axis), color='c', label='Initial Box1D')
        plt.plot(fitted_box(spatial_axis), color='k', label='Fitted Box1D')
        plt.plot(ccd_section_median, label='Median Along Disp.')
        # plt.plot(pseudo_derivative, color='g', label='Pseudo Derivative')
        # plt.axvline(None, color='r', label='Detected Edges')
        # -1 to make it zero-based.
        plt.axvline(low_lim - 1, color='r', label='Detected Edges')
        plt.axvline(high_lim - 1, color='r')
        # for peak in peaks:
        #     plt.axvline(peak, color='r')
        plt.legend(loc='best')
        plt.show()

    # plt.imshow(master_flat.data[low_lim:high_lim, :])
    # plt.axvline(low_lim, color='r')
    # plt.axvline(high_lim, color='r')
    # plt.show()

    return slit_trim_section


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
                sun_set_time (str): Sun set time in the format
                  'YYYY-MM-DDTHH:MM:SS.SS'
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

            time_first_frame, time_last_frame = Time(min(date_obs)), Time(
                max(date_obs))

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

            return (twilight_evening,
                    twilight_morning,
                    sun_set_time,
                    sun_rise_time)


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
        new_x_axis = [i for i in range(len(clipped_profile))
                      if not clipped_profile.mask[i]]

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

        clipped_final_profile = clipped_final_profile[
            ~clipped_final_profile.mask]

        background_level = np.abs(np.max(clipped_final_profile) -
                                  np.min(clipped_final_profile))

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
            plt.plot(new_x_axis,
                     clipped_final_profile,
                     color='r',
                     label='Sigma Clip Data')

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
                log.debug('Discarding peak: {:.3f}'.format(val))

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

        profile_model = []
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
                profile_model.append(fitted_gaussian)
        if plots:
            plt.plot(median_profile, color='b')
            for profile in profile_model:
                plt.plot(profile(range(len(median_profile))), color='r')
            plt.show()

        # plt.imshow(ccd.data, clim=(50, 200), cmap='gray')
        # for peak in selected_peaks:
        #     plt.axhline(peak, color='r')
        # plt.show()

        if profile_model == []:
            log.error("Impossible to identify targets.")
            return profile_model
        else:
            return profile_model


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
    if overscan_region is not None:
        log.debug(
            'Applying overscan Correction: {:s}'.format(overscan_region))
        ccd = ccdproc.subtract_overscan(ccd=ccd,
                                        median=True,
                                        fits_section=overscan_region,
                                        add_keyword=add_keyword)

        ccd.header['GSP_OVER'] = (overscan_region, 'Overscan region')
    else:
        log.debug("Overscan region is None, returning the original data.")
        # ccd.header['GSP_OVER'] = ('none', 'Overscan region')

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
    if trim_section is not None:
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
            log.warning('Unrecognized trim type')
    else:
        log.info("{:s} trim section is not "
                 "defined.".format(trim_type.capitalize()))
        log.debug("Trim section is None, returning the same data.")

    return ccd


def interpolate(spectrum, interpolation_size):
    """Creates an interpolated version of the input spectrum

    This method creates an interpolated version of the input array, it is
    used mainly for a spectrum but it can also be used with any
    unidimensional array, assuming you are happy with the interpolation_size
    attribute defined for this class. The reason for doing interpolation is
    that it allows to find the lines and its respective center more
    precisely. The default interpolation size is 200 (two hundred) points.

    Args:
        spectrum (array): an uncalibrated spectrum or any unidimensional
          array.
        interpolation_size (int): Number of points to interpolate. (points added
          between two existing ones)
    Returns:
        Two dimensional array containing x-axis and interpolated array.
            The x-axis preserves original pixel values.


    """
    x_axis = range(spectrum.size)
    first_x = x_axis[0]
    last_x = x_axis[-1]

    new_x_axis = np.linspace(first_x,
                             last_x,
                             spectrum.size * interpolation_size)

    tck = scipy.interpolate.splrep(x_axis, spectrum, s=0)
    new_spectrum = scipy.interpolate.splev(new_x_axis, tck, der=0)
    return [new_x_axis, new_spectrum]


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

        ccd.header['GSP_COSM'] = ("LACosmic", "Cosmic ray rejection method")
        log.info("Cosmic rays rejected with LACosmic")
        if mask_only:
            return ccd.mask
        else:
            return ccd
    else:
        log.debug('Skipping cosmic ray rejection for image of datatype: '
                  '{:s}'.format(ccd.header['OBSTYPE']))
        return ccd


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
        log.debug('Normalizing by mean')
        master.data /= master.data.mean()

        master.header['GSP_NORM'] = ('mean', 'Flat normalization method')

    elif method == 'simple' or method == 'full':
        log.debug('Normalizing flat by {:s} model'.format(method))

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
            log.warning('This part of the code was left here for '
                        'experimental purposes only')
            log.warning('This procedure takes a lot to process, you might '
                        'want to see other method such as "simple" or '
                        '"mean".')
            for i in range(x_size):
                fit = model_fitter(model_init, x_axis, master.data[i])
                master.data[i] = master.data[i] / fit(x_axis)
            # master.header.add_history('Flat Normalized by full model')
            master.header['GSP_NORM'] = ('full', 'Flat normalization method')

    # write normalized flat to a file
    # master.write(norm_name, clobber=True)
    write_fits(ccd=master,
               full_path=norm_name,
               parent_file=name)

    return master, norm_name


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
                'destination': '-d --proc-path',
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
        log.debug('Default value for {:s} is {:s}'.format(
            arg_name[key],
            str(args.__getattribute__(key))))


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


def read_fits(full_path, technique='Unknown'):
    assert os.path.isfile(full_path)
    ccd = CCDData.read(full_path, unit=u.adu)

    all_keys = [key for key in ccd.header.keys()]

    ccd.header.set('GSP_VERS',
                   value=__version__,
                   comment='Goodman Spectroscopic Pipeline Version')

    if 'GSP_ONAM' not in all_keys:
        ccd.header.set('GSP_ONAM',
                       value=os.path.basename(full_path),
                       comment='Original file name')
    if 'GSP_ONAM' not in all_keys:
        ccd.header.set('GSP_PNAM',
                       value=os.path.basename(full_path),
                       comment='Parent file name')
        
    ccd.header.set('GSP_FNAM',
                   value=os.path.basename(full_path),
                   comment='Current file name')

    ccd.header.set('GSP_PATH',
                   value=os.path.dirname(full_path),
                   comment='Location at moment of reduce')

    if 'GSP_TECH' not in all_keys:
        ccd.header.set('GSP_TECH',
                       value=technique,
                       comment='Observing technique')

    if 'GSP_DATE' not in all_keys:
        ccd.header.set('GSP_DATE',
                       value=time.strftime("%Y-%m-%d"),
                       comment='Processing date')

    if 'GSP_OVER' not in all_keys:
        ccd.header.set('GSP_OVER',
                       value='none',
                       comment='Overscan region')

    if 'GSP_TRIM' not in all_keys:
        ccd.header.set('GSP_TRIM',
                       value='none',
                       comment='Trim section')

    if 'GSP_SLIT' not in all_keys:
        ccd.header.set('GSP_SLIT',
                       value='none',
                       comment='Slit trim section, slit illuminated area only')

    if 'GSP_BIAS' not in all_keys:
        ccd.header.set('GSP_BIAS',
                       value='none',
                       comment='Master bias image')

    if 'GSP_FLAT' not in all_keys:
        ccd.header.set('GSP_FLAT',
                       value='none',
                       comment='Master flat image')

    if 'GSP_NORM' not in all_keys:
        ccd.header.set('GSP_NORM',
                       value='none',
                       comment='Flat normalization method')

    if 'GSP_COSM' not in all_keys:
        ccd.header.set('GSP_COSM',
                       value='none',
                       comment='Cosmic ray rejection method')

    if 'GSP_WRMS' not in all_keys:
        ccd.header.set('GSP_WRMS',
                       value='none',
                       comment='Wavelength solution RMS Error')

    if 'GSP_WPOI' not in all_keys:
        ccd.header.set('GSP_WPOI',
                       value='none',
                       comment='Number of points used to '
                               'calculate wavelength solution')

    if 'GSP_WREJ' not in all_keys:
        ccd.header.set('GSP_WREJ',
                       value='none',
                       comment='Number of points rejected')
    if '' not in all_keys:
        ccd.header.add_blank('-- Goodman Spectroscopic Pipeline --',
                             before='GSP_VERS')
        
        ccd.header.add_blank('-- GSP END --', after='GSP_WREJ')

    ccd.header.set('BUNIT', after='CCDSUM')
    # ccd.header.set('', value='', comment='')
    # ccd.header.set('', value='', comment='')
    # ccd.header.set('', value='', comment='')
    # ccd.header.set('', value='', comment='')
    # ccd.header.set('', value='', comment='')
    # ccd.header.set('', value='', comment='')
    return ccd


def remove_background_by_median(ccd):
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


def save_extracted(ccd, destination):
    """Save extracted spectrum using `e` as prefix

    Args:
        ccd (object): CCDData instance
        destination (str): Path where the file will be saved.

    Returns:
        CCDData instance

    """
    assert isinstance(ccd, CCDData)
    assert os.path.isdir(destination)

    file_name = ccd.header['GSP_FNAM']
    new_file_name = 'e' + file_name
    log.info("Saving uncalibrated(w) extracted spectrum to file: "
             "{:s}".format(new_file_name))
    full_path = os.path.join(destination, new_file_name)
    write_fits(ccd=ccd, full_path=full_path, parent_file=file_name)
    return ccd


def search_comp_group(object_group, comp_groups, reference_data):
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
        reference_data (object): Instance of
          `goodman.pipeline.core.ReferenceData` contains all information
          related to the reference lamp library.

    Returns:

    """
    log.debug('Finding a suitable comparison lamp group')

    object_confs = object_group.groupby(['grating',
                                         'cam_targ',
                                         'grt_targ',
                                         'filter',
                                         'filter2']
                                        ).size().reset_index()
    # .rename(columns={0: 'count'})

    for comp_group in comp_groups:

        if ((comp_group['grating'] == object_confs.iloc[0]['grating']) &
                (comp_group['cam_targ'] == object_confs.iloc[0]['cam_targ']) &
                (comp_group['grt_targ'] == object_confs.iloc[0]['grt_targ']) &
                (comp_group['filter'] == object_confs.iloc[0]['filter']) &
                (comp_group['filter2'] == object_confs.iloc[0]['filter2']
                 )).all():
            if reference_data.check_comp_group(comp_group) is not None:
                # print(comp_group)
                log.debug('Found a matching comparison lamp group')
                return comp_group

    raise NoMatchFound


def trace(ccd, model, trace_model, model_fitter, sampling_step, nsigmas=2):
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
        model_fitter (object): An astropy.modeling.fitting.Fitter instance. Will
          fit the sampled points to construct the trace model
        sampling_step (int): Step for sampling the spectrum.
        nsigmas (int): Number of stddev to each side of the mean to be used for
          searching the trace.

    Returns:
        An astropy.modeling.Model instance, that defines the trace of the
          spectrum.

    """
    # TODO (simon): modify this function in a way that does not requiere to
    # TODO (simon): parse objects other than ccd. Models should be instantiated
    # TODO (simon): (created) inside this function.
    spatial_length, dispersion_length = ccd.data.shape

    sampling_axis = range(0, dispersion_length, sampling_step)
    sample_values = []

    model_stddev = model.stddev.value
    model_mean = model.mean.value

    sample_center = float(model_mean)
    lower_limit = None
    upper_limit = None

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

        if np.abs(sample_peak + lower_limit - model_mean) < \
                nsigmas * model_stddev:
            sample_center = int(sample_peak + lower_limit)
        else:
            # print(np.abs(sample_peak + lower_limit - model_mean),
            #       nsigmas * model_stddev)
            sample_center = float(model_mean)

    fitted_trace = model_fitter(trace_model, sampling_axis, sample_values)

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


def trace_targets(ccd, target_list, sampling_step=5, pol_deg=2, plots=False):
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
        target_list (list): List of single target profiles.
        sampling_step (int): Frequency of sampling in pixels
        pol_deg (int): Polynomial degree for fitting the trace
        plots (bool): If True will show plots (debugging)

    Returns:
        all_traces (list): List that contains traces that are
            astropy.modeling.Model instance

    """

    # added two assert for debugging purposes
    assert isinstance(ccd, CCDData)
    assert all([isinstance(profile, Model) for profile in target_list])

    # Initialize model fitter
    model_fitter = fitting.LevMarLSQFitter()

    # Initialize the model to fit the traces
    trace_model = models.Polynomial1D(degree=pol_deg)

    # List that will contain all the Model instances corresponding to traced
    # targets
    all_traces = []

    for profile in target_list:
        single_trace = trace(ccd=ccd,
                             model=profile,
                             trace_model=trace_model,
                             model_fitter=model_fitter,
                             sampling_step=sampling_step,
                             nsigmas=10)
        if 0 < single_trace.c0.value < ccd.shape[0]:
            log.debug('Adding trace to list')
            all_traces.append([single_trace, profile])
        else:
            # print(profile)
            # plt.plot(profile(range(ccd.shape[0])))
            # plt.plot(np.median(ccd.data, axis=1))
            # plt.show()
            # print(single_trace.c0, ccd.shape[0])
            log.error("Unable to trace target.")
            log.error('Trace is out of boundaries. Center: '
                      '{:.4f}'.format(single_trace.c0.value))
        # print(single_trace)
    if plots:
        plt.title('Traces')
        plt.imshow(ccd.data)
        for strace, prof in all_traces:
            plt.plot(strace(range(ccd.data.shape[1])), color='r')
        plt.show()
    return all_traces


def write_fits(ccd,
               full_path,
               combined=False,
               parent_file=None,
               overwrite=True):
    assert isinstance(ccd, CCDData)
    assert os.path.isdir(os.path.dirname(full_path))

    # Original File Name
    # This should be set only once.
    if combined:
        ccd.header.set('GSP_ONAM',
                       value=os.path.basename(full_path))

        ccd.header.set('GSP_PNAM',
                       value='combined')

    # Parent File Name
    if not combined and parent_file is not None:
        ccd.header.set('GSP_PNAM',
                       value=os.path.basename(parent_file))

    # Current File Name
    ccd.header.set('GSP_FNAM', value=os.path.basename(full_path))

    # write to file
    ccd.write(full_path, overwrite=overwrite)
    assert os.path.isfile(full_path)
    return ccd


# classes definition


class CriticalError(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)


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
        self.gain = None
        self.rdnoise = None
        self.roi = None
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

    def __repr__(self):
        if self.is_empty:
            return str("Empty Data Container")
        else:

            class_info = str("{:s}\n"
                             "Full Path: {:s}\n"
                             "Instrument: {:s}\n"
                             "Technique: {:s}".format(str(self.__class__),
                                                      self.full_path,
                                                      self.instrument,
                                                      self.technique))
            # print(class_info)

            if all([self.gain, self.rdnoise, self.roi]):
                class_info += str("\nGain: {:.2f}\n"
                                  "Readout Noise: {:.2f}\n"
                                  "ROI: {:s}".format(self.gain,
                                                     self.rdnoise,
                                                     self.roi))
            class_info += str("\nIs Empty: {:s}\n".format(str(self.is_empty)))

            group_info = "\n Data Grouping Information\n"
            if self.technique == 'Spectroscopy':

                group_info += "COMP Group:\n"
                group_info += self._get_group_repr(self.comp_groups)
                group_info += "OBJECT Group\n"
                group_info += self._get_group_repr(self.object_groups)
                group_info += "OBJECT + COMP Group:\n"
                group_info += self._get_group_repr(self.spec_groups)

                class_info += group_info
            return class_info

    @staticmethod
    def _get_group_repr(group):
        """Converts the file names in each group to string

        This class has a __repr__ method and in this method the file names
        contained in the different groups gets formatted as a string for
        displaying in a readable way.
        """
        group_str = ""
        if group is not None:
            for i in range(len(group)):
                group_str += " Group {:d}\n".format(i + 1)
                for _file in group[i]['file']:
                    group_str += "  {:s}\n".format(_file)
            return group_str
        else:

            return "  Group is Empty\n"

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
        comp_group = spec_group[spec_group.obstype == 'COMP']
        self.add_comp_group(comp_group=comp_group)
        # print(comp_group)

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

    def set_readout(self, gain, rdnoise, roi):
        self.gain = gain
        self.rdnoise = rdnoise
        self.roi = roi

    # def reset(self):
    #     """Resets the class as it was first initialized"""
    #     self.__init__(path=self.full_path,
    #                   instrument=self.instrument,
    #                   technique=self.technique)


class NoMatchFound(Exception):
    def __init__(self):
        Exception.__init__(self, 'Did not find a match')


class NoTargetException(Exception):
    """Exception to be raised when no target is identified"""
    def __init__(self):
        Exception.__init__(self, 'No targets identified.')


class NotEnoughLinesDetected(Exception):
    def __init__(self):
        Exception.__init__(self, 'Not enough lines detected.')


class ReferenceData(object):
    """Contains spectroscopic reference lines values and filename to templates.

    This class stores:
        - file names for reference fits spectrum
        - file names for CSV tables with reference lines and relative
          intensities
        - line positions only for the elements used in SOAR comparison lamps
    """

    def __init__(self, reference_dir):
        """Init method for the ReferenceData class

        This methods uses ccdproc.ImageFileCollection on the reference_dir to
        capture all possible reference lamps. Also defines dictionaries
        containg line lists for several elements used in lamps.

        Args:
            reference_dir (str): full path to the reference data directory
        """
        self.log = logging.getLogger(__name__)
        self.reference_dir = reference_dir
        reference_collection = ccdproc.ImageFileCollection(self.reference_dir)
        self.ref_lamp_collection = reference_collection.summary.to_pandas()
        # print(self.ref_lamp_collection)
        self.lines_pixel = None
        self.lines_angstrom = None
        self._ccd = None

    def get_reference_lamp(self, header):
        """Finds a suitable template lamp from the catalog

        Args:
            header (object): FITS header of image we are looking a a reference
                lamp.

        Returns:
            full path to best matching reference lamp.

        """
        filtered_collection = self.ref_lamp_collection[
            # TODO (simon): User input for OBJECT keyword can differ from
            # TODO Reference Library
            (self.ref_lamp_collection['object'] == header['object']) &
            # TODO (simon): Wavemode can be custom (GRT_TARG, CAM_TARG, GRATING)
            (self.ref_lamp_collection['wavmode'] == header['wavmode'])]

        # print(filtered_collection)
        if filtered_collection.empty:
            raise NotImplementedError("It was not possible to find any lamps "
                                      "that match")
        elif len(filtered_collection) == 1:
            self.log.info(
                "Lamp: {:s}".format(filtered_collection.file.to_string(
                    index=False)))
            full_path = os.path.join(self.reference_dir,
                                     filtered_collection.file.to_string(
                                         index=False))
            self._ccd = CCDData.read(full_path, unit=u.adu)
            self._recover_lines()
            return self._ccd
        else:
            raise NotImplementedError

    def lamp_exists(self, object_name, grating, grt_targ, cam_targ):
        """

        Args:
            object_name:
            grating:
            grt_targ:
            cam_targ:

        Returns:

        """
        filtered_collection = self.ref_lamp_collection[
            (self.ref_lamp_collection['object'] == object_name) &
            (self.ref_lamp_collection['grating'] == grating) &
            (self.ref_lamp_collection['grt_targ'] == grt_targ) &
            (self.ref_lamp_collection['cam_targ'] == cam_targ)
            ]

        if filtered_collection.empty:
            return False
        elif len(filtered_collection) == 1:
            return True
        else:
            raise NotImplementedError

    def check_comp_group(self, comp_group):
        """

        Args:
            comp_group:

        Returns:

        """
        lamps = comp_group.groupby(['object',
                                    'grating',
                                    'grt_targ',
                                    'cam_targ']).size().reset_index(
        ).rename(columns={0: 'count'})

        for i in lamps.index:
            if self.lamp_exists(
                    object_name=lamps.iloc[i]['object'],
                    grating=lamps.iloc[i]['grating'],
                    grt_targ=lamps.iloc[i]['grt_targ'],
                    cam_targ=lamps.iloc[i]['cam_targ']):
                new_group = comp_group[
                    (comp_group['object'] == lamps.iloc[i]['object']) &
                    (comp_group['grating'] == lamps.iloc[i]['grating']) &
                    (comp_group['grt_targ'] == lamps.iloc[i]['grt_targ']) &
                    (comp_group['cam_targ'] == lamps.iloc[i]['cam_targ'])]
                # print(new_group.file)
                return new_group
            else:
                # print(lamps.iloc[i])
                # print(comp_group.file)
                self.log.warning("The target's comparison lamps do not have "
                                 "reference lamps.")
                self.log.debug("In this case a compatible lamp will be "
                               "obtained from all the lamps obtained in the "
                               "data or present in the files.")
        return None

    def _recover_lines(self):
        self.log.info("Recovering line information from reference Lamp.")
        self.lines_pixel = []
        self.lines_angstrom = []
        pixel_keys = self._ccd.header['GSP_P*']
        for pixel_key in pixel_keys:
            if re.match(r'GSP_P\d{3}', pixel_key) is not None:
                angstrom_key = re.sub('GSP_P', 'GSP_A', pixel_key)
                assert pixel_key[-3:] == angstrom_key[-3:]
                assert angstrom_key in self._ccd.header
                if int(self._ccd.header[angstrom_key]) != 0:
                    self.lines_pixel.append(float(self._ccd.header[pixel_key]))
                    self.lines_angstrom.append(
                        float(self._ccd.header[angstrom_key]))
                else:
                    self.log.debug(
                        "File: {:s}".format(self._ccd.header['GSP_FNAM']))
                    self.log.warning(
                        "Ignoring keywords: {:s}={:f}, {:s}={:f}".format(
                            pixel_key,
                            self._ccd.header[pixel_key],
                            angstrom_key,
                            self._ccd.header[angstrom_key]))

    def _validate_lines(self):
        """Calls all available validation methods

        Returns:
            True if none of the validation fails.
        """
        assert len(self.lines_pixel) == len(self.lines_angstrom)
        if not self._order_validation(self.lines_pixel):
            return False
        if not self._order_validation(self.lines_angstrom):
            return False
        # if not self._validate_line_existence():
        #     return False
        self._validate_line_existence()
        return True

    @staticmethod
    def _order_validation(lines_array):
        """Checks that the array of lines only increases."""
        previous = None
        for line_value in lines_array:
            # print(line_value)
            if previous is not None:
                try:
                    assert line_value > previous
                    previous = line_value
                except AssertionError:
                    print("Error: Line {:f} is not larger "
                          "than {:f}".format(line_value, previous))
                    return False
            else:
                previous = line_value
        return True

    def _load_nist_list(self, **kwargs):
        """Load all csv files from strong lines in nist."""
        nist_path = kwargs.get(
            'path',
            os.path.join(os.path.dirname(sys.modules['pipeline'].__file__),
                         'data/nist_list'))
        assert os.path.isdir(nist_path)
        nist_files = glob.glob(os.path.join(nist_path, "*.txt"))
        for nist_file in nist_files:
            key = os.path.basename(nist_file)[22:-4]
            nist_data = pandas.read_csv(nist_file, names=['intensity',
                                                          'air_wavelength',
                                                          'spectrum',
                                                          'reference'])
            self.nist[key] = nist_data

    def _validate_line_existence(self):
        """Check if a line actually exists in any table of NIST
        Notes:
            It does not actually check NIST, it loads six csv tables
            from NIST's strong lines for the elements used in lamps.
            Ar,  Cu, Fe, He, Hg, Ne. It does not work perfect so far
            so the it is not checking existence actually but if it finds it
            it will get the value at "spectrum" column in NIST tables which
            correspond to the source of the line for instance Ne I, Ar II.
            """

        lamp_elements = []
        lamp_name = self._ccd.header['OBJECT']
        if len(lamp_name) % 2 == 0:
            for element_index in range(0, len(lamp_name), 2):
                element = lamp_name[element_index:element_index + 2]
                lamp_elements.append(element)

        if self.nist is None:
            self.nist = {}
            self._load_nist_list()

        self.spectrum = list(self.lines_angstrom)
        for i in range(len(self.lines_angstrom)):
            for element in lamp_elements:
                line_info = self.nist[element][
                    self.nist[element].air_wavelength == self.lines_angstrom[i]]
                if line_info.empty:
                    # print(self.lines_angstrom[i], 'no-info')
                    self.spectrum[i] = ''
                else:
                    self.spectrum[i] = line_info['spectrum'].to_string(
                        index=False)


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
                 (self.modes_data_frame['grttarg'] == grating_targ) &
                 (self.modes_data_frame['ob_filter'] == blocking_filter))]
            if _mode.empty:
                central_wavelength = get_central_wavelength(
                    grating=grating,
                    grt_ang=grating_targ,
                    cam_ang=camera_targ)
                return 'Custom_{:d}nm'.format(int(round(central_wavelength)))
            else:
                # print('%s %s' %(grating, _mode.wavmode))
                # print(_mode['wavmode'].to_string)
                return _mode['wavmode'].to_string(index=False)
