from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import calendar
import ccdproc
import collections
import datetime
import glob
import logging
import math
import numpy as np
import os
import pandas
import random
import re
import scipy
import subprocess
import sys
import time

from astropy.utils import iers
iers.Conf.iers_auto_url.set('ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all')
from astroplan import Observer
from astropy import units as u
from astropy.io import fits
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
from astropy.coordinates import EarthLocation
from astropy.modeling import (models, fitting, Model)
from astropy.stats import sigma_clip
from astropy.time import Time
from astroscrappy import detect_cosmics
from ccdproc import CCDData, ImageFileCollection
from matplotlib import pyplot as plt
from scipy import signal
from threading import Timer

from . import check_version

__version__ = __import__('goodman_pipeline').__version__

log = logging.getLogger(__name__)


def astroscrappy_lacosmic(ccd, red_path=None, save_mask=False):

    mask, ccd.data = detect_cosmics(ccd.data)

    ccd.header['GSP_COSM'] = ('LACosmic',
                              "Cosmic ray rejection method")
    log.info("Cosmic rays rejected using astroscrappy's lacosmic")

    if save_mask and red_path is not None:
        mask_ccd = ccd.copy()
        mask_ccd.mask = mask
        new_file_name = 'crmask_' + mask_ccd.header['GSP_FNAM']
        mask_ccd.header['GSP_FNAM'] = new_file_name
        log.info("Saving binary mask of cosmic rays to "
                 "{:s}".format(new_file_name))
        write_fits(ccd=mask_ccd,
                   full_path=os.path.join(red_path, new_file_name))

    return ccd


def add_wcs_keys(ccd):
    """Adds generic keyword for linear wavelength solution to the header

    Linear wavelength solutions require a set of standard fits keywords. Later
    on they will be updated accordingly.
    The main goal of putting them here is to have consistent and nicely ordered
    headers.

    Notes:
        This does NOT add a WCS solution, just the keywords.

    Args:
        ccd (CCDData) A :class:~astropy.nddata.CCDData` instance with no wcs
          keywords.

    Returns:
        ccd (CCDData) A :class:`~astropy.nddata.CCDData` instance with modified
          header with added WCS keywords

    """
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


def add_linear_wavelength_solution(ccd, x_axis, reference_lamp, crpix=1):
    """Add wavelength solution to the new FITS header

    Defines FITS header keyword values that will represent the wavelength
    solution in the header so that the image can be read in any other
    astronomical tool. (e.g. IRAF)

    Args:
        ccd (CCDData) Instance of :class:`~astropy.nddata.CCDData`
        x_axis (ndarray): Linearized x-axis in angstrom
        reference_lamp (str): Name of lamp used to get wavelength solution.
        crpix (int): reference pixel for defining wavelength solution.
        Default 1. For most cases 1 should be fine.


    Returns:
        ccd (CCDData) A :class:`~astropy.nddata.CCDData` instance with
          linear wavelength solution on it.

    """
    assert crpix > 0
    new_crpix = crpix
    new_crval = x_axis[new_crpix - crpix]
    new_cdelt = x_axis[new_crpix] - x_axis[new_crpix - crpix]

    ccd.header.set('BANDID1', 'spectrum - background none, weights none, '
                              'clean no')
    ccd.header.set('WCSDIM', 1)
    ccd.header.set('CTYPE1', 'LINEAR  ')
    ccd.header.set('CRVAL1', new_crval)
    ccd.header.set('CRPIX1', new_crpix)
    ccd.header.set('CDELT1', new_cdelt)
    ccd.header.set('CD1_1', new_cdelt)
    ccd.header.set('LTM1_1', 1.)
    ccd.header.set('WAT0_001', 'system=equispec')
    ccd.header.set('WAT1_001', 'wtype=linear label=Wavelength units=angstroms')
    ccd.header.set('DC-FLAG', 0)
    ccd.header.set('DCLOG1', 'REFSPEC1 = {:s}'.format(reference_lamp))

    return ccd


def bias_subtract(ccd, master_bias, master_bias_name):
    """Subtract bias from file.

    Wrapper for :func:`~ccdproc.subtract_bias`. The main goal is to have a
    consistent API for apps using the Goodman Pipeline as a library.

    Args:
        ccd (CCDData): A file to be bias-subtracted
        master_bias (CCDData):
        master_bias_name (str): Full path to master bias file, this is added to
        the bias-subtracted ccd under `GSP_BIAS`.

    Returns:
        A bias-subtracted file.
    """
    ccd = ccdproc.subtract_bias(ccd=ccd, master=master_bias, add_keyword=False)
    log.info("Bias subtracted")
    ccd.header.set('GSP_BIAS',
                   value=os.path.basename(master_bias_name),
                   comment="Master Bias Image")
    return ccd


def bin_reference_data(wavelength, intensity, serial_binning):
    """Bins a 1D array

    This method reduces the size of an unbinned array by binning.
    The function to combine data is `numpy.mean`.

    Args:
        wavelength (array): Wavelength axis
        intensity (array): Intensity
        serial_binning (int): Serial Binning is the binning in the
        dispersion axis.

    Returns:
        Binned wavelength and intensity arrays.

    """
    if serial_binning != 1:
        b_wavelength = ccdproc.block_reduce(wavelength,
                                            serial_binning,
                                            np.mean)
        b_intensity = ccdproc.block_reduce(intensity,
                                           serial_binning,
                                           np.mean)
        return b_wavelength, b_intensity
    else:
        return wavelength, intensity


def call_cosmic_rejection(ccd,
                          image_name,
                          out_prefix,
                          red_path,
                          keep_files=False,
                          prefix='c',
                          method='dcr',
                          save=False):
    """Call for the appropriate cosmic ray rejection method

    There are four options when dealing with cosmic ray rejection in this
    pipeline, The default option is called ``default`` and it will choose the
    rejection method based on the binning of the image. Note that there are only
    two *real* methdos: ``dcr`` and ``lacosmic``.

    For ``binning 1x1`` the choice will be ``dcr`` for ``binning 2x2`` and
    ``3x3`` will be ``lacosmic``.

    The method ``dcr`` is a program written in C by Wojtek
    Pych (http://users.camk.edu.pl/pych/DCR/) that works very well for
    spectroscopy the only negative aspect is that integration with python was
    difficult and not natively (through subprocess).

    The method `lacosmic` is well known but there are different implementations,
    we started using :func:`~ccdproc.cosmicray_lacosmic` but later we shifted
    towards ``astroscrappy.detect_cosmics``. The LACosmic method was developed
    by Pieter G. van Dokkum. See <http://www.astro.yale.edu/dokkum/lacosmic/>

    There is also the option of skipping cosmic ray removal by using ``none``.

    Args:
        ccd (CCCData): a :class:`~astropy.nddata.CCDData` instance.
        image_name (str): Science image name.
        out_prefix (str): Partial prefix to be added to the image name. Related
          to previous processes and not cosmic ray rejection.
        red_path (str): Path to reduced data directory.
        keep_files (bool): If True, the original file and the cosmic ray mask
          will not be deleted. Default is False.
        prefix (str): Cosmic ray rejection related prefix to be added to image
          name.
        method (str): Method to use for cosmic ray rejection. There are four
          options: `default`, `dcr`, `lacosmic` and `none`.
        save (bool): Disables by default saving the images

    Returns:
        :class:`~astropy.nddata.CCDData` instance and `out_prefix` which is the
          prefix added to the image name.

    Raises:
        NotImplementedError if the `method` argument is not `dcr`, `lacosmic`
        nor `none`.

    """
    log.debug("Cosmic ray rejection method from input is '{:s}'".format(method))

    binning, _ = [int(i) for i in ccd.header['CCDSUM'].split()]

    if method == 'default':
        if binning == 1:
            method = 'dcr'
            log.info('Setting cosmic ray rejection method to:'
                     ' {:s}'.format(method))
        elif binning == 2:
            method = 'lacosmic'
            log.info('Setting cosmic ray rejection method to:'
                     ' {:s}'.format(method))
        elif binning == 3:
            method = 'lacosmic'
            log.info('Setting cosmic ray rejection method to:'
                     ' {:s}'.format(method))

    if ccd.header['OBSTYPE'] == 'COMP' and method != 'none':
        log.info("Changing cosmic ray rejection method from '{:s}' to 'none'"
                 " for comparison lamp. Prefix 'c' will be added "
                 "anyway.".format(method))

        method = 'none'
        log.debug("Cosmic ray rejection changed to 'none' for this file: "
                  "{:s}".format(ccd.header['GSP_FNAM']))
        out_prefix = prefix + out_prefix

    if method == 'dcr':
        log.warning('DCR does apply the correction to images if you want '
                    'the mask use --keep-cosmic-files')

        if not os.path.isfile(os.path.join(red_path, 'dcr.par')):
            _create = GenerateDcrParFile()
            _instrument = ccd.header['INSTCONF']
            _binning, _ = ccd.header['CCDSUM'].split()

            _create(instrument=_instrument, binning=_binning, path=red_path)

        full_path = os.path.join(red_path, out_prefix + image_name)

        ccd.header.set('GSP_COSM',
                       value="DCR",
                       comment="Cosmic ray rejection method")

        write_fits(ccd=ccd, full_path=full_path)
        log.info('Saving image: {:s}'.format(full_path))

        in_file = out_prefix + image_name

        # This is to return the prefix that will be used by dcr
        # Not to be used by dcr_cosmicray_rejection
        out_prefix = prefix + out_prefix

        ccd = dcr_cosmicray_rejection(data_path=red_path,
                                      in_file=in_file,
                                      prefix=prefix,
                                      keep_cosmic_files=keep_files,
                                      save=save)
        return ccd, out_prefix

    elif method == 'lacosmic':
        ccd = astroscrappy_lacosmic(ccd=ccd,
                                    red_path=red_path,
                                    save_mask=keep_files)

        out_prefix = prefix + out_prefix
        full_path = os.path.join(red_path, out_prefix + image_name)

        if save:
            log.info('Saving image: {:s}'.format(full_path))
            write_fits(ccd=ccd, full_path=full_path)
        return ccd, out_prefix

    elif method == 'none':
        full_path = os.path.join(red_path, out_prefix + image_name)
        if save:
            log.info('Saving image: {:s}'.format(full_path))
            write_fits(ccd=ccd, full_path=full_path)

        return ccd, out_prefix

    else:
        log.error('Unrecognized Cosmic Method {:s}'.format(method))
        raise NotImplementedError


def create_master_bias(bias_files,
                       raw_data,
                       reduced_data,
                       technique):
    """Create Master Bias

    Given a :class:`~pandas.DataFrame` object that contains a list of compatible bias.
    This function creates the master flat using ccdproc.combine using median
    and 3-sigma clipping.

    Args:
        bias_files (list): List of all bias files to be combined. They have
        to be compatible with each other as no check is done in this method.
        raw_data (str): Full path to raw data location.
        reduced_data (str): Full path to were reduced data will reside.
        technique (str): Name of observing technique. Imaging or
        Spectroscopy.

    Returns:
        master_bias (object):
        master_bias_name (str):

    """
    assert isinstance(bias_files, list)

    master_bias_list = []
    log.info('Creating master bias')
    for image_file in bias_files:
        image_full_path = os.path.join(raw_data, image_file)
        ccd = read_fits(image_full_path, technique=technique)
        log.debug('Loading bias image: ' + image_full_path)
        master_bias_list.append(ccd)

    # combine bias for spectroscopy
    log.info("Combining {} images to create master bias".format(
        len(master_bias_list)))
    master_bias = ccdproc.combine(master_bias_list,
                                  method='median',
                                  sigma_clip=True,
                                  sigma_clip_low_thresh=3.0,
                                  sigma_clip_high_thresh=3.0,
                                  add_keyword=False)

    # add name of images used to create master bias
    for n in range(len(bias_files)):
        master_bias.header['GSP_IC{:02d}'.format(n + 1)] = (
            bias_files[n],
            'Image used to create master bias')

    master_bias_name = "master_bias_{}_{}_{}_R{:05.2f}_G{:05.2f}.fits".format(
        master_bias.header['INSTCONF'].upper(),
        technique[0:2].upper(),
        "x".join(master_bias.header['CCDSUM'].split()),
        float(master_bias.header['RDNOISE']),
        float(master_bias.header['GAIN'])
    )

    write_fits(ccd=master_bias,
               full_path=os.path.join(reduced_data, master_bias_name),
               combined=True,
               overwrite=True)

    log.info('Created master bias: ' + master_bias_name)
    return master_bias, master_bias_name


def create_master_flats(flat_files,
                        raw_data,
                        reduced_data,
                        technique,
                        overscan_region,
                        trim_section,
                        master_bias_name,
                        new_master_flat_name,
                        saturation_threshold,
                        ignore_bias=False):
    """Creates master flats

    Using a list of compatible flat images it combines them using median and
    1-sigma clipping. Also it apply all previous standard calibrations to
    each image.

    Args:
        flat_files (list): List of files previously filtered, there is no
        compatibility check in this function and is assumed the files are
        combinables.
        raw_data (str): Full path to raw data.
        reduced_data (str): Full path to reduced data. Where reduced data
        should be stored.
        technique (str): Observing technique. Imaging or Spectroscopy.
        overscan_region (str): Defines the area to be used to estimate the
        overscan region for overscan correction. Should be in the format.
        `[x1:x2.,y1:y2]`.
        trim_section (str):Defines the area to be used after trimming
        unusable selected parts (edges). In the format `[x1:x2.,y1:y2]`.
        master_bias_name (str): Master bias name, can be full path or not.
        If it is a relative path, the path will be ignored and will define
        the full path as `raw_path` + `basename`.
        new_master_flat_name (str): Name of the file to save new master
        flat. Can be absolute path or not.
        saturation_threshold (int): Saturation threshold, defines the percentage of
        pixels above saturation level allowed for flat field images.
        ignore_bias (bool): Flag to create master bias without master bias.


    Returns:
        The master flat :class:`~astropy.nddata.CCDData` instance and the
        name of under which the master flat was stored. If it can't build
        the master flat it will return None, None.

    """
    cleaned_flat_list = []
    master_flat_list = []

    if os.path.isabs(os.path.dirname(new_master_flat_name)):
        master_flat_name = new_master_flat_name
    else:
        master_flat_name = os.path.join(
            reduced_data, os.path.basename(new_master_flat_name))

    if not ignore_bias:
        if os.path.isabs(os.path.dirname(master_bias_name)) and \
                os.path.exists(master_bias_name):
            master_bias = read_fits(master_bias_name, technique=technique)
        else:
            master_bias_name = os.path.join(reduced_data,
                                            os.path.basename(master_bias_name))
            master_bias = read_fits(master_bias_name, technique=technique)

        master_bias = image_trim(ccd=master_bias,
                                 trim_section=trim_section,
                                 trim_type='trimsec')

    log.info('Creating Master Flat')
    for flat_file in flat_files:
        if os.path.isabs(flat_file):
            image_full_path = flat_file
        else:
            image_full_path = os.path.join(raw_data, flat_file)
        log.debug('Loading flat image: ' + image_full_path)
        ccd = read_fits(image_full_path, technique=technique)

        if ignore_bias and technique == 'Spectroscopy':
            ccd = image_overscan(ccd, overscan_region=overscan_region)

        ccd = image_trim(ccd=ccd,
                         trim_section=trim_section,
                         trim_type='trimsec')

        if not ignore_bias:
            ccd = ccdproc.subtract_bias(ccd,
                                        master_bias,
                                        add_keyword=False)
            ccd.header['GSP_BIAS'] = (
                os.path.basename(master_bias_name),
                'Master bias image')

        else:
            log.warning('Ignoring bias on request')
        if is_file_saturated(ccd=ccd,
                             threshold=saturation_threshold):
            log.warning('Removing saturated image {:s}. '
                        'Use --saturation_threshold to change saturation_threshold '
                        'level'.format(flat_file))
            continue
        else:
            cleaned_flat_list.append(flat_file)
            master_flat_list.append(ccd)

    if master_flat_list != []:
        log.info("Combining {} images to create master flat".format(
            len(master_flat_list)))
        master_flat = ccdproc.combine(master_flat_list,
                                      method='median',
                                      sigma_clip=True,
                                      sigma_clip_low_thresh=1.0,
                                      sigma_clip_high_thresh=1.0,
                                      add_keyword=False)

        # add name of images used to create master bias
        for n in range(len(cleaned_flat_list)):
            master_flat.header['GSP_IC{:02d}'.format(n + 1)] = (
                cleaned_flat_list[n],
                'Image used to create master flat')

        write_fits(ccd=master_flat,
                   full_path=master_flat_name,
                   combined=True)

        log.info('Created Master Flat: ' + master_flat_name)
        return master_flat, master_flat_name
    else:
        log.error('Empty flat list. Check that they do not exceed the '
                  'saturation_threshold limit.')
        return None, None


def cross_correlation(reference,
                      compared,
                      slit_size,
                      serial_binning,
                      mode='full',
                      plot=False):
    """Do cross correlation of two 1D spectra

    It convolves the reference lamp depending on the slit size of the new_array
    that corresponds with a comparison lamp.
    If the slit is larger than 3 arcseconds the reference lamp is convolved with
    a `~astropy.convolution.Box1DKernel` because spectral lines look more like a
    block than a line. And if it is smaller or equal to 3 it will
    use a `~astropy.convolution.Gaussian1DKernel` ponderated by the binning.
    All reference lamp are unbinned, or binning is 1x1.

    Args:
        reference (array): Reference array.
        compared (array): Array to be matched. A new reference lamp.
        slit_size (float): Slit width in arcseconds
        serial_binning (int): Binning in the spectral axis
        mode (str): Correlation mode for `scipy.signal.correlate`.
        plot (bool): Switch debugging plots on or off.

    Returns:
        correlation_value (int): Shift value in pixels.

    """
    cyaxis2 = compared
    if slit_size > 3:

        box_width = slit_size / (0.15 * serial_binning)

        log.debug('BOX WIDTH: {:f}'.format(box_width))
        box_kernel = Box1DKernel(width=box_width)
        max_before = np.max(reference)
        cyaxis1 = convolve(reference, box_kernel)
        max_after = np.max(cyaxis1)
        cyaxis1 *= max_before / max_after

    else:
        kernel_stddev = slit_size / (0.15 * serial_binning)

        gaussian_kernel = Gaussian1DKernel(stddev=kernel_stddev)
        cyaxis1 = convolve(reference, gaussian_kernel)
        cyaxis2 = convolve(compared, gaussian_kernel)

    ccorr = signal.correlate(cyaxis1, cyaxis2, mode=mode)

    max_index = np.argmax(ccorr)

    x_ccorr = np.linspace(-int(len(ccorr) / 2.),
                          int(len(ccorr) / 2.),
                          len(ccorr))

    correlation_value = x_ccorr[max_index]
    if plot:  # pragma: no cover
        plt.ion()
        plt.title('Cross Correlation')
        plt.xlabel('Lag Value')
        plt.ylabel('Correlation Value')
        plt.plot(x_ccorr, ccorr)
        plt.draw()
        plt.pause(2)
        plt.clf()
        plt.ioff()
    return correlation_value


def classify_spectroscopic_data(path, search_pattern):
    """Classify data by grouping them by a set of keywords.

    This function uses :class:`~ccdproc.ImageFileCollection`. First it creates a
    collection of information regarding the images located in ``path`` that
    match the pattern ``search_pattern``. The information obtained are all
    keywords listed in the list ``keywords``.
    The :class:`~ccdproc.ImageFileCollection` object is translated into
    :class:`~pandas.DataFrame` and then is used much like an SQL database to
    select and filter values and in that way put them in groups that are
    :class:`~pandas.DataFrame` instances.
    The keywords retrieved are:
    - ``date``
    - ``slit``
    - ``date-obs``
    - ``obstype``
    - ``object``
    - ``exptime``
    - ``obsra``
    - ``obsdec``
    - ``grating``
    - ``cam_targ``
    - ``grt_targ``
    - ``filter``
    - ``filter2``
    - ``gain``
    - ``rdnoise``.

    Then all data is grouped by matching the following keywords:

    - ``slit``
    - ``radeg``
    - ``decdeg``
    - ``grating``
    - ``cam_targ``
    - ``grt_targ``
    - ``filter``
    - ``filter2``
    - ``gain``
    - ``rdnoise``

    And finally,  every group is classified as: a *comparison lamp-only* group,
    an *object-only* group or a *group of object and comparison lamps*. The
    comparison lamps present in the last group (``COMP`` + ``OBJECT``) are also
    added in the first one (``COMP``-only).


    Args:
        path (str): Path to data location
        search_pattern (str): Prefix to match files.

    Returns:
        Instance of :class:`goodman_pipeline.core.core.NightDataContainer`

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
                'rdnoise',
                'lamp_hga',
                'lamp_ne',
                'lamp_ar',
                'lamp_fe',
                'lamp_cu'
                ]

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

        if any([value in ['COMP', 'ARC'] for value in group_obstype]) and \
                len(group_obstype) == 1:
            log.debug('Adding COMP group')
            data_container.add_comp_group(comp_group=spec_group)
        elif any([value in ['OBJECT', 'SPECTRUM'] for value in group_obstype]) \
                and len(group_obstype) == 1:
            log.debug('Adding OBJECT group')
            data_container.add_object_group(object_group=spec_group)
        else:
            log.debug('Adding OBJECT-COMP group')
            data_container.add_spec_group(spec_group=spec_group)

    return data_container


def combine_data(image_list, dest_path, prefix=None, output_name=None,
                 method="median",
                 save=False):
    """Combine a list of :class:`~astropy.nddata.CCDData` instances.

    Args:
        image_list (list): Each element should be an instance of
          :class:`~astropy.nddata.CCDData`
        dest_path (str): Path to where the new image should saved
        prefix (str): Prefix to add to the image file name
        output_name (str): Alternatively a file name can be parsed, this will
          ignore `prefix`.
        method (str): Method for doing the combination, this goes straight to
          the call of `ccdproc.combine` function.
        save (bool): If True will save the combined images. If False it will
          ignore `prefix` or `output_name`.

    Returns:
        A combined image as a :class:`~astropy.nddata.CCDData` object.

    """
    assert len(image_list) > 1

    combined_full_path = os.path.join(dest_path, 'combined_file.fits')

    if output_name is not None:
        combined_full_path = os.path.join(dest_path, output_name)
    elif prefix is not None:
        sample_image_name = random.choice(image_list).header['GSP_FNAM']
        splitted_name = sample_image_name.split('_')
        splitted_name[0] = re.sub('_', '', prefix)
        splitted_name[1] = 'comb'
        splitted_name[-1] = re.sub('.fits', '', splitted_name[-1])

        combined_base_name = "_".join(splitted_name)

        number = len(glob.glob(
            os.path.join(dest_path,
                         combined_base_name + "*.fits")))

        combined_full_path = os.path.join(
            dest_path,
            combined_base_name + "_{:03d}.fits".format(number + 1))

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
        new_image_name = '_' + image_name
        if os.path.isfile(os.path.join(dest_path, image_name)):
            write_fits(image_list[i],
                       full_path=os.path.join(dest_path,
                                              new_image_name))
            log.info("Deleting file {}".format(image_name))
            os.unlink(os.path.join(dest_path, image_name))
        else:
            log.error("File {} does not exists".format(
                os.path.join(dest_path,
                             image_name)))
        combined_image.header.set("GSP_IC{:02d}".format(i + 1),
                                  value=new_image_name,
                                  comment='Image used to create combined')

    if save:
        write_fits(combined_image,
                   full_path=combined_full_path,
                   combined=True)
        log.info("Saved combined file to {}".format(combined_full_path))

    return combined_image


def convert_time(in_time):
    """Converts time to seconds since epoch

    Args:
        in_time (str): time obtained from header's keyword DATE-OBS

    Returns:
        time in seconds since epoch

    """
    return calendar.timegm(time.strptime(in_time, "%Y-%m-%dT%H:%M:%S.%f"))


def dcr_cosmicray_rejection(data_path, in_file, prefix,
                            keep_cosmic_files=False, save=True):
    """Runs an external code for cosmic ray rejection

    DCR was created by Wojtek Pych and the code can be obtained from
    http://users.camk.edu.pl/pych/DCR/ and is written in C. Contrary to
    ccdproc's LACosmic it actually applies the correction, and also doesn't
    update the mask attribute since it doesn't work with :class:`~astropy.nddata.CCDData` instances.

    The binary takes three positional arguments, they are: 1. input image,
    2. output image and 3. cosmic rays image. Also it needs that a dcr.par file
    is located in the directory. All this is implemented in this function, if
    `delete` is True it will remove the original image and the cosmic rays
    image. The removal of the original image is absolutely safe when used in the
    context of the goodman pipeline, however if you want to implement it
    somewhere else, be careful.

    Notes:
        This function operates an external code therefore it doesn't return
        anything natively, instead it creates a new image. A workaround has been
        created that loads the new image and deletes the file.

    Args:
        data_path (str): Data location
        in_file (str): Name of the file to have its cosmic rays removed
        prefix (str): Prefix to add to the file with the cosmic rays removed
        keep_cosmic_files (bool): True for deleting the input and cosmic ray
          file.
        save (bool): Toggles the option of saving the image.

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

    log.debug('DCR command:')
    log.debug(command)
    # get the current working directory to go back to it later in case the
    # the pipeline has not been called from the same data directory.
    cwd = os.getcwd()

    # move to the directory were the data is, dcr is expecting a file dcr.par
    os.chdir(data_path)

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
    def kill_process(process):  # pragma: no cover
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
    if not keep_cosmic_files and stderr == b'' and b'USAGE:' not in stdout \
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

    # recovers the saved file and returns the :class:`~astropy.nddata.CCDData`
    # instance
    if os.path.isfile(full_path_out):
        ccd = CCDData.read(full_path_out, unit=u.adu)
        if not save:
            log.warning("Removing file because the attribute 'save' "
                        "is set to False")
            os.unlink(full_path_out)
        return ccd


def define_trim_section(sample_image, technique):
    """Get the initial trim section

    The initial trim section is usually defined in the header with the
    keyword ``TRIMSEC`` but in the case of Goodman HTS this does not work well.
    In particular for spectroscopy where is more likely to have combined
    binning and so on.

    Args:
        sample_image (str): Full path to sample image.
        technique (str): The name of the technique, the options are:
            Imaging or Spectroscopy.

    Returns:
        The trim section in the format ``[x1:x2, y1:y2]``

    """

    assert os.path.isabs(os.path.dirname(sample_image))
    assert os.path.isfile(sample_image)

    trim_section = None
    # TODO (simon): Consider binning and possibly ROIs for trim section
    log.warning('Determining trim section. Assuming you have only one '
                'kind of data in this folder')

    ccd = read_fits(sample_image, technique=technique)

    # serial binning - dispersion binning
    # parallel binning - spatial binning
    spatial_length, dispersion_length = ccd.data.shape
    serial_binning, \
    parallel_binning = [int(x) for x
                        in ccd.header['CCDSUM'].split()]

    # Trim section is valid for Blue and Red Camera Binning 1x1 and
    # Spectroscopic ROI
    if technique == 'Spectroscopy':

        # left
        low_lim_spectral = int(np.ceil(51. / serial_binning))

        # right
        high_lim_spectral = int(4110 / serial_binning)

        # bottom
        low_lim_spatial = 2

        # top
        # t = int(1896 / parallel_binning)
        # TODO (simon): Need testing
        # trim_section = '[{:d}:{:d},{:d}:{:d}]'.format(l, r, b, t)
        trim_section = '[{:d}:{:d},{:d}:{:d}]'.format(
            low_lim_spectral,
            high_lim_spectral,
            low_lim_spatial,
            spatial_length)

    elif technique == 'Imaging':
        trim_section = ccd.header['TRIMSEC']

    log.info('Trim Section: %s', trim_section)
    return trim_section


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
        ccd (CCDData): Instance of :class:`~astropy.nddata.CCDData` containing a
          2D spectrum
        target_trace (object): Instance of astropy.modeling.Model, a low order
          polynomial that defines the trace of the spectrum in the ccd object.
        spatial_profile (Model): Instance of :class:`~astropy.modeling.Model`,
          a Gaussian model previously fitted to the spatial profile of the 2D
          spectrum contained in the ccd object.
        extraction_name (str): Extraction type, can be `fractional` or
          `optimal` though the optimal extraction is not implemented yet.

    Returns:
        ccd (CCDData): Instance of :class:`~astropy.nddata.CCDData` containing a
        1D spectrum. The attribute 'data' is replaced by the 1D array resulted
        from the extraction process.

    Raises:
        NotImplementedError: When `extraction_name` is `optimal`.

    """
    assert isinstance(ccd, CCDData)
    assert isinstance(target_trace, Model)

    if spatial_profile.__class__.name == 'Gaussian1D':
        target_fwhm = spatial_profile.fwhm
    elif spatial_profile.__class__.name == 'Moffat1D':
        target_fwhm = spatial_profile.fwhm
    else:
        raise NotImplementedError

    if extraction_name == 'fractional':
        extracted, background, bkg_info = extract_fractional_pixel(
            ccd=ccd,
            target_trace=target_trace,
            target_fwhm=target_fwhm,
            extraction_width=2)

        background_1, background_2 = bkg_info

        if background_1 is not None:
            log.info('Background extraction zone 1: {:s}'.format(background_1))
            extracted.header.set('GSP_BKG1', value=background_1)
        else:
            log.info("Background extraction zone 1: 'none'")

        if background_2 is not None:
            log.info('Background extraction zone 2: {:s}'.format(background_2))
            extracted.header.set('GSP_BKG2', value=background_2)
        else:
            log.info("Background extraction zone 2: 'none'")

        return extracted

    elif extraction_name == 'optimal':
        raise NotImplementedError


def extract_fractional_pixel(ccd, target_trace, target_fwhm, extraction_width,
                             background_spacing=3):
    """Performs an spectrum extraction using fractional pixels.

    Args:
        ccd (CCDData) Instance of :class:`~astropy.nddata.CCDData` that
        contains a 2D spectrum.
        target_trace (object):  Instance of astropy.modeling.models.Model that
          defines the trace of the target on the image (ccd).
        target_fwhm (float): FWHM value for the spatial profile
          fitted to the target.
        extraction_width (int): Width of the extraction area as a function of
          `target_fwhm`. For instance if `extraction_with` is set to 1 the
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

    apnum1 = None

    background_info_1 = None
    background_info_2 = None

    non_background_sub = []
    extracted_spectrum = []
    background_list = []

    if ccd.header['OBSTYPE'] not in ['OBJECT', 'SPECTRUM']:
        log.debug("No background subtraction for OBSTYPE = "
                  "{:s}".format(ccd.header['OBSTYPE']))

    for i in disp_axis:

        # this defines the extraction limit for every column
        low_limit = trace_points[i] - 0.5 * extraction_width * target_fwhm
        high_limit = trace_points[i] + 0.5 * extraction_width * target_fwhm

        if apnum1 is None:
            # TODO (simon): add secondary targets

            apnum1 = '{:d} {:d} {:.2f} {:.2f}'.format(1,
                                                      1,
                                                      low_limit,
                                                      high_limit)

            ccd.header.set('APNUM1',
                           value=apnum1,
                           comment="Aperture in first column")

            ccd.header.set('GSP_EXTR',
                           value="{:.2f}:{:.2f}".format(low_limit,
                                                        high_limit))

            log.info("Extraction aperture in first column: {:s}".format(
                ccd.header['GSP_EXTR']))

        column_sum = fractional_sum(data=ccd.data,
                                    index=i,
                                    low_limit=low_limit,
                                    high_limit=high_limit)

        non_background_sub.append(column_sum)

        if ccd.header['OBSTYPE'] in ['OBJECT', 'SPECTRUM']:
            # background limits

            # background_spacing is the distance from the edge of the target's
            # limits defined by `int_low_limit` and
            # `int_high_limit` in stddev units

            background_width = high_limit - low_limit

            # define pixel values for background subtraction
            # low_background_zone
            high_1 = low_limit - background_spacing * target_fwhm
            low_1 = high_1 - background_width

            # High background zone
            low_2 = high_limit + background_spacing * target_fwhm
            high_2 = low_2 + background_width

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
                log.error("Invalid Zone 1: [{}:{}]".format(low_1, high_1))

            if high_2 < spat_length and neighbouring_target_condition:
                background_2 = fractional_sum(data=ccd.data,
                                              index=i,
                                              low_limit=low_2,
                                              high_limit=high_2)
            else:
                log.error("Invalid Zone 2: [{}:{}]".format(low_2, high_2))

            # background = 0
            if background_1 is not None and background_2 is None:
                background = background_1
                if background_info_1 is None:
                    background_info_1 = "{:.2f}:{:.2f} column {:d}".format(
                        low_1, high_1, i+1)
            elif background_1 is None and background_2 is not None:
                background = background_2
                if background_info_2 is None:
                    background_info_2 = "{:.2f}:{:.2f} column {:d}".format(
                        low_2, high_2, i+1)
            else:
                background = np.mean([background_1, background_2])
                if background_info_1 is None:
                    background_info_1 = "{:.2f}:{:.2f} column {:d}".format(
                        low_1, high_1, i+1)
                if background_info_2 is None:
                    background_info_2 = "{:.2f}:{:.2f} column {:d}".format(
                        low_2, high_2, i+1)

            # actual background subtraction
            background_subtracted_column_sum = column_sum - background

            # append column value to list
            extracted_spectrum.append(background_subtracted_column_sum)
            background_list.append(background)
        else:
            extracted_spectrum.append(column_sum)

    new_ccd = ccd.copy()
    new_ccd.data = np.asarray(extracted_spectrum)
    if new_ccd.header['NAXIS'] != 1:
        for i in range(int(new_ccd.header['NAXIS']), 1, -1):
            new_ccd.header.remove(keyword="NAXIS{:d}".format(i))
        new_ccd.header.set('NAXIS', value=1)
    return new_ccd, np.asarray(background_list), [background_info_1,
                                                  background_info_2]


def extract_optimal():
    """Placeholder for optimal extraction method.

    Raises:
        NotImplementedError
    """
    raise NotImplementedError


def evaluate_wavelength_solution(clipped_differences):
    """Calculates Root Mean Square Error for the wavelength solution.

    Args:
        clipped_differences (ndarray): Numpy masked array of differences
          between reference line values in angstrom and the value calculated
          using the model of the wavelength solution.

    Returns:
        Root Mean Square Error, number of points and number of points
          rejected in the calculation of the wavelength solution.

    """
    n_points = len(clipped_differences)
    n_rejections = np.ma.count_masked(clipped_differences)
    square_differences = []

    for i in range(len(clipped_differences)):
        if clipped_differences[i] is not np.ma.masked:
            square_differences.append(clipped_differences[i] ** 2)

    rms_error = np.sqrt(
        np.sum(square_differences) / len(square_differences))

    log.info('Wavelength solution RMS Error : {:.3f}'.format(
        rms_error))

    return rms_error, n_points, n_rejections


def fix_keywords(path, pattern='*.fits'):
    """Fix FITS uncompliance of some keywords

    Uses automatic header fixing by :class:`~astropy.nddata.CCDData`. Note that
    this only fixes FITS compliance.

    Args:
        path (str): Path to raw data
        pattern (str): Search pattern for listing file in path.

    """
    file_list = glob.glob(os.path.join(path, pattern))
    for _file in file_list:
        log.info("Fixing file {:s}".format(_file))
        ccd = CCDData.read(_file, unit='adu')

        ccd.write(_file, overwrite=True)
        log.info("Fix succeeded!")


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

    return column_sum


def get_best_flat(flat_name, path):
    """Look for matching master flat

    Given a basename for master flats defined as a combination of key parameters
    extracted from the header of the image that we want to flat field, this
    function will find the name of the files that matches the base name and then
    will choose the first. Ideally this should go further as to check signal,
    time gap, etc.
    After it identifies the file it will load it using
    :class:`~astropy.nddata.CCDData` and return it along the filename.
    In the case it fails it will return None instead of master_flat and another
    None instead of master_flat_name.

    Args:
        flat_name (str): Full path of master flat basename. Ends in '*.fits'
          for using glob.
        path (str): Location to look for flats.

    Returns:
        master_flat (object): A :class:`~astropy.nddata.CCDData` instance.
        master_flat_name (str): Full path to the chosen master flat.

    """
    flat_list = glob.glob(os.path.join(path, flat_name))
    log.debug('Flat base name {:s}'.format(flat_name))
    log.debug('Matching master flats found: {:d}'.format(len(flat_list)))
    if len(flat_list) > 0:
        master_flat_name = flat_list[0]
        # if len(flat_list) == 1:
        #     master_flat_name = flat_list[0]
        # else:
        #     master_flat_name = flat_list[0]
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

    grating_frequency = float(grating) / u.mm
    grt_ang = float(grt_ang) * u.deg
    cam_ang = float(cam_ang) * u.deg

    alpha = grt_ang.to(u.rad)
    beta = cam_ang.to(u.rad) - grt_ang.to(u.rad)

    # central_wavelength = (1e6 / grating_frequency) * \
    #                      (np.sin(alpha * np.pi / 180.) +
    #                       np.sin(beta * np.pi / 180.))
    central_wavelength = (np.sin(alpha) + np.sin(beta)) / grating_frequency
    central_wavelength = central_wavelength.to(u.angstrom)
    log.debug('Found {:.3f} as central wavelength'.format(central_wavelength))

    return central_wavelength


def get_lines_in_lamp(ccd, plots=False):
    """Identify peaks in a lamp spectrum

    Uses `scipy.signal.argrelmax` to find peaks in a spectrum i.e emission
    lines, then it calls the recenter_lines method that will recenter them
    using a "center of mass", because, not always the maximum value (peak)
    is the center of the line.

    Args:
        ccd (CCDData): Lamp `ccdproc.CCDData` instance.
        plots (bool): Wether to plot or not.

    Returns:
        lines_candidates (list): A common list containing pixel values at
            approximate location of lines.

    """
    if isinstance(ccd, CCDData):
        lamp_data = ccd.data
        lamp_header = ccd.header
        raw_pixel_axis = range(len(lamp_data))
    else:
        log.error('Error receiving lamp')
        return None

    no_nan_lamp_data = np.asarray(np.nan_to_num(lamp_data))

    filtered_data = np.where(
        np.abs(no_nan_lamp_data > no_nan_lamp_data.min() +
               0.03 * no_nan_lamp_data.max()),
        no_nan_lamp_data,
        None)

    # replace None to zero and convert it to an array
    none_to_zero = [0 if it is None else it for it in filtered_data]
    filtered_data = np.array(none_to_zero)

    _upper_limit = no_nan_lamp_data.min() + 0.03 * no_nan_lamp_data.max()
    slit_size = np.float(re.sub('[a-zA-Z"_*]', '', lamp_header['slit']))

    serial_binning, parallel_binning = [
        int(x) for x in lamp_header['CCDSUM'].split()]

    new_order = int(round(float(slit_size) / (0.15 * serial_binning)))
    log.debug('New Order:  {:d}'.format(new_order))

    peaks = signal.argrelmax(filtered_data, axis=0, order=new_order)[0]

    if slit_size >= 5.:

        lines_center = recenter_broad_lines(
            lamp_data=no_nan_lamp_data,
            lines=peaks,
            order=new_order)
    else:
        # lines_center = peaks
        lines_center = recenter_lines(no_nan_lamp_data, peaks)

    if plots:  # pragma: no cover
        plt.close('all')
        fig, ax = plt.subplots()
        fig.canvas.set_window_title('Lines Detected')

        mng = plt.get_current_fig_manager()
        mng.window.showMaximized()

        ax.set_title('Lines detected in Lamp\n'
                     '{:s}'.format(lamp_header['OBJECT']))
        ax.set_xlabel('Pixel Axis')
        ax.set_ylabel('Intensity (counts)')

        # Build legends without data to avoid repetitions
        ax.plot([], color='k', label='Comparison Lamp Data')

        ax.plot([], color='k', linestyle=':',
                label='Spectral Line Detected')

        ax.axhline(_upper_limit, color='r')

        for line in peaks:
            ax.axvline(line, color='k', linestyle=':')

        ax.plot(raw_pixel_axis, no_nan_lamp_data, color='k')
        ax.legend(loc='best')
        plt.tight_layout()
        plt.show()

    return lines_center


def get_overscan_region(sample_image, technique):
    """Get the right overscan region for spectroscopy

    It works for the following ROI:
        Spectroscopic 1x1
        Spectroscopic 2x2
        Spectroscopic 3x3

    The limits where measured on a Spectroscopic 1x1 image and then divided
    by the binning size. This was checked
    that it actually works as expected.

    Notes:
        The regions are 1-based i.e. different to Python convention.
        For Imaging there is no overscan region.

    Args:
          sample_image (str): Full path to randomly chosen image.
          technique (str): Observing technique, either `Spectroscopy` or
          `Imaging`


    Returns:
        overscan_region (str) Region for overscan in the format
          '[min:max,:]' where min is the starting point and max is the end
           point of the overscan region.

    """

    assert os.path.isabs(os.path.dirname(sample_image))
    assert os.path.isfile(sample_image)

    log.debug('Overscan Sample File ' + sample_image)
    ccd = CCDData.read(sample_image, unit=u.adu)

    # Image height - spatial direction
    spatial_length, dispersion_length = ccd.data.shape

    # Image width - spectral direction
    # w = ccd.data.shape[1]

    # Take the binnings
    serial_binning, parallel_binning = \
        [int(x) for x in ccd.header['CCDSUM'].split()]

    if technique == 'Spectroscopy':
        log.info('Overscan regions has been tested for ROI '
                 'Spectroscopic 1x1, 2x2 and 3x3')

        # define l r b and t to avoid local variable might be
        # defined before assignment warning
        low_lim_spectral, \
        high_lim_spectral, \
        low_lim_spatial, \
        high_lim_spatial = [None] * 4
        if ccd.header['INSTCONF'] == 'Red':
            # for red camera it is necessary to eliminate the first
            # rows/columns (depends on the point of view) because
            # they come with an abnormal high signal. Usually the
            # first 5 pixels. In order to find the corresponding
            # value for the subsequent binning divide by the
            # binning size.
            # The numbers 6 and 49 where obtained from visual
            # inspection

            # left
            low_lim_spectral = int(np.ceil(6. / serial_binning))
            # right
            high_lim_spectral = int(49. / serial_binning)
            # bottom
            low_lim_spatial = 1
            # top
            high_lim_spatial = spatial_length
        elif ccd.header['INSTCONF'] == 'Blue':
            # 16 is the length of the overscan region with no
            # binning.

            # left
            low_lim_spectral = 1
            # right
            high_lim_spectral = int(16. / serial_binning)
            # bottom
            low_lim_spatial = 1
            # top
            high_lim_spatial = spatial_length

        overscan_region = '[{:d}:{:d},{:d}:{:d}]'.format(
            low_lim_spectral,
            high_lim_spectral,
            low_lim_spatial,
            high_lim_spatial)

    elif technique == 'Imaging':
        log.warning("Imaging mode doesn't have overscan "
                    "region. Use bias instead.")
        overscan_region = None
    else:
        overscan_region = None
    log.info('Overscan Region: %s', overscan_region)
    return overscan_region


def get_slit_trim_section(master_flat):
    """Find the slit edges to trim all data

    Using a master flat, ideally with good signal to noise ratio, this function
    will identify the edges of the slit projected into the detector. Having this
    done will allow to reduce the overall processing time and also reduce the
    introduction of artifacts due to non-illuminated regions in the detectors,
    such as NaNs -INF +INF, etc.

    Args:
        master_flat (CCDData): A :class:`~astropy.nddata.CCDData` instance.

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

    low_lim = int(np.max([1 + offset, l_lim + 1]))

    high_lim = int(np.min([h_lim, len(ccd_section_median) - offset]))

    # define the slit trim section as (IRAF)
    # convert o 1-based
    slit_trim_section = '[1:{:d},{:d}:{:d}]'.format(y,
                                                    low_lim,
                                                    high_lim)
    log.debug("Slit Trim Section: {:s}".format(slit_trim_section))

    # debugging plots that have to be manually turned on
    if False:  # pragma: no cover
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

    return slit_trim_section


def get_spectral_characteristics(ccd, pixel_size, instrument_focal_length):
    """Calculates some Goodman's specific spectroscopic values.

    From the header value for Grating, Grating Angle and Camera Angle it is
    possible to estimate what are the wavelength values at the edges as well
    as in the center. It was necessary to add offsets though, since the
    formulas provided are slightly off. The values are only an estimate.

    Args:
        ccd (CCDData): Lamp `ccdproc.CCDData` instance
        pixel_size (float): Pixel size in microns
        instrument_focal_length (float): Instrument focal length

    Returns:
        spectral_characteristics (dict): Contains the following parameters:
            center: Center Wavelength
            blue: Blue limit in Angstrom
            red: Red limit in Angstrom
            alpha: Angle
            beta: Angle
            pix1: Pixel One
            pix2: Pixel Two


    """
    # TODO (simon): find a definite solution for this, this only work
    # TODO (simon): (a little) for one configuration
    blue_correction_factor = -50 * u.angstrom
    red_correction_factor = -37 * u.angstrom

    grating_frequency = float(re.sub('[A-Za-z_-]',
                                     '',
                                     ccd.header['GRATING'])) / u.mm

    grating_angle = float(ccd.header['GRT_ANG']) * u.deg
    camera_angle = float(ccd.header['CAM_ANG']) * u.deg

    # serial binning - dispersion binning
    # parallel binning - spatial binning
    serial_binning, parallel_binning = [
        int(x) for x in ccd.header['CCDSUM'].split()]

    pixel_count = len(ccd.data)
    # Calculations
    # TODO (simon): Check whether is necessary to remove the
    # TODO (simon): slit_offset variable
    alpha = grating_angle.to(u.rad)
    beta = camera_angle.to(u.rad) - grating_angle.to(u.rad)

    center_wavelength = (np.sin(alpha) + np.sin(beta)) / grating_frequency

    limit_angle = np.arctan(pixel_count * ((pixel_size * serial_binning) / instrument_focal_length) / 2)

    blue_limit = ((np.sin(alpha) + np.sin(beta - limit_angle.to(u.rad))) / grating_frequency).to(u.angstrom) + blue_correction_factor

    red_limit = ((np.sin(alpha) + np.sin(beta + limit_angle.to(u.rad))) / grating_frequency).to(u.angstrom) + red_correction_factor

    pixel_one = 0
    pixel_two = 0
    log.debug(
        'Center Wavelength : {:.3f} Blue Limit : '
        '{:.3f} Red Limit : {:.3f} '.format(center_wavelength.to(u.angstrom),
                                            blue_limit,
                                            red_limit))

    spectral_characteristics = {'center': center_wavelength,
                                'blue': blue_limit,
                                'red': red_limit,
                                'alpha': alpha,
                                'beta': beta,
                                'pix1': pixel_one,
                                'pix2': pixel_two}
    return spectral_characteristics


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


def identify_targets(ccd,
                     fit_model,
                     background_threshold,
                     nfind=3,
                     plots=False):
    """Identify Spectroscopic Targets

    Wrapper to the class `IdentifySpectroscopicTargets`.

    Args:
        ccd (CCDData): Image containing spectra
        fit_model (str): Name of the model to be fitted `moffat` or `gaussian`.
        background_threshold (int): Number of background levels for target
        discrimination.
        nfind (int): Maximum number of targets passing the background threshold
        to be returned, they are order from most intense peak to least intense.
        plots (bool): Flat for plotting results.

    Returns:
        identified_targets (list): List of models successfully fitted.

    """
    identify = IdentifySpectroscopicTargets()

    identified_targets = identify(ccd=ccd,
                                  nfind=nfind,
                                  background_threshold=background_threshold,
                                  model_name=fit_model,
                                  plots=plots)

    return identified_targets


def identify_technique(target, obstype, slit, grating, wavmode, roi):
    """Identify whether is Imaging or Spectroscopic data

    Args:
        target (str): Target name as in the keyword `OBJECT` this is useful in
        Automated aquisition mode, such as AEON.
        obstype (str): Observation type as in `OBSTYPE`
        slit (str): Value of `SLIT` keyword.
        grating (str): Value of `GRATING` keyword.
        wavmode (str): Value of `WAVMODE` keyword.
        roi (str): Value of `ROI` keyword.

    Returns:
        Observing technique as a string. Either `Imaging` or `Spectroscopy`.

    """
    if 'Spectroscopic' in roi or \
            obstype in ['ARC', 'SPECTRUM', 'COMP'] or \
            slit not in ['NO_MASK', '<NO MASK>'] or \
            grating not in ['NO_GRATING', '<NO GRATING>'] or \
            '_SP_' in target:
        technique = 'Spectroscopy'
    elif 'Imaging' in roi or \
            obstype in ['EXPOSE'] or\
            wavmode == 'IMAGING' or '_IM_' in target:
        technique = 'Imaging'
    else:
        technique = 'Unknown'
    return technique


def image_overscan(ccd, overscan_region, add_keyword=False):
    """Apply overscan correction to data

    Uses ccdproc.subtract_overscan to perform the task.

    Notes:
        The overscan_region argument uses FITS convention, just like IRAF,
        therefore is 1 based. i.e. it starts in 1 not 0.

    Args:
        ccd (CCDData) A :class:`~astropy.nddata.CCDData` instance to be
          overscan corrected.
        overscan_region (str): The overscan region in the format `[x1:x2,y1:y2]`
          where x is the spectral axis and y is the spatial axis.
        add_keyword (bool): Tells ccdproc whether to add a keyword or not.
          Default False.

    Returns:
        ccd (CCDData) Overscan corrected :class:`~astropy.nddata.CCDData`
          instance

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


def image_trim(ccd, trim_section, trim_type='trimsec', add_keyword=False):
    """Trim image to a given section

    Notes:
        The overscan_region argument uses FITS convention, just like IRAF,
        therefore is 1 based. i.e. it starts in 1 not 0.

    Args:
        ccd (CCDData) A :class:`~astropy.nddata.CCDData` instance.
        trim_section (str): The trimming section in the format `[x1:x2,y1:y2]`
          where x is the spectral axis and y is the spatial axis.
        trim_type (str): trimsec or slit trim.
        add_keyword (bool): Tells ccdproc whether to add a keyword or not.
          Default False.

    Returns:
        ccd (CCDData) Trimmed :class:`~astropy.nddata.CCDData` instance

    """
    if trim_section is not None:
        ccd = ccdproc.trim_image(ccd=ccd,
                                 fits_section=trim_section,
                                 add_keyword=add_keyword)
        if trim_type == 'trimsec':
            ccd.header['GSP_TRIM'] = (trim_section, 'Trim section from TRIMSEC')
        elif trim_type == 'slit':
            ccd.header['GSP_SLIT'] = (trim_section,
                                      'Slit trim section, slit illuminated '
                                      'area only.')
        else:
            log.warning('Unrecognized trim type: {}'.format(trim_type))
            ccd.header['GSP_TRIM'] = (trim_section,
                                      'Image trimmed by unrecognized method: '
                                      '{:s}'.format(trim_type))
    else:
        log.info("{:s} trim section is not "
                 "defined.".format(trim_type.capitalize()))
        log.debug("Trim section is None, returning the same data.")

    return ccd


def interpolate(spectrum, interpolation_size):
    """Creates an interpolated version of the input spectrum

    This method creates an interpolated version of the input array, it is
    used mainly for a spectrum but it can also be used with any
    unidimensional array. The reason for doing interpolation is
    that it allows to find the lines and its respective center more
    precisely.

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


def is_file_saturated(ccd, threshold):
    """Detects a saturated image

    It counts the number of pixels above the saturation_threshold level, then finds
    which percentage they represents and if it is above the threshold it
    will return True. The percentage threshold can be set using the command
    line argument ``--saturation_threshold``.

    Args:
        ccd (CCDData): Image to be tested for saturation_threshold
        threshold (float): Percentage of saturated pixels allowed. Default 1.

    Returns:
        True for saturated and False for non-saturated

    """

    saturation_values = SaturationValues()

    pixels_above_saturation = np.count_nonzero(
        ccd.data[np.where(
            ccd.data > saturation_values.get_saturation_value(
                ccd=ccd))])

    total_pixels = np.count_nonzero(ccd.data)

    saturated_percent = (pixels_above_saturation * 100.) / total_pixels

    if saturated_percent >= float(threshold):
        log.warning(
            "The current image has more than {:.2f} percent "
            "of pixels above saturation_threshold level".format(float(threshold)))
        return True
    else:
        return False


def linearize_spectrum(data, wavelength_solution, plots=False):
    """Produces a linearized version of the spectrum

    Storing wavelength solutions in a FITS header is not simple at all for
    non-linear solutions therefore is easier for the final user and for the
    development code to have the spectrum linearized. It first finds a
    spline representation of the data, then creates a linear wavelength axis
    (angstrom) and finally it resamples the data from the spline
    representation to the linear wavelength axis.

    It also applies a median filter of kernel size three to smooth the
    linearized spectrum. Sometimes the splines produce funny things when
    the original data is too steep.

    Args:
        data (Array): The non-linear spectrum
        wavelength_solution (object): Mathematical model representing the
        wavelength solution.
        plots (bool): Whether to show the plots or not.

    Returns:
        linear_data (list): Contains two elements: Linear wavelength axis
        and the smoothed linearized data itself.

    """
    pixel_axis = range(len(data))
    if any(np.isnan(data)):
        log.error("there are nans")
        sys.exit(0)

    if wavelength_solution is not None:
        x_axis = wavelength_solution(pixel_axis)
        try:  # pragma: no cover
            plt.imshow(data)
            plt.show()
        except TypeError:
            pass
        new_x_axis = np.linspace(x_axis[0], x_axis[-1], len(data))
        tck = scipy.interpolate.splrep(x_axis, data, s=0)
        linearized_data = scipy.interpolate.splev(new_x_axis,
                                                  tck,
                                                  der=0)

        smoothed_linearized_data = signal.medfilt(linearized_data)
        if plots:  # pragma: no cover
            fig6 = plt.figure(6)
            plt.xlabel('Wavelength (Angstrom)')
            plt.ylabel('Intensity (Counts)')
            fig6.canvas.set_window_title('Linearized Data')

            plt.plot(x_axis,
                     data,
                     color='k',
                     label='Data')

            plt.plot(new_x_axis,
                     linearized_data,
                     color='r',
                     linestyle=':',
                     label='Linearized Data')

            plt.plot(new_x_axis,
                     smoothed_linearized_data,
                     color='m',
                     alpha=0.5,
                     label='Smoothed Linearized Data')

            fig6.tight_layout()
            plt.legend(loc=3)
            plt.show()

            fig7 = plt.figure(7)
            plt.xlabel('Pixels')
            plt.ylabel('Angstroms')
            fig7.canvas.set_window_title('Wavelength Solution')
            plt.plot(x_axis, color='b', label='Non linear wavelength-axis')
            plt.plot(new_x_axis, color='r', label='Linear wavelength-axis')
            fig7.tight_layout()
            plt.legend(loc=3)
            plt.show()

        linear_data = [new_x_axis, smoothed_linearized_data]
        return linear_data


def name_master_flats(header,
                      technique,
                      reduced_data,
                      sun_set,
                      sun_rise,
                      evening_twilight,
                      morning_twilight,
                      target_name='',
                      get=False):
    """Defines the name of a master flat or what master flat is compatible
    with a given data

    Given the header of a flat image this method will look for certain
    keywords that are unique to a given instrument configuration therefore
    they are used to discriminate compatibility.

    It can be used to define a master flat's name when creating it or find
    a base name to match existing master flat files thus finding a
    compatible one for a given non-flat image.

    Args:
        header (object): Fits header. Instance of
        :class:`~astropy.io.fits.header.Header`
        technique (str): Observing technique, either Spectroscopy or
        Imaging.
        reduced_data (str): Full path to reduced data directory
        sun_set (str): Sunset time formatted as "%Y-%m-%dT%H:%M:%S.%f"
        sun_rise (str): Sunrise time formatted as "%Y-%m-%dT%H:%M:%S.%f"
        evening_twilight (str): End of evening twilight formatted as
        "%Y-%m-%dT%H:%M:%S.%f"
        morning_twilight (str): Start of morning twilight in the format
        "%Y-%m-%dT%H:%M:%S.%f"
        target_name (str): Optional science target name to be added to the
            master flat name.
        get (bool): This option is used when trying to find a suitable
            master flat for a given data.

    Returns:
        A master flat name, or basename to find a match among existing
        files.

    """
    master_flat_name = os.path.join(reduced_data, 'master_flat')
    sunset = datetime.datetime.strptime(sun_set,
                                        "%Y-%m-%dT%H:%M:%S.%f")

    sunrise = datetime.datetime.strptime(sun_rise,
                                         "%Y-%m-%dT%H:%M:%S.%f")

    afternoon_twilight = datetime.datetime.strptime(evening_twilight,
                                                    "%Y-%m-%dT%H:%M:%S.%f")

    morning_twilight = datetime.datetime.strptime(morning_twilight,
                                                  "%Y-%m-%dT%H:%M:%S.%f")

    date_obs = datetime.datetime.strptime(header['DATE-OBS'],
                                          "%Y-%m-%dT%H:%M:%S.%f")
    if target_name != '':
        target_name = '_' + target_name
    if not get:
        # TODO (simon): There must be a pythonic way to do this
        if afternoon_twilight < date_obs < morning_twilight:
            dome_sky = '_night'
        elif (sunset < date_obs < afternoon_twilight) or \
                (morning_twilight < date_obs < sunrise):
            dome_sky = '_sky'
        else:
            dome_sky = '_dome'
    else:
        dome_sky = '*'

    if technique == 'Spectroscopy':
        if header['GRATING'] != '<NO GRATING>':
            flat_grating = '_' + re.sub('[A-Za-z_ ]',
                                        '',
                                        header['GRATING'])

            # self.spec_mode is an instance of SpectroscopicMode
            spectroscopic_mode = SpectroscopicMode()
            wavmode = spectroscopic_mode(header=header)
        else:
            flat_grating = '_no_grating'
            wavmode = ''

        flat_slit = re.sub('[A-Za-z_ ]',
                           '',
                           header['SLIT'])

        filter2 = header['FILTER2']
        if filter2 == '<NO FILTER>':
            filter2 = ''
        else:
            filter2 = '_' + filter2

        master_flat_name += target_name \
                            + flat_grating \
                            + wavmode \
                            + filter2 \
                            + '_' \
                            + flat_slit \
                            + dome_sky \
                            + '.fits'

    elif technique == 'Imaging':
        flat_filter = re.sub('[- ]', '_', header['FILTER'])
        flat_filter = re.sub('[<> ]', '', flat_filter)
        master_flat_name += '_' + flat_filter + dome_sky + '.fits'

    return master_flat_name


def normalize_master_flat(master, name, method='simple', order=15):
    """Master flat normalization method

    This function normalize a master flat in three possible ways:
     *mean*: simply divide the data by its mean

     *simple*: Calculates the median along the spatial axis in order to obtain
     the dispersion profile. Then fits a
     :class:`~astropy.modeling.polynomial.Chebyshev1D` model and apply this to
     all the data.

     *full*: This is for experimental purposes only because it takes a lot of
     time to process. It will fit a model to each line along the dispersion axis
     and then divide it by the fitted model. I do not recommend this method
     unless you have a good reason as well as a very powerful computer.

    Args:
        master (CCDData): Master flat. Has to be a
          :class:`~astropy.nddata.CCDData` instance.
        name (str): Full path of master flat prior to normalization.
        method (str): Normalization method, 'mean', 'simple' or 'full'.
        order (int): Order of the polynomial to be fitted.

    Returns:
        master (CCDData):  The normalized master flat.
         :class:`~astropy.nddata.CCDData` instance.

    """
    assert isinstance(master, CCDData)
    master = master.copy()

    # define new name, base path and full new name
    new_name = 'norm_' + os.path.basename(name)
    path = os.path.dirname(name)
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
            master.header['GSP_NORM'] = ('full', 'Flat normalization method')

    # write normalized flat to a file
    write_fits(ccd=master,
               full_path=norm_name,
               parent_file=name)

    return master, norm_name


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

    # RIGHT ASCENSION conversion
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
    """Read fits files while adding important information to the header

    It is necessary to record certain data to the image header so that's the
    reason for this wrapper of :meth:`~astropy.nddata.CCDData.read` to exist.
    It will add the following keywords. In most cases, if the keyword already
    exist it will skip it except for `GSP_FNAM`, `GSP_PATH` and `BUNIT`.
    GSP_VERS: Goodman Spectroscopic Pipeline version number
    GSP_ONAM: Original File name
    GSP_PNAM: Parent file name or name of the file from which this one
    originated after some process or just a copy.
    GSP_FNAM: Current file name.
    GSP_PATH: Path to file at the moment of reading.
    GSP_TECH: Observing technique. `Spectroscopy` or `Imaging`.
    GSP_DATE: Date of first reading.
    GSP_OVER: Overscan region.
    GSP_TRIM: Trim section (region).
    GSP_SLIT: Slit trim section, obtained from the slit illuminated area.
    GSP_BIAS: Master bias image used. Default `none`.
    GSP_FLAT: Master flat image used. Default `none`.
    GSP_SCTR: Science target file
    GSP_NORM: Flat normalization method.
    GSP_COSM: Cosmic ray rejection method.
    GSP_EXTR: Extraction window at first column
    GSP_BKG1: First background extraction zone
    GSP_BKG2: Second background extraction zone
    GSP_WRMS: Wavelength solution RMS Error.
    GSP_WPOI: Number of points used to calculate the wavelength solution
    Error.
    GSP_WREJ: Number of points rejected.


    Args:
        full_path (str): Full path to file.
        technique (str): Observing technique. 'Imaging' or 'Spectroscopy'.

    Returns:
        Instance of :class:`~astropy.nddata.CCDData` corresponding to the file
          from `full_path`.

    """
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
    if 'GSP_PNAM' not in all_keys:
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

    if 'GSP_TMOD' not in all_keys:
        ccd.header.set('GSP_TMOD',
                       value='none',
                       comment='Model name used to fit trace')

    if 'GSP_EXTR' not in all_keys:
        ccd.header.set('GSP_EXTR',
                       value='none',
                       comment='Extraction window at first column')

    if 'GSP_BKG1' not in all_keys:
        ccd.header.set('GSP_BKG1',
                       value='none',
                       comment='First background extraction zone')

    if 'GSP_BKG2' not in all_keys:
        ccd.header.set('GSP_BKG2',
                       value='none',
                       comment='Second background extraction zone')

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

    return ccd


def recenter_broad_lines(lamp_data, lines, order):
    """Recenter broad lines

    Notes:
        This method is used to recenter broad lines only, there is a special
        method for dealing with narrower lines.

    Args:
        lamp_data (ndarray): numpy.ndarray instance. It contains the lamp
            data.
        lines (list): A line list in pixel values.
        order (float): A rough estimate of the FWHM of the lines in pixels
            in the data. It is calculated using the slit size divided by the
            pixel scale multiplied by the binning.

    Returns:
        A list containing the recentered line positions.

    """
    # TODO (simon): use slit size information for a square function
    # TODO (simon): convolution
    new_line_centers = []
    gaussian_kernel = Gaussian1DKernel(stddev=2.)
    lamp_data = convolve(lamp_data, gaussian_kernel)
    for line in lines:
        lower_index = max(0, int(line - order))
        upper_index = min(len(lamp_data), int(line + order))
        lamp_sample = lamp_data[lower_index:upper_index]
        x_axis = np.linspace(lower_index, upper_index, len(lamp_sample))
        line_max = np.max(lamp_sample)

        gaussian_model = models.Gaussian1D(amplitude=line_max,
                                           mean=line,
                                           stddev=order)

        fit_gaussian = fitting.LevMarLSQFitter()
        fitted_gaussian = fit_gaussian(gaussian_model, x_axis, lamp_sample)
        new_line_centers.append(fitted_gaussian.mean.value)

    return new_line_centers


def recenter_lines(data, lines, plots=False):
    """Finds the centroid of an emission line

    For every line center (pixel value) it will scan left first until the
    data stops decreasing, it assumes it is an emission line and then will
    scan right until it stops decreasing too. Defined those limits it will
    use the line data in between and calculate the centroid.

    Notes:
        This method is used to recenter relatively narrow lines only, there
        is a special method for dealing with broad lines.

    Args:
        data (ndarray): numpy.ndarray instance. or the data attribute of a
            :class:`~astropy.nddata.CCDData` instance.
        lines (list): A line list in pixel values.
        plots (bool): If True will plot spectral line as well as the input
            center and the recentered value.

    Returns:
        A list containing the recentered line positions.

    """
    new_center = []
    x_size = data.shape[0]
    median = np.median(data)
    for line in lines:
        # TODO (simon): Check if this definition is valid, so far is not
        # TODO (cont..): critical
        left_limit = 0
        right_limit = 1
        condition = True
        left_index = int(line)

        while condition and left_index - 2 > 0:

            if (data[left_index - 1] > data[left_index]) and \
                    (data[left_index - 2] > data[left_index - 1]):

                condition = False
                left_limit = left_index

            elif data[left_index] < median:
                condition = False
                left_limit = left_index

            else:
                left_limit = left_index

            left_index -= 1

        # id right limit
        condition = True
        right_index = int(line)
        while condition and right_index + 2 < x_size - 1:

            if (data[right_index + 1] > data[right_index]) and \
                    (data[right_index + 2] > data[right_index + 1]):

                condition = False
                right_limit = right_index

            elif data[right_index] < median:
                condition = False
                right_limit = right_index

            else:
                right_limit = right_index
            right_index += 1
        index_diff = [abs(line - left_index), abs(line - right_index)]

        sub_x_axis = range(line - min(index_diff),
                           (line + min(index_diff)) + 1)

        sub_data = data[line - min(index_diff):(line + min(index_diff)) + 1]
        centroid = np.sum(sub_x_axis * sub_data) / np.sum(sub_data)

        # checks for asymmetries
        differences = [abs(data[line] - data[left_limit]),
                       abs(data[line] - data[right_limit])]

        if max(differences) / min(differences) >= 2.:
            if plots:  # pragma: no cover
                plt.axvspan(line - 1, line + 1, color='g', alpha=0.3)
            new_center.append(line)
        else:
            new_center.append(centroid)
    if plots:  # pragma: no cover
        fig, ax = plt.subplots(1, 1)
        fig.canvas.set_window_title('Lines Detected in Lamp')
        ax.axhline(median, color='b')

        ax.plot(range(len(data)),
                data,
                color='k',
                label='Lamp Data')

        for line in lines:

            ax.axvline(line + 1,
                       color='k',
                       linestyle=':',
                       label='First Detected Center')

        for center in new_center:

            ax.axvline(center,
                       color='k',
                       linestyle='.-',
                       label='New Center')

        plt.show()
    return new_center


def record_trace_information(ccd, trace_info):
    """Adds trace information to fits header

    Notes:
        Example of trace_info.
        OrderedDict([('GSP_TMOD', ['Polynomial1D', 'Model name used to fit trace']),
                     ('GSP_TORD', [2, 'Degree of the model used to fit target trace']),
                     ('GSP_TC00', [80.92244303468138, 'Parameter c0']),
                     ('GSP_TC01', [0.0018921968204536187, 'Parameter c1']),
                     ('GSP_TC02', [-7.232545448865748e-07, 'Parameter c2']),
                     ('GSP_TERR', [0.18741058188097284, 'RMS error of target trace'])])

    Args:
        ccd (CCDData): CCDData instance to have trace info recorded into its
        header.
        trace_info (OrderedDict): Ordered Dictionary with a set of fits keywords
        associated to a list of values corresponding to value and comment.

    Returns:
        ccd (CCDData): Same CCDData instance with the header modified.
    """
    last_keyword = None
    for info_key in trace_info:
        info_value, info_comment = trace_info[info_key]
        log.debug(
            "Adding trace information: "
            "{:s} = {:s} / {:s}".format(info_key,
                                        str(info_value),
                                        info_comment))

        if last_keyword is None:
            ccd.header.set(info_key,
                           value=info_value,
                           comment=info_comment)
            last_keyword = info_key
        else:
            ccd.header.set(info_key,
                           value=info_value,
                           comment=info_comment,
                           after=last_keyword)
            last_keyword = info_key

    return ccd


def save_extracted(ccd, destination, prefix='e', target_number=1):
    """Save extracted spectrum while adding a prefix.

    Args:
        ccd (CCDData) :class:`~astropy.nddata.CCDData` instance
        destination (str): Path where the file will be saved.
        prefix (str): Prefix to be added to images. Default `e`.
        target_number (int): Secuential number of spectroscopic target.

    Returns:
        :class:`~astropy.nddata.CCDData` instance of the image just recorded.
          although is not really necessary.

    """
    assert isinstance(ccd, CCDData)
    assert os.path.isdir(destination)

    file_name = ccd.header['GSP_FNAM']
    if target_number > 0:
        new_suffix = '_target_{:d}.fits'.format(target_number)
        file_name = re.sub('.fits', new_suffix, file_name)

    if ccd.header['OBSTYPE'] in ['COMP', 'ARC']:
        extraction_region = re.sub(':','-', ccd.header['GSP_EXTR'])
        file_name = re.sub('.fits', '_{:s}.fits'.format(extraction_region),
                           file_name)
        new_file_name = prefix + file_name

    else:
        new_file_name = prefix + file_name

    log.info("Saving uncalibrated(w) extracted spectrum to file: "
             "{:s}".format(new_file_name))
    full_path = os.path.join(destination, new_file_name)
    ccd = write_fits(ccd=ccd, full_path=full_path, parent_file=file_name)
    return ccd


def search_comp_group(object_group, comp_groups, reference_data):
    """Search for a suitable comparison lamp group

    In case a science target was observed without comparison lamps, usually
    right before or right after, this function will look for a compatible set
    obtained at a different time or pointing.

    Notes:
        This methodology is not recommended for radial velocity studies.

    Args:
        object_group (DataFrame): A :class:`~pandas.DataFrame` instances
          containing a group of images for a given scientific target.
        comp_groups (list): A list in which every element is a
          :class:`~pandas.DataFrame`
          that contains information regarding groups of comparison lamps.
        reference_data (ReferenceData): Instance of
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
                log.debug('Found a matching comparison lamp group')
                return comp_group

    raise NoMatchFound


def setup_logging(debug=False, generic=False):  # pragma: no cover
    """configures logging

    Notes:
        Logging file name is set to default 'goodman_log.txt'.
        If --debug is activated then the format of the message is different.
    """

    log_filename = 'goodman_log.txt'

    if '--debug' in sys.argv or debug:
        log_format = '[%(asctime)s][%(levelname)8s]: %(message)s ' \
                     '[%(module)s.%(funcName)s:%(lineno)d]'
        logging_level = logging.DEBUG
    else:
        log_format = '[%(asctime)s][%(levelname).1s]: %(message)s'
        logging_level = logging.INFO

    date_format = '%H:%M:%S'

    formatter = logging.Formatter(fmt=log_format,
                                  datefmt=date_format)

    logging.basicConfig(level=logging_level,
                        format=log_format,
                        datefmt=date_format)

    log = logging.getLogger(__name__)

    file_handler = logging.FileHandler(filename=log_filename)
    file_handler.setFormatter(fmt=formatter)
    file_handler.setLevel(level=logging_level)
    log.addHandler(file_handler)

    if not generic:
        log.info("Starting Goodman HTS Pipeline Log")
        log.info("Local Time    : {:}".format(
            datetime.datetime.now()))
        log.info("Universal Time: {:}".format(
            datetime.datetime.utcnow()))

        try:
            latest_release = check_version.get_last()

            if "dev" in __version__:
                log.warning("Running Development version: {:s}".format(__version__))
                log.info("Latest Release: {:s}".format(latest_release))
            elif check_version.am_i_updated(__version__):
                if __version__ == latest_release:
                    log.info("Pipeline Version: {:s} (latest)".format(__version__))
                else:
                    log.warning("Current Version: {:s}".format(__version__))
                    log.info("Latest Release: {:s}".format(latest_release))
            else:
                log.warning("Current Version '{:s}' is outdated.".format(
                    __version__))
                log.info("Latest Release: {:s}".format(latest_release))
        except ConnectionRefusedError:
            log.error('Unauthorized GitHub API Access reached maximum')
            log.info("Current Version: {:s}".format(__version__))


def trace(ccd,
          model,
          trace_model,
          model_fitter,
          sampling_step,
          nfwhm=1,
          plots=False):
    """Find the trace of a spectrum

    This function is called by the `trace_targets` function, the difference is
    that it only takes single models only not `CompoundModels` so this function
    is called for every single target. `CompoundModels` are a bit tricky when
    you need each model separated so all `CompoundModels` have been removed.

    Notes:
        This method forces the trace to go withing a rectangular region of
        center `model.mean.value` and width `2 * nsigmas`, this is for allowing
        the tracing of low SNR targets. The assumption is valid since the
        spectra are always well aligned to the detectors's pixel columns.
        (dispersion axis)

    Args:
        ccd (CCDData) A :class:`~astropy.nddata.CCDData` instance, 2D image.
        model (Model): An astropy.modeling.Model instance that contains
          information regarding the target to be traced.
        trace_model (object): An astropy.modeling.Model instance, usually a low
          order polynomial.
        model_fitter (Fitter): An astropy.modeling.fitting.Fitter instance. Will
          fit the sampled points to construct the trace model
        sampling_step (int): Step for sampling the spectrum.
        nfwhm (int): Number of fwhm to each side of the mean to be used for
          searching the trace.
        plots (bool): Toggles debugging plot

    Returns:
        An `astropy.modeling.Model` instance, that defines the trace of the
        spectrum.

    """
    assert isinstance(ccd, CCDData)
    assert isinstance(model, Model)
    assert isinstance(trace_model, Model)

    spatial_length, dispersion_length = ccd.data.shape

    sampling_axis = range(0, dispersion_length, sampling_step)
    sample_values = []

    if model.__class__.name == 'Gaussian1D':
        model_fwhm = model.fwhm
        model_mean = model.mean.value
    elif model.__class__.name == 'Moffat1D':
        model_fwhm = model.fwhm
        model_mean = model.x_0.value
    else:
        raise NotImplementedError

    sample_center = float(model_mean)
    lower_limit_list = []
    upper_limit_list = []
    lower_limit = None
    upper_limit = None

    for point in sampling_axis:

        lower_limit = np.max([0, int(sample_center - nfwhm * model_fwhm)])
        upper_limit = np.min([int(sample_center + nfwhm * model_fwhm),
                              spatial_length])

        lower_limit_list.append(lower_limit)
        upper_limit_list.append(upper_limit)

        sample = ccd.data[lower_limit:upper_limit, point:point + sampling_step]
        sample_median = np.median(sample, axis=1)

        try:
            sample_peak = np.argmax(sample_median)
        except ValueError:  # pragma: no cover
            log.error('Nfwhm {}'.format(nfwhm))
            log.error('Model Stddev {}'.format(model_fwhm))
            log.error('sample_center {}'.format(sample_center))
            log.error('sample {}'.format(sample))
            log.error('sample_median {}'.format(sample_median))
            log.error('lower_limit {}'.format(lower_limit))
            log.error('upper_limit {}'.format(upper_limit))
            log.error('point {}'.format(point))
            log.error('point + sampling_step {}'.format(point + sampling_step))
            log.error("Spatial length: {}, Dispersion length {}".format(
                spatial_length,
                dispersion_length))
            sys.exit()

        sample_values.append(sample_peak + lower_limit)

        if np.abs(sample_peak + lower_limit - model_mean)\
                < nfwhm * model_fwhm:

            sample_center = int(sample_peak + lower_limit)

        else:
            sample_center = float(model_mean)

    trace_model.c2.fixed = True

    fitted_trace = model_fitter(trace_model, sampling_axis, sample_values)

    sampling_differences = [
        (fitted_trace(sampling_axis[i]) - sample_values[i]) ** 2
        for i in range(len(sampling_axis))]

    rms_error = np.sqrt(
        np.sum(np.array(sampling_differences))/len(sampling_differences))

    log.debug("RMS Error of unclipped trace differences {:.3f}".format(
        rms_error))

    clipped_values = sigma_clip(sampling_differences,
                                sigma=2,
                                maxiters=3,
                                cenfunc=np.ma.median)
    if np.ma.is_masked(clipped_values):
        _sampling_axis = list(sampling_axis)
        _sample_values = list(sample_values)

        sampling_axis = []
        sample_values = []
        for i in range(len(clipped_values)):
            if clipped_values[i] is not np.ma.masked:
                sampling_axis.append(_sampling_axis[i])
                sample_values.append(_sample_values[i])

        log.debug("Re-fitting the trace for a better trace.")

        trace_model.c2.fixed = False

        fitted_trace = model_fitter(trace_model, sampling_axis, sample_values)

        sampling_differences = [
            (fitted_trace(sampling_axis[i]) - sample_values[i]) ** 2 for i in
            range(len(sampling_axis))]

        rms_error = np.sqrt(
            np.sum(np.array(sampling_differences)) / len(sampling_differences))

        log.debug(
            "RMS Error after sigma-clipping trace differences {:.3f}".format(
                rms_error))

    trace_info = collections.OrderedDict()

    trace_info['GSP_TMOD'] = [fitted_trace.__class__.__name__,
                              'Model name used to fit trace']

    trace_info['GSP_TORD'] = [fitted_trace.degree,
                              'Degree of the model used to fit target trace']

    for i in range(fitted_trace.degree + 1):
        trace_info['GSP_TC{:02d}'.format(i)] = [
            fitted_trace.__getattribute__('c{:d}'.format(i)).value,
            'Parameter c{:d}'.format(i)]

    trace_info['GSP_TERR'] = [rms_error, 'RMS error of target trace']

    log.info("Target tracing RMS error: {:.3f}".format(rms_error))

    if plots:  # pragma: no cover
        z1 = np.mean(ccd.data) - 0.5 * np.std(ccd.data)
        z2 = np.median(ccd.data) + np.std(ccd.data)
        fig, ax = plt.subplots()
        fig.canvas.set_window_title(ccd.header['GSP_FNAM'])

        mng = plt.get_current_fig_manager()
        mng.window.showMaximized()

        ax.set_title("Tracing information\n{:s}\n"
                     "RMS Error {:.2f}".format(ccd.header['GSP_FNAM'],
                                               rms_error))
        ax.imshow(ccd.data, clim=(z1, z2), cmap='gray')
        ax.plot(sampling_axis,
                sample_values,
                color='b',
                marker='o',
                alpha=0.4,
                label='Sampling points')

        sampling_axis_limits = range(0, dispersion_length, sampling_step)

        low_span = fitted_trace(sampling_axis_limits) - (fitted_trace(sampling_axis_limits) - np.mean(lower_limit_list))
        up_span = fitted_trace(sampling_axis_limits) + (np.mean(upper_limit_list) - fitted_trace(sampling_axis_limits))

        ax.fill_between(sampling_axis_limits,
                        low_span,
                        up_span,
                        where=up_span > low_span,
                        facecolor='g',
                        interpolate=True,
                        alpha=0.3,
                        label='Aproximate extraction window')
        ax.plot(fitted_trace(range(dispersion_length)),
                color='r',
                linestyle='--',
                label='Fitted Trace Model')
        # plt.plot(model(range(spatial_length)))
        ax.legend(loc='best')
        plt.tight_layout()
        if plt.isinteractive():
            plt.draw()
            plt.pause(2)
        else:
            plt.show()

    return fitted_trace, trace_info


def trace_targets(ccd, target_list, sampling_step=5, pol_deg=2, nfwhm=5,
                  plots=False):
    """Find the trace of the target's spectrum on the image

    This function defines a low order polynomial that trace the location of the
    spectrum. The attributes pol_deg and sampling_step define the polynomial
    degree and the spacing in pixels for the samples. For every sample a
    gaussian model is fitted and the center (mean) is recorded and since
    spectrum traces vary smoothly this value is used as a new center for the
    base model used to fit the spectrum profile.

    Notes:
        This doesn't work for extended sources. Also this calls for the function
        `trace` for doing the actual trace, the difference is that this method
        is at a higher level.

    Args:
        ccd (CCDData) Instance of :class:`~astropy.nddata.CCDData`
        target_list (list): List of single target profiles.
        sampling_step (int): Frequency of sampling in pixels
        pol_deg (int): Polynomial degree for fitting the trace
        plots (bool): If True will show plots (debugging)
        nfwhm (int): Number of fwhm from spatial profile center to search for
        a target. default 10.

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

        single_trace, trace_info = trace(ccd=ccd,
                                         model=profile,
                                         trace_model=trace_model,
                                         model_fitter=model_fitter,
                                         sampling_step=sampling_step,
                                         nfwhm=nfwhm,
                                         plots=plots)

        if 0 < single_trace.c0.value < ccd.shape[0]:
            log.debug('Adding trace to list')
            all_traces.append([single_trace, profile, trace_info])
        else:
            log.error("Unable to trace target.")
            log.error('Trace is out of boundaries. Center: '
                      '{:.4f}'.format(single_trace.c0.value))

    if plots:  # pragma: no cover
        z1 = np.mean(ccd.data) - 0.5 * np.std(ccd.data)
        z2 = np.median(ccd.data) + np.std(ccd.data)
        fig, ax = plt.subplots()
        fig.canvas.set_window_title(ccd.header['GSP_FNAM'])

        mng = plt.get_current_fig_manager()
        mng.window.showMaximized()

        ax.set_title("Trace(s) for {:s}".format(ccd.header['GSP_FNAM']))
        ax.imshow(ccd.data, clim=(z1, z2), cmap='gray')
        ax.plot([], color='r', label='Trace(s)')
        for strace, prof, trace_info in all_traces:
            ax.plot(strace(range(ccd.data.shape[1])), color='r')
        ax.legend(loc='best')
        plt.tight_layout()
        plt.show()
    return all_traces


def validate_ccd_region(ccd_region, regexp='^\[\d*:\d*,\d*:\d*\]$'):
    compiled_reg_exp = re.compile(regexp)
    if not compiled_reg_exp.match(ccd_region):
        raise SyntaxError("ccd regions must be defined in the format "
                          "'[x1:x2,y1:y2]'")
    else:
        return True


def write_fits(ccd,
               full_path,
               combined=False,
               parent_file=None,
               overwrite=True):
    """Write fits while adding information to the header.

    This is a wrapper for allowing to save files while being able to add
    information into the header. Mostly for historical reasons.

    Args:
        ccd (CCDData) A :class:`~astropy.nddata.CCDData` instance to be saved
          to fits.
        full_path (str): Full path of file.
        combined (bool): True if `ccd` is the result of combining images.
        parent_file (str): Name of the file from which ccd originated. If
          combined is True this will be set to `combined`.
        overwrite (bool): Overwrite files, default True.

    Returns:
        :class:`~astropy.nddata.CCDData` instance.

    """
    assert isinstance(ccd, CCDData)
    if os.path.isabs(full_path) and not os.path.isdir(os.path.dirname(full_path)):
        log.error("Directory {} does not exist. Creating it right now."
                  "".format(os.path.dirname(full_path)))
        os.mkdir(os.path.dirname(full_path))

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
    ccd.header.set('GSP_PATH', value=os.path.dirname(full_path))

    # write to file
    log.info("Saving FITS file to {:s}".format(os.path.basename(full_path)))
    ccd.write(full_path, overwrite=overwrite)
    assert os.path.isfile(full_path)
    return ccd


# classes definition


class GenerateDcrParFile(object):
    """Creates dcr.par file based on lookup table

    `dcr` parameters depend heavily on binning, this class generates a file
    using the default format. The lookup table considers camera and binning.

    """
    _format = [
        "THRESH  = {:.1f} // Threshold (in STDDEV)",
        "XRAD    = {:d}   // x-radius of the box (size = 2 * radius)",
        "YRAD    = {:d}   // y-radius of the box (size = 2 * radius)",
        "NPASS   = {:d}   // Maximum number of cleaning passes",
        "DIAXIS  = {:d}   // Dispersion axis: 0 - no dispersion, 1 - X, 2 - Y",
        "LRAD    = {:d}   // Lower radius of region for replacement statistics",
        "URAD    = {:d}   // Upper radius of region for replacement statistics",
        "GRAD    = {:d}   // Growing radius",
        "VERBOSE = {:d}   // Verbose level [0,1,2]",
        "END"]

    _columns = ['parameter',
                'red-1',
                'red-2',
                'red-3',
                'blue-1',
                'blue-2',
                'blue-3']

    _lookup = [
        ['thresh', 3.0, 4.0, 3.0, 3.0, 3.0, 3.0],
        ['xrad', 9, 7, 9, 8, 9, 9],
        ['yrad', 9, 9, 9, 8, 9, 9],
        ['npass', 5, 5, 5, 5, 5, 5],
        ['diaxis', 0, 0, 0, 0, 0, 0],
        ['lrad', 1, 1, 1, 1, 1, 1],
        ['urad', 3, 3, 3, 3, 3, 3],
        ['grad', 1, 0, 1, 1, 1, 1],
        ['verbose', 1, 1, 1, 1, 1, 1]
    ]

    def __init__(self, par_file_name='dcr.par'):
        """

        Args:
            par_file_name:
        """
        self._file_name = par_file_name
        self._df = pandas.DataFrame(self._lookup, columns=self._columns)
        self._binning = "{:s}-{:s}"
        self._data_format = "\n".join(self._format)

    def __call__(self, instrument='Red', binning='1', path='default'):
        """

        Args:
            instrument (str): Instrument from INSTCONF keyword
            binning (str): Serial (dispersion) Binning from the header.
            path (str): Directory where to save the file.

        """
        assert any([instrument == option for option in ['Red', 'Blue']])
        b = self._binning.format(instrument.lower(), binning)
        self._data_format = self._data_format.format(
            self._df[b][self._df.parameter == 'thresh'].values[0],
            int(self._df[b][self._df.parameter == 'xrad'].values[0]),
            int(self._df[b][self._df.parameter == 'yrad'].values[0]),
            int(self._df[b][self._df.parameter == 'npass'].values[0]),
            int(self._df[b][self._df.parameter == 'diaxis'].values[0]),
            int(self._df[b][self._df.parameter == 'lrad'].values[0]),
            int(self._df[b][self._df.parameter == 'urad'].values[0]),
            int(self._df[b][self._df.parameter == 'grad'].values[0]),
            int(self._df[b][self._df.parameter == 'verbose'].values[0]))
        self._create_file(path=path)

    def _create_file(self, path):
        """Creates `dcr.par` file

        Args:
            path (str): Path to where to save the `dcr.par` file.

        """
        if os.path.isdir(path):
            full_path = os.path.join(path, self._file_name)
        else:
            full_path = os.path.join(os.getcwd(), self._file_name)

        with open(full_path, 'w') as dcr_par:
            dcr_par.write(self._data_format)


class NightDataContainer(object):
    """This class is designed to be the organized data container. It doesn't
    store image data but a list of :class:`~pandas.DataFrame` objects. Also it
    stores critical variables such as sunrise and sunset times.

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

        # comp_groups will store :class:`~pandas.DataFrame` (groups) that
        # contain only OBSTYPE == COMP, they should be requested only when
        # needed, for the science case when for every science target is
        # observed with comparison lamps and quartz (if)
        self.comp_groups = None

        # object_groups will store :class:`~pandas.DataFrame` (groups) with only
        # OBSTYPE == OBJECT this is the case when the observer takes comparison
        # lamps only at the beginning or end of the night.
        self.object_groups = None

        # spec_groups will store :class:`~pandas.DataFrame` (groups) with a set
        # of OBJECT and COMP, this is usually the case for radial velocity
        # studies.
        self.spec_groups = None

        """Time reference points"""
        self.sun_set_time = None
        self.sun_rise_time = None
        self.evening_twilight = None
        self.morning_twilight = None

    def __repr__(self):
        """Produces a nice summary of the information contained"""
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

            if all([self.gain, self.rdnoise, self.roi]):
                class_info += str("\nGain: {:.2f}\n"
                                  "Readout Noise: {:.2f}\n"
                                  "ROI: {:s}".format(self.gain,
                                                     self.rdnoise,
                                                     self.roi))
            class_info += str("\nIs Empty: {:s}\n".format(str(self.is_empty)))

            group_info = "\nData Grouping Information\n"
            group_info += "BIAS Group:\n"
            group_info += self._get_group_repr(self.bias)
            group_info += "Day FLATs Group:\n"
            group_info += self._get_group_repr(self.day_flats)
            group_info += "Dome FLATs Group:\n"
            group_info += self._get_group_repr(self.dome_flats)
            group_info += "Sky FLATs Group:\n"
            group_info += self._get_group_repr(self.sky_flats)
            if self.technique == 'Spectroscopy':

                group_info += "COMP Group:\n"
                group_info += self._get_group_repr(self.comp_groups)
                group_info += "OBJECT Group\n"
                group_info += self._get_group_repr(self.object_groups)
                group_info += "OBJECT + COMP Group:\n"
                group_info += self._get_group_repr(self.spec_groups)
                # group_info += self._get_group_repr(self.data_groups)

                class_info += group_info
            elif self.technique == 'Imaging':
                group_info += "DATA Group:\n"
                group_info += self._get_group_repr(self.data_groups)

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
                if len(group) == 1:
                    group_str += "Files in Group\n"
                else:
                    group_str += "Files in Group {:d}\n".format(i + 1)
                for _file in group[i]['file']:
                    group_str += "  {:s}\n".format(_file)
            return group_str
        else:

            return "  Group is Empty\n"

    def add_bias(self, bias_group):
        """Adds a bias group

        Args:
            bias_group (DataFrame): A :class:`~pandas.DataFrame` Contains a set
              of keyword values of grouped image metadata

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
            day_flats (DataFrame): A :class:`~pandas.DataFrame` Contains a set
              of keyword values of grouped image metadata

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
            data_group (DataFrame): A :class:`~pandas.DataFrame` Contains a set
              of keyword values of grouped image metadata

        """

        if self.data_groups is None:
            self.data_groups = [data_group]
        else:
            self.data_groups.append(data_group)
        if self.data_groups is not None:
            self.is_empty = False

    def add_comp_group(self, comp_group):
        """Adds a comp-only group

        All comparison lamps groups are added here. The ones that may have been
        taken in the afternoon (isolated) or along science target. This will
        act as a pool of comparison lamp groups for eventual science targets
        taken without comparison lamps.

        Args:
            comp_group (DataFrame): A :class:`~pandas.DataFrame` Contains a set
              of keyword values of grouped image metadata

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
            object_group (DataFrame): A :class:`~pandas.DataFrame` Contains a
              set of keyword values of grouped image metadata

        """

        if self.object_groups is None:
            self.object_groups = [object_group]
        else:
            self.object_groups.append(object_group)
        if self.object_groups is not None:
            self.is_empty = False

    def add_spec_group(self, spec_group):
        """Adds a data group containing object and comp

        The comparison lamp groups are also added to a general pool of
        comparison lamps.

        Args:
            spec_group (DataFrame): A :class:`~pandas.DataFrame` Contains a set
              of keyword values of grouped image metadata

        """

        if self.spec_groups is None:
            self.spec_groups = [spec_group]
        else:
            self.spec_groups.append(spec_group)
        if self.spec_groups is not None:
            self.is_empty = False
        comp_group = spec_group[spec_group.obstype == 'COMP']
        self.add_comp_group(comp_group=comp_group)

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
        """Set Gain, Read noise and ROI.

        Args:
            gain (float): Gain from header
            rdnoise (float): Read noise from header.
            roi (str): ROI from header.

        """
        self.gain = gain
        self.rdnoise = rdnoise
        self.roi = roi


class NoMatchFound(Exception):  # pragma: no cover
    """Exception for when no match is found."""
    def __init__(self, message="No match found"):
        Exception.__init__(self, message)


class NoTargetException(Exception):  # pragma: no cover
    """Exception to be raised when no target is identified"""
    def __init__(self):
        Exception.__init__(self, 'No targets identified.')


class NotEnoughLinesDetected(Exception):  # pragma: no cover
    """Exception for when there are no lines detected."""
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

        This methods uses ccdproc.ImageFileCollection on the `reference_dir` to
        capture all possible reference lamps. The reference lamps have a list
        of lines detected on the data registered to the header as GSP_P??? where
        ??? are numbers from 001 to 999. Also the pixel values are stored in
        keywords of the form GSP_A???.

        Args:
            reference_dir (str): full path to the reference data directory
        """
        self.log = logging.getLogger(__name__)
        self.reference_dir = reference_dir
        reference_collection = ccdproc.ImageFileCollection(self.reference_dir)
        self.ref_lamp_collection = reference_collection.summary.to_pandas()
        self.lines_pixel = None
        self.lines_angstrom = None
        self._ccd = None
        self.nist = {}
        self.lamp_status_keywords = [
            'LAMP_HGA',
            'LAMP_NE',
            'LAMP_AR',
            'LAMP_FE',
            'LAMP_CU',
            'LAMP_QUA',
            'LAMP_QPE',
            'LAMP_BUL',
            'LAMP_DOM',
            'LAMP_DPE']

    def get_reference_lamp(self, header):
        """Finds a suitable template lamp from the catalog

        Args:
            header (Header): FITS header of image we are looking a reference
                lamp.

        Returns:
            full path to best matching reference lamp.

        """

        if all([keyword in [hkey for hkey in header.keys()] for keyword in self.lamp_status_keywords]):
            self.log.info("Searching matching reference lamp")
            filtered_collection = self.ref_lamp_collection[(
                (self.ref_lamp_collection['lamp_hga'] == header['LAMP_HGA']) &
                (self.ref_lamp_collection['lamp_ne'] == header['LAMP_NE']) &
                (self.ref_lamp_collection['lamp_ar'] == header['LAMP_AR']) &
                (self.ref_lamp_collection['lamp_fe'] == header['LAMP_FE']) &
                (self.ref_lamp_collection['lamp_cu'] == header['LAMP_CU']) &
                (self.ref_lamp_collection['wavmode'] == header['wavmode']))]
            if filtered_collection.empty:
                error_message = "Unable to find a match for: "\
                                "LAMP_HGA = {}, "\
                                "LAMP_NE = {}, "\
                                "LAMP_AR = {}, "\
                                "LAMP_FE = {}, "\
                                "LAMP_CU = {}, "\
                                "WAVMODE = {} ".format(header['LAMP_HGA'],
                                                       header['LAMP_NE'],
                                                       header['LAMP_AR'],
                                                       header['LAMP_FE'],
                                                       header['LAMP_CU'],
                                                       header['WAVMODE'])
                self.log.error(error_message)
                raise NoMatchFound(error_message)
        else:
            filtered_collection = self.ref_lamp_collection[
                (self.ref_lamp_collection['object'] == header['object']) &
                # TODO (simon): Wavemode can be custom (GRT_TARG, CAM_TARG, GRATING)
                (self.ref_lamp_collection['wavmode'] == re.sub(' ', '_', header['wavmode']).upper())]
            if filtered_collection.empty:
                error_message = "Unable to find matching "\
                                "reference lamp for: "\
                                "OBJECT = {}, "\
                                "WAVMODE = {}".format(header['OBJECT'],
                                                      header['WAVMODE'])
                self.log.error(error_message)
                raise NoMatchFound(error_message)

        if len(filtered_collection) == 1:
            self.log.info(
                "Reference Lamp Found: {:s}"
                "".format(filtered_collection.file.to_string(index=False)))
            full_path = os.path.join(self.reference_dir,
                                     filtered_collection.file.to_string(
                                         index=False))
            self._ccd = CCDData.read(full_path, unit=u.adu)
            self._recover_lines()
            return self._ccd
        else:
            raise NotImplementedError(
                "Found {} matches".format(len(filtered_collection)))

    def lamp_exists(self, header):
        """Checks whether a matching lamp exist or not

        Args:
            object_name (str): Name of the lamp from 'OBJECT' keyword.
            grating (str): Grating from 'GRATING' keyword.
            grt_targ (float): Grating target from keyword 'GRT_TARG'.
            cam_targ (float): Camera target from keyword 'CAM_TARG'.

        Returns:
            True of False depending if a single matching lamp exist.

        Raises:
            NotImplementedError if there are more than one lamp found.

        """
        filtered_collection = self.ref_lamp_collection[
            (self.ref_lamp_collection['lamp_hga'] == header['LAMP_HGA']) &
            (self.ref_lamp_collection['lamp_ne'] ==  header['LAMP_NE']) &
            (self.ref_lamp_collection['lamp_ar'] ==  header['LAMP_AR']) &
            (self.ref_lamp_collection['lamp_cu'] ==  header['LAMP_CU']) &
            (self.ref_lamp_collection['lamp_fe'] ==  header['LAMP_FE']) &
            (self.ref_lamp_collection['grating'] ==  header['GRATING']) &
            (self.ref_lamp_collection['grt_targ'] == header['GRT_TARG']) &
            (self.ref_lamp_collection['cam_targ'] == header['CAM_TARG'])]

        if filtered_collection.empty:
            return False
        elif len(filtered_collection) == 1:
            return True
        else:
            raise NotImplementedError

    def check_comp_group(self, comp_group):
        """Check if comparison lamp group has matching reference lamps


        Args:
            comp_group (DataFrame): A :class:`~pandas.DataFrame` instance that
              contains meta-data for a group of comparison lamps.

        Returns:

        """
        lamps = comp_group.groupby(['grating',
                                    'grt_targ',
                                    'cam_targ',
                                    'lamp_hga',
                                    'lamp_ne',
                                    'lamp_ar',
                                    'lamp_fe',
                                    'lamp_cu']).size().reset_index(
        ).rename(columns={0: 'count'})

        # for the way the input is created this should run only once but the
        # for loop has been left in case this happens.
        for i in lamps.index:
            pseudo_header = fits.Header()

            # pseudo_header.set('OBJECT', value=lamps.iloc[i]['object'])
            pseudo_header.set('GRATING', value=lamps.iloc[i]['grating'])
            pseudo_header.set('GRT_TARG', value=lamps.iloc[i]['grt_targ'])
            pseudo_header.set('CAM_TARG', value=lamps.iloc[i]['cam_targ'])
            pseudo_header.set('LAMP_HGA', value=lamps.iloc[i]['lamp_hga'])
            pseudo_header.set('LAMP_NE', value=lamps.iloc[i]['lamp_ne'])
            pseudo_header.set('LAMP_AR', value=lamps.iloc[i]['lamp_ar'])
            pseudo_header.set('LAMP_FE', value=lamps.iloc[i]['lamp_fe'])
            pseudo_header.set('LAMP_CU', value=lamps.iloc[i]['lamp_cu'])

            if self.lamp_exists(header=pseudo_header):
                new_group = comp_group[
                    (comp_group['grating'] == lamps.iloc[i]['grating']) &
                    (comp_group['grt_targ'] == lamps.iloc[i]['grt_targ']) &
                    (comp_group['cam_targ'] == lamps.iloc[i]['cam_targ']) &
                    (comp_group['lamp_hga'] == lamps.iloc[i]['lamp_hga']) &
                    (comp_group['lamp_ne'] == lamps.iloc[i]['lamp_ne']) &
                    (comp_group['lamp_ar'] == lamps.iloc[i]['lamp_ar']) &
                    (comp_group['lamp_fe'] == lamps.iloc[i]['lamp_fe']) &
                    (comp_group['lamp_cu'] == lamps.iloc[i]['lamp_cu'])]
                return new_group
            else:
                self.log.warning("The target's comparison lamps do not have "
                                 "reference lamps.")
                self.log.debug("In this case a compatible lamp will be "
                               "obtained from all the lamps obtained in the "
                               "data or present in the files.")
                self.log.debug("Using the full set of comparison lamps "
                               "for extraction.")
                return comp_group
        return None

    def _recover_lines(self):
        """Read lines from the reference lamp's header."""
        self.log.info("Recovering line information from reference Lamp.")
        self.lines_pixel = []
        self.lines_angstrom = []
        pixel_keys = self._ccd.header['GSP_P*']
        for pixel_key in pixel_keys:
            if re.match(r'GSP_P\d{3}', pixel_key) is not None:
                angstrom_key = re.sub('GSP_P', 'GSP_A', pixel_key)
                assert pixel_key[-3:] == angstrom_key[-3:]
                assert angstrom_key in self._ccd.header
                if int(float(self._ccd.header[angstrom_key])) != 0:
                    self.lines_pixel.append(float(self._ccd.header[pixel_key]))
                    self.lines_angstrom.append(
                        float(self._ccd.header[angstrom_key]))
                else:
                    self.log.debug(
                        "File: {:s}".format(self._ccd.header['GSP_FNAM']))
                    self.log.debug(
                        "Ignoring keywords: {:s}={:f}, {:s}={:f}".format(
                            pixel_key,
                            self._ccd.header[pixel_key],
                            angstrom_key,
                            float(self._ccd.header[angstrom_key])))

    @staticmethod
    def _order_validation(lines_array):
        """Checks that the array of lines only increases."""
        previous = None
        for line_value in lines_array:
            if previous is not None:
                try:
                    assert line_value > previous
                    previous = line_value
                except AssertionError:
                    log.error("Error: Line {:f} is not larger "
                              "than {:f}".format(line_value, previous))
                    return False
            else:
                previous = line_value
        return True

    def _load_nist_list(self, **kwargs):
        """Load all csv files from strong lines in NIST."""
        nist_path = kwargs.get(
            'path',
            os.path.join(os.path.dirname(
                sys.modules['goodman_pipeline'].__file__),
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


class SaturationValues(object):
    """Contains a complete table of readout modes and 50% half well

    """

    def __init__(self, ccd=None):
        """Defines a :class:`~pandas.DataFrame` with saturation_threshold information

        Both, Red and Blue cameras have tabulated saturation_threshold values depending
        on the readout configurations. It defines a :class:`~pandas.DataFrame`
        object.

        Notes:
            For the purposes of this documentation *50% full well* is the same
            as ``saturation_threshold level`` though they are not the same thing.

        Args:
            ccd (CCDData): Image to be tested for saturation_threshold

        """
        self.log = logging.getLogger(__name__)
        self.__saturation = None
        columns = ['camera',
                   'read_rate',
                   'analog_attn',
                   'gain',
                   'read_noise',
                   'half_full_well',
                   'saturates_before']

        saturation_table = [['Blue', 50, 0, 0.25, 3.33, 279600, True],
                            ['Blue', 50, 2, 0.47, 3.35, 148723, True],
                            ['Blue', 50, 3, 0.91, 3.41, 76813, True],
                            ['Blue', 100, 0, 0.56, 3.69, 124821, True],
                            ['Blue', 100, 2, 1.06, 3.72, 65943, True],
                            ['Blue', 100, 3, 2.06, 3.99, 33932, False],
                            ['Blue', 200, 0, 1.4, 4.74, 49928, False],
                            ['Blue', 200, 2, 2.67, 5.12, 26179, False],
                            ['Blue', 400, 0, 5.67, 8.62, 12328, False],
                            ['Red', 100, 3, 1.54, 3.45, 66558, True],
                            ['Red', 100, 2, 3.48, 5.88, 29454, False],
                            ['Red', 344, 3, 1.48, 3.89, 69257, True],
                            ['Red', 344, 0, 3.87, 7.05, 26486, False],
                            ['Red', 750, 2, 1.47, 5.27, 69728, True],
                            ['Red', 750, 2, 1.45, 5.27, 69728, True],
                            ['Red', 750, 0, 3.77, 8.99, 27188, False],
                            ['Red', 750, 0, 3.78, 8.99, 27188, False]]
        self._sdf = pandas.DataFrame(saturation_table,
                                     columns=columns)

        if ccd is not None:
            self.get_saturation_value(ccd=ccd)

    @property
    def saturation_value(self):
        """Saturation value in counts

        In fact the value it returns is the 50% of full potential well,
        Some configurations reach digital saturation_threshold before 50% of full
        potential well, they are specified in the last column:
        ``saturates_before``.

        Returns:
            None if the value has not been defined

        """
        if self.__saturation is None:
            self.log.error('Saturation value not set')
            return None
        else:
            return self.__saturation

    def get_saturation_value(self, ccd):
        """Defines the saturation_threshold level

        Args:
            ccd (CCDData): Image to be tested for saturation_threshold

        Returns:
            The saturation_threshold value or None

        """
        hfw = self._sdf.half_full_well[
            (self._sdf.camera == ccd.header['INSTCONF']) &
            (self._sdf.gain == ccd.header['GAIN']) &
            (self._sdf.read_noise == ccd.header['RDNOISE'])]

        if hfw.empty:
            self.log.critical('Unable to obtain saturation_threshold level')
            self.__saturation = None
            return None
        else:
            self.__saturation = float(hfw.to_string(index=False))
            self.log.debug("Set saturation_threshold level as {:.0f}".format(
                self.__saturation))
            return self.__saturation


class SpectroscopicMode(object):

    def __init__(self):
        """Init method for the Spectroscopic Mode

        This method defines a :class:`~pandas.DataFrame` instance that contains
        all the current standard wavelength modes for Goodman HTS.

        """
        self.log = logging.getLogger(__name__)
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
            header (Header): FITS header.
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

            return self.get_mode(grating=grating,
                                 camera_targ=camera_targ,
                                 grating_targ=grating_targ,
                                 blocking_filter=blocking_filter)
        elif not all(x is None for x in (
                grating, camera_targ, grating_targ, blocking_filter)):
            grating = re.sub('[A-Za-z_-]', '', grating)

            return self.get_mode(grating=grating,
                                 camera_targ=camera_targ,
                                 grating_targ=grating_targ,
                                 blocking_filter=blocking_filter)
        else:
            raise SyntaxError("Either a fits header or grating, camera angle, "
                              "grating angle and order blocking filter are "
                              "required.")

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
        if any(grat == grating for grat in ('1800', '2100', '2400')):
            central_wavelength = get_central_wavelength(grating=grating,
                                                        grt_ang=grating_targ,
                                                        cam_ang=camera_targ)
            central_wavelength.to(u.nm)
            return 'Custom_{:d}nm'.format(int(round(central_wavelength.value)))

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
                central_wavelength.to(u.nm)
                return 'Custom_{:d}nm'.format(int(round(
                    central_wavelength.value)))
            else:
                return _mode['wavmode'].to_string(index=False)

    def get_cam_grt_targ_angle(self, grating, mode):
        """Get the camera and grating target values grating and mode

        Args:
            grating (float): Grating frequency in lines/mm (unitless value)
            mode (str): Name of the grating's mode for which the camera and
              grating target values are required.

        Returns:
            Camera and grating target values. None and None if no such values
            exists.

        """
        if any(grat == str(grating) for grat in ('1800', '2100', '2400')):
            self.log.warning("Grating {:s} does not define "
                             "modes.".format(str(grating)))
            return None, None
        else:
            angle = self.modes_data_frame[
                ((self.modes_data_frame['grating_freq'] == str(grating)) &
                 (self.modes_data_frame['wavmode'] == mode))]
            if angle.empty:
                self.log.error("No data")
                return None, None
            else:
                return (angle['camtarg'].to_string(index=False),
                        angle['grttarg'].to_string(index=False))


class IdentifySpectroscopicTargets(object):

    def __init__(self):
        self.nfind = 1
        self.plots = False
        self.background_threshold = 3
        self.profile_model = []
        self.model_name = None
        self.ccd = None
        self.slit_size = None
        self.serial_binning = None
        self.order = None
        self.file_name = None
        self.background_model = None
        self.background_level = None
        self.spatial_profile = None
        self.all_peaks = None
        self.selected_peaks = None

    def __call__(self,
                 ccd,
                 nfind=3,
                 background_threshold=3,
                 model_name='gaussian',
                 plots=False):
        assert isinstance(ccd, CCDData)
        assert ccd.header['OBSTYPE'] in ['OBJECT', 'SPECTRUM'], \
            "Can't search for targets in files with" \
            " OBSTYPE = {}".format(ccd.header['OBSTYPE'])
        self.file_name = ccd.header['GSP_FNAM']

        log.info('Searching spectroscopic targets in file: {:s}'
                 ''.format(self.file_name))

        self.ccd = ccd
        self.nfind = nfind
        self.plots = plots
        self.model_name = model_name
        self.background_threshold = background_threshold

        self.slit_size = re.sub('[a-zA-Z"_*]', '', self.ccd.header['SLIT'])
        log.debug('Slit size: {:s}'.format(self.slit_size))
        self.serial_binning = int(self.ccd.header['CCDSUM'].split()[0])
        log.debug('Serial binning: {:d}'.format(self.serial_binning))

        self.order = int(round(float(self.slit_size) / (0.15 * self.serial_binning)))

        if self.plots:  # pragma: no cover
            z1 = np.mean(self.ccd.data) - 0.5 * np.std(self.ccd.data)
            z2 = np.median(self.ccd.data) + np.std(self.ccd.data)

            plt.switch_backend('Qt5Agg')
            fig, ax = plt.subplots()
            fig.canvas.set_window_title(self.file_name)

            mng = plt.get_current_fig_manager()
            mng.window.showMaximized()

            ax.set_title(self.file_name)
            ax.imshow(ccd.data, clim=(z1, z2), cmap='gray')
            ax.set_xlabel('Dispersion Axis (x)')
            ax.set_ylabel('Spatial Axis (y)')
            fig.tight_layout()
            plt.show()

        self.spatial_profile = np.median(ccd.data, axis=1)

        # assert all([self.spatial_profile, self.file_name])

        self.fit_background()

        self.subtract_background()

        self.get_peaks()

        self.filter_peaks()

        self.fit_model()

        if self.profile_model == []:
            log.error("Impossible to identify targets.")
        else:
            log.info('Identified {:d} target{:s}'.format(
                len(self.profile_model),
                ['s' if len(self.profile_model) > 1 else ''][0]))

        return self.profile_model

    def fit_background(self, spatial_profile=None, file_name=None, plots=False):
        """

        Args:
            spatial_profile :
            file_name (String):
            plots:

        Returns:

        """
        if spatial_profile is None and self.spatial_profile is not None:
            spatial_profile = self.spatial_profile
        else:
            raise NotImplementedError
        if file_name is None and self.file_name is not None:
            file_name = self.file_name
        else:
            raise NotImplementedError

        log.info('Fitting Linear1D model to spatial profile to detect '
                 'background shape')

        clipped_profile = sigma_clip(spatial_profile, sigma=2, maxiters=5)

        linear_model = models.Linear1D(slope=0,
                                       intercept=np.median(spatial_profile))

        linear_fitter = fitting.LinearLSQFitter()

        # the fitters do not support masked arrays so we need to have a new
        # array without the masked (clipped) elements.
        new_profile = clipped_profile[~clipped_profile.mask]

        # also the indexes are different
        new_x_axis = np.array([i for i in range(len(clipped_profile))
                               if not clipped_profile.mask[i]])

        self.background_model = linear_fitter(linear_model, new_x_axis, new_profile)

        if plots or self.plots:  # pragma: no cover
            fig, ax = plt.subplots()
            fig.canvas.set_window_title(file_name)

            mng = plt.get_current_fig_manager()
            mng.window.showMaximized()

            ax.set_title('Background Fitting Model Defined')
            ax.plot(spatial_profile, color='k', label='Median profile')
            ax.plot(linear_model(range(len(spatial_profile))),
                    color='r',
                    label='Background Linear Model')
            ax.set_xlabel("Spatial Axis (Pixels)")
            ax.set_ylabel("Median Intensity")
            ax.legend(loc='best')
            plt.tight_layout()
            plt.show()

            fig, ax = plt.subplots()
            fig.canvas.set_window_title(file_name)

            mng = plt.get_current_fig_manager()
            mng.window.showMaximized()

            ax.set_title('Background Fitted Model')
            ax.plot(spatial_profile, color='k', label='Median profile')
            ax.plot(self.background_model(range(len(spatial_profile))),
                    color='r',
                    label='Fitted Background Linear Model')
            ax.set_xlabel("Spatial Axis (Pixels)")
            ax.set_ylabel("Median Intensity")
            ax.legend(loc='best')
            plt.tight_layout()
            plt.show()

        return self.background_model

    def subtract_background(self, spatial_profile=None, background_model=None,
                                     file_name=None, plots=False):
        """

        Args:
            spatial_profile:
            background_model:
            file_name:
            plots:

        Returns:

        """

        if not all([spatial_profile, background_model, file_name]):
            if self.spatial_profile is not None:
                spatial_profile = self.spatial_profile
            if self.background_model is not None:
                background_model = self.background_model
            if self.file_name is None:
                file_name = ''
            else:
                file_name = self.file_name

        log.info('Subtracting background shape and level spatial profile for '
                 'better target identification')

        background_array = background_model(range(len(spatial_profile)))

        background_subtracted = spatial_profile - background_array

        background_subtracted[background_subtracted < 0] = 0

        self.spatial_profile = background_subtracted.copy()

        clipped_final_profile = sigma_clip(self.spatial_profile, sigma=3, maxiters=3)

        new_x_axis = [i for i in range(len(clipped_final_profile)) if
                      not clipped_final_profile.mask[i]]

        clipped_final_profile = clipped_final_profile[
            ~clipped_final_profile.mask]

        self.background_level = np.abs(np.max(clipped_final_profile) -
                                       np.min(clipped_final_profile))
        log.debug('New background level after subtraction was found to be '
                  '{:.2f}'.format(self.background_level))

        if plots or self.plots:  # pragma: no cover
            plt.ioff()
            plt.close()

            fig, ax = plt.subplots()
            fig.canvas.set_window_title(file_name)

            mng = plt.get_current_fig_manager()
            mng.window.showMaximized()

            ax.set_title('Median Along Dispersion Axis (spatial)')
            ax.plot(background_subtracted, label='Background Subtracted Data')
            ax.plot(new_x_axis,
                    clipped_final_profile,
                    color='r',
                    label='Sigma Clipped Data')

            ax.axhline(self.background_level, color='m', label='Min-Max Difference')
            ax.set_xlabel("Spatial Axis (Pixels)")
            ax.set_ylabel("Median Intensity")
            plt.legend(loc='best')
            plt.tight_layout()
            if plt.isinteractive():
                plt.draw()
                plt.pause(5)
            else:
                plt.show()

        return self.spatial_profile, self.background_level

    def get_peaks(self,
                  spatial_profile=None,
                  order=None,
                  file_name=None,
                  plots=False):
        """

        Args:
            spatial_profile: Background subtracted profile
            order:
            file_name:
            plots:

        Returns:

        """

        if not all([spatial_profile, order, file_name]):
            if self.spatial_profile is not None:
                spatial_profile = self.spatial_profile
            if self.order is not None:
                order = self.order
            if self.file_name is None:
                file_name = ''
            else:
                file_name = self.file_name

        log.info("Finding all peaks in spatial profile")

        spatial_profile = signal.medfilt(spatial_profile, kernel_size=1)
        _upper_limit = spatial_profile.min() + 0.03 * spatial_profile.max()

        filtered_profile = np.where(np.abs(
            spatial_profile > spatial_profile.min() + 0.03 * spatial_profile.max()),
            spatial_profile,
            None)

        none_to_zero_prof = [0 if it is None else it for it in filtered_profile]

        filtered_profile = np.array(none_to_zero_prof)

        # order *= 2
        self.all_peaks = signal.argrelmax(filtered_profile,
                                          axis=0,
                                          order=order)[0]
        log.debug("Found {:d} peaks".format(len(self.all_peaks)))

        if plots or self.plots:  # pragma: no cover
            plt.ioff()

            fig, ax = plt.subplots()
            fig.canvas.set_window_title(file_name)

            ax.set_title('All detected Peaks')

            mng = plt.get_current_fig_manager()
            mng.window.showMaximized()
            for peak in self.all_peaks:
                ax.axvline(peak, color='r', alpha=0.7)
            ax.plot(spatial_profile, label='Background subtracted profile')
            ax.axhline(_upper_limit, color='g', label='Peak Detection Threshold')
            ax.plot([], color='r', label='Peak location')

            ax.set_xlabel("Spatial Axis (Pixels)")
            ax.set_ylabel("Background subtracted median intensity")
            ax.legend(loc='best')
            plt.tight_layout()
            plt.show()

        return self.all_peaks

    def filter_peaks(self,
                     spatial_profile=None,
                     detected_peaks=None,
                     nfind=None,
                     background_threshold=None,
                     file_name=None,
                     plots=False):
        """

        Args:
            spatial_profile:
            detected_peaks:
            nfind:
            background_threshold:
            file_name:
            plots:

        Returns:

        """
        if not all([spatial_profile,
                    detected_peaks,
                    nfind,
                    background_threshold,
                    file_name]):

            if self.spatial_profile is not None:
                spatial_profile = self.spatial_profile
            if self.all_peaks is not None:
                detected_peaks = self.all_peaks
            if self.nfind is not None:
                nfind = self.nfind
            if self.background_threshold is not None:
                background_threshold = self.background_threshold
            if self.file_name is None:
                file_name = ''
            else:
                file_name = self.file_name
        else:
            raise NotImplementedError

        log.info("Selecting the {:d} most intense peaks out of {:d} found"
                 "".format(nfind, len(detected_peaks)))

        peak_data_values = [spatial_profile[i] for i in detected_peaks]

        sorted_values = np.sort(peak_data_values)[::-1]

        detection_limit = spatial_profile.min() + 0.03 * spatial_profile.max()

        n_strongest_values = sorted_values[:nfind]

        self.selected_peaks = []
        log.info("Validating peaks by setting threshold {:d} times the "
                 "background level {:.2f}".format(background_threshold,
                                                  detection_limit))
        log.debug('Intensity threshold set to: {:.2f}'
                  ''.format(background_threshold * detection_limit))
        for peak_value in n_strongest_values:
            index = np.where(peak_data_values == peak_value)[0]
            if peak_value > background_threshold * detection_limit:
                self.selected_peaks.append(detected_peaks[index[0]])
                log.info(
                    'Selecting peak: Centered: {:.1f} Intensity {:.3f}'.format(
                        self.selected_peaks[-1], peak_value))
            else:
                log.debug('Discarding peak: Center {:.1f} Intensity {:.3f} '
                          'Reason: Below intensity threshold ({:.2f})'
                          ''.format(detected_peaks[index[0]],
                                    peak_value,
                                    background_threshold * detection_limit))

        if plots or self.plots:  # pragma: no cover
            plt.ioff()

            fig, ax = plt.subplots()
            fig.canvas.set_window_title(file_name)

            mng = plt.get_current_fig_manager()
            mng.window.showMaximized()

            ax.plot(spatial_profile, label='Background subtracted profile')
            ax.axhline(detection_limit, color='g', label='Upper limit for peak detection')
            ax.axhline(background_threshold * detection_limit,
                       color='m',
                       label="Intensity Threshold")
            for peak in self.selected_peaks:
                ax.axvline(peak, color='r', label='Peak location')
            ax.set_xlabel("Spatial Axis (Pixels)")
            ax.set_ylabel("Background subtracted median intensity")
            ax.legend(loc='best')
            plt.tight_layout()
            plt.show()

        return self.selected_peaks

    def fit_model(self,
                  spatial_profile=None,
                  selected_peaks=None,
                  order=None,
                  model_name=None,
                  file_name=None,
                  plots=False):
        if not all([spatial_profile,
                    selected_peaks,
                    order,
                    model_name,
                    file_name]):

            if self.spatial_profile is not None:
                spatial_profile = self.spatial_profile
            if self.all_peaks is not None:
                selected_peaks = self.selected_peaks
            if self.order is not None:
                order = self.order
            if self.model_name is not None:
                model_name = self.model_name
            if self.file_name is None:
                file_name = ''
            else:
                file_name = self.file_name
        else:
            raise NotImplementedError

        fitter = fitting.LevMarLSQFitter()

        if model_name == 'gaussian':
            self.profile_model = self._fit_gaussian(
                fitter=fitter,
                spatial_profile=spatial_profile,
                selected_peaks=selected_peaks,
                order=order,
                file_name=file_name,
                plots=plots or self.plots)
            return self.profile_model

        if model_name == 'moffat':
            self.profile_model = self._fit_moffat(
                fitter=fitter,
                spatial_profile=spatial_profile,
                selected_peaks=selected_peaks,
                order=order,
                file_name=file_name,
                plots=plots or self.plots)
            return self.profile_model

    @staticmethod
    def _fit_gaussian(fitter,
                      spatial_profile,
                      selected_peaks,
                      order,
                      file_name,
                      plots):
        log.info("Fitting 'Gaussian1D' to spatial profile of targets.")
        profile_model = []
        for peak in selected_peaks:
            peak_value = spatial_profile[peak]
            gaussian = models.Gaussian1D(amplitude=peak_value,
                                         mean=peak,
                                         stddev=order).rename(
                'Gaussian_{:}'.format(peak))

            fitted_gaussian = fitter(gaussian,
                                     range(len(spatial_profile)),
                                     spatial_profile)

            # this ensures the profile returned are valid
            if (fitted_gaussian.stddev.value > 0) and \
                    (fitted_gaussian.stddev.value < 4 * order):
                profile_model.append(fitted_gaussian)
                log.info(
                    "Recording target centered at: {:.2f}, stddev: {:.2f}"
                    "".format(fitted_gaussian.mean.value,
                              fitted_gaussian.stddev.value))
            else:
                log.error("Discarding target with stddev: {:.3f}".format(
                    fitted_gaussian.stddev.value))

        if plots:  # pragma: no cover
            fig, ax = plt.subplots()
            fig.canvas.set_window_title(file_name)

            mng = plt.get_current_fig_manager()
            mng.window.showMaximized()
            ax.set_title('Successfully fitted profiles')

            ax.plot(spatial_profile, color='k', label='Median Profile')
            for profile in profile_model:
                ax.plot(profile(range(len(spatial_profile))),
                        label=profile.name)

            ax.set_xlabel("Spatial Axis (Pixels)")
            ax.set_ylabel("Median Intensity")
            ax.legend(loc='best')
            plt.tight_layout()
            plt.show()

        return profile_model

    @staticmethod
    def _fit_moffat(fitter,
                    spatial_profile,
                    selected_peaks,
                    order,
                    file_name,
                    plots):
        log.info("Fitting 'Moffat1D' to spatial profile of targets.")
        profile_model = []
        for peak in selected_peaks:
            peak_value = spatial_profile[peak]
            moffat = models.Moffat1D(amplitude=peak_value,
                                       x_0=peak,
                                       gamma=order).rename(
                'Moffat_{:}'.format(peak))

            fitted_moffat = fitter(moffat,
                                   range(len(spatial_profile)),
                                   spatial_profile)

            # this ensures the profile returned are valid
            if (fitted_moffat.fwhm > 0.5 * order) and \
                    (fitted_moffat.fwhm < 4 * order):
                profile_model.append(fitted_moffat)
                log.info(
                    "Recording target centered at: {:.2f}, fwhm: {:.2f}"
                    "".format(fitted_moffat.x_0.value,
                              fitted_moffat.fwhm))
            else:
                log.error("Discarding target centered at: {:.3f}".format(
                    fitted_moffat.x_0.value))
                if fitted_moffat.fwhm < 0:
                    log.error("Moffat model FWHM is negative")
                elif 0 <= fitted_moffat.fwhm < 0.5 * order:
                    log.error("Moffat model FWHM is too small: {:.3f}, most "
                              "likely is an artifact".format(fitted_moffat.fwhm))
                else:
                    log.error("Moffat model FWHM too large: {:.3f}"
                              "".format(fitted_moffat.fwhm))

        if plots:  # pragma: no cover
            fig, ax = plt.subplots()
            fig.canvas.set_window_title(file_name)

            mng = plt.get_current_fig_manager()
            mng.window.showMaximized()
            ax.set_title('Successfully fitted profiles')

            ax.plot(spatial_profile, color='k', label='Median Profile')
            for profile in profile_model:
                ax.plot(profile(range(len(spatial_profile))),
                        label=profile.name)

            ax.set_xlabel("Spatial Axis (Pixels)")
            ax.set_ylabel("Median Intensity")
            ax.legend(loc='best')
            plt.tight_layout()
            plt.show()

        return profile_model
