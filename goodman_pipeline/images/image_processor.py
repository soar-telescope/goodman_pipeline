from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import ccdproc
import datetime
import glob
import logging
import numpy as np
import matplotlib.pyplot as plt
import random
import re
import os

from astropy import units as u
from ccdproc import CCDData
from ..core import (astroscrappy_lacosmic,
                    call_cosmic_rejection,
                    combine_data,
                    get_best_flat,
                    get_slit_trim_section,
                    image_overscan,
                    image_trim,
                    normalize_master_flat,
                    read_fits,
                    write_fits)

from ..core import SaturationValues, SpectroscopicMode


class ImageProcessor(object):
    """Image processing class

    This class contains methods for performing CCD image reduction for
    spectroscopy and imaging.

    """

    def __init__(self, args, data_container):
        """Initialization method for ImageProcessor class

        Args:
            args (Namespace): argparse instance
            data_container (DataFrame): Contains relevant information of the
            night and the data itself.
        """
        # TODO (simon): Check how inheritance could be used here.
        self.log = logging.getLogger(__name__)
        self.args = args
        self.instrument = data_container.instrument
        self.technique = data_container.technique
        self.gain = data_container.gain
        self.rdnoise = data_container.rdnoise
        self.roi = data_container.roi
        self.bias = data_container.bias
        self.day_flats = data_container.day_flats
        self.dome_flats = data_container.dome_flats
        self.sky_flats = data_container.sky_flats
        self.data_groups = data_container.data_groups
        self.spec_groups = data_container.spec_groups
        self.sun_set = data_container.sun_set_time
        self.sun_rise = data_container.sun_rise_time
        self.morning_twilight = data_container.morning_twilight
        self.evening_twilight = data_container.evening_twilight
        self.pixel_scale = 0.15 * u.arcsec
        self.queue = None
        self.trim_section = self.define_trim_section(
            technique=self.technique)
        self.overscan_region = self.get_overscan_region()
        self.saturation_values = SaturationValues()
        self.spec_mode = SpectroscopicMode()
        self.master_bias = None
        self.master_bias_name = None
        self.out_prefix = None

    def __call__(self):
        """Call method for ImageProcessor class

        This method manages the image processing by calling the appropriate
        methods.

        """

        for group in [self.bias,
                      self.day_flats,
                      self.dome_flats,
                      self.sky_flats,
                      self.data_groups,
                      self.spec_groups]:
            if group is not None:
                for sub_group in group:
                    group_obstype = sub_group.obstype.unique()

                    if len(group_obstype) == 1 and \
                        group_obstype[0] == 'BIAS' and \
                            not self.args.ignore_bias:

                        self.log.debug('Creating Master Bias')
                        self.create_master_bias(sub_group)
                    elif len(group_obstype) == 1 and group_obstype[0] == 'FLAT':
                        self.log.debug('Create Master FLATS')
                        self.create_master_flats(sub_group)
                    else:
                        self.log.debug('Process Data Group')
                        if self.technique == 'Spectroscopy':
                            self.process_spectroscopy_science(sub_group)
                        else:
                            self.log.info('Processing Imaging Science Data')
                            self.process_imaging_science(sub_group)

    def _is_file_saturated(self, ccd):
        """Detects a saturated image

        It counts the number of pixels above the saturation level, then finds
        which percentage they represents and if it is above the threshold it
        will return True. The percentage threshold can be set using the command
        line argument ``--saturation``.

        Args:
            ccd (CCDData): Image to be tested for saturation

        Returns:
            True for saturated and False for non-saturated

        """

        pixels_above_saturation = np.count_nonzero(
            ccd.data[np.where(
                ccd.data > self.saturation_values.get_saturation_value(
                    ccd=ccd))])

        total_pixels = np.count_nonzero(ccd.data)

        saturated_percent = (pixels_above_saturation * 100.) / total_pixels

        if saturated_percent >= float(self.args.saturation_threshold):
            self.log.warning(
                "The current image has more than {:.2f} percent "
                "of pixels above saturation level".format(
                float(self.args.saturation_threshold)))
            return True
        else:
            return False

    def define_trim_section(self, technique=None):
        """Get the initial trim section

        The initial trim section is usually defined in the header with the
        keyword ``TRIMSEC`` but in the case of Goodman HTS this does not work well.
        In particular for spectroscopy where is more likely to have combined
        binning and so on.

        Args:
            technique (str): The name of the technique, the options are:
                Imaging or Spectroscopy.

        Returns:
            The trim section in the format ``[x1:x2, y1:y2]``

        """

        assert technique is not None
        trim_section = None
        # TODO (simon): Consider binning and possibly ROIs for trim section
        self.log.warning('Determining trim section. Assuming you have only one '
                         'kind of data in this folder')
        for group in [self.bias,
                      self.day_flats,
                      self.dome_flats,
                      self.sky_flats,
                      self.data_groups]:
            if group is not None:
                # print(self.bias[0])
                image_list = group[0]['file'].tolist()
                sample_image = os.path.join(self.args.raw_path,
                                            random.choice(image_list))
                ccd = read_fits(sample_image, technique=technique)

                # serial binning - dispersion binning
                # parallel binngin - spatial binning
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

                    #top
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

                self.log.info('Trim Section: %s', trim_section)
                return trim_section

    def get_overscan_region(self):
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


        Returns:
            overscan_region (str) Region for overscan in the format
              '[min:max,:]' where min is the starting point and max is the end
               point of the overscan region.

        """
        for group in [self.bias,
                      self.day_flats,
                      self.dome_flats,
                      self.sky_flats,
                      self.data_groups]:
            if group is not None:
                # 'group' is a list
                image_list = group[0]['file'].tolist()
                sample_image = os.path.join(self.args.raw_path,
                                            random.choice(image_list))
                self.log.debug('Overscan Sample File ' + sample_image)
                ccd = CCDData.read(sample_image, unit=u.adu)

                # Image height - spatial direction
                spatial_length, dispersion_length = ccd.data.shape

                # Image width - spectral direction
                # w = ccd.data.shape[1]

                # Take the binnings
                serial_binning, parallel_binning = \
                    [int(x) for x in ccd.header['CCDSUM'].split()]

                if self.technique == 'Spectroscopy':
                    self.log.info('Overscan regions has been tested for ROI '
                                  'Spectroscopic 1x1, 2x2 and 3x3')

                    # define l r b and t to avoid local variable might be
                    # defined before assignment warning
                    low_lim_spectral,\
                        high_lim_spectral,\
                        low_lim_spatial,\
                        high_lim_spatial = [None] * 4
                    if self.instrument == 'Red':
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
                    elif self.instrument == 'Blue':
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

                elif self.technique == 'Imaging':
                    self.log.warning("Imaging mode doesn't have overscan "
                                     "region. Use bias instead.")
                    if self.bias is None:
                        self.log.warning('Bias are needed for Imaging mode')
                    overscan_region = None
                else:
                    overscan_region = None
                self.log.info('Overscan Region: %s', overscan_region)
                return overscan_region

    def create_master_bias(self, bias_group):
        """Create Master Bias

        Given a :class:`~pandas.DataFrame` object that contains a list of compatible bias.
        This function creates the master flat using ccdproc.combine using median
        and 3-sigma clipping.

        Args:
            bias_group (object): :class:`~pandas.DataFrame` instance that contains a list
                of bias images compatible with each other.

        """
        bias_file_list = bias_group.file.tolist()
        default_bias_name = os.path.join(self.args.red_path, 'master_bias.fits')
        search_bias_name = re.sub('.fits', '*.fits', default_bias_name)
        n_bias = len(glob.glob(search_bias_name))
        if n_bias > 0:
            self.master_bias_name = re.sub('.fits',
                                           '_{:d}.fits'.format(n_bias + 1),
                                           default_bias_name)

            self.log.info('New name for master bias: ' + self.master_bias_name)
        else:
            self.master_bias_name = default_bias_name

        master_bias_list = []
        self.log.info('Creating master bias')
        for image_file in bias_file_list:
            image_full_path = os.path.join(self.args.raw_path, image_file)
            ccd = read_fits(image_full_path, technique=self.technique)
            self.log.debug('Loading bias image: ' + image_full_path)
            if self.technique == 'Spectroscopy':
                self.log.debug(
                    'Overscan Region: {:s}'.format(str(self.overscan_region)))
                ccd = image_overscan(ccd=ccd,
                                     overscan_region=self.overscan_region)
            ccd = image_trim(ccd,
                             trim_section=self.trim_section,
                             trim_type='trimsec')
            master_bias_list.append(ccd)

        # combine bias for spectroscopy
        self.master_bias = ccdproc.combine(master_bias_list,
                                           method='median',
                                           sigma_clip=True,
                                           sigma_clip_low_thresh=3.0,
                                           sigma_clip_high_thresh=3.0,
                                           add_keyword=False)

        # add name of images used to create master bias
        for n in range(len(bias_file_list)):
            self.master_bias.header['GSP_IC{:02d}'.format(n + 1)] = (
                bias_file_list[n],
                'Image used to create master bias')

        write_fits(ccd=self.master_bias,
                   full_path=self.master_bias_name,
                   combined=True)

        self.log.info('Created master bias: ' + self.master_bias_name)
        return self.master_bias, self.master_bias_name

    def create_master_flats(self, flat_group, target_name=''):
        """Creates master flats

        Using a list of compatible flat images it combines them using median and
        1-sigma clipping. Also it apply all previous standard calibrations to
        each image.

        Args:
            flat_group (DataFrame): :class:`~pandas.DataFrame` instance.
              Contains a list of compatible flat images
            target_name (str): Science target name. This is used in some science
                case uses only.


        Returns:
            The master flat :class:`~astropy.nddata.CCDData` instance and the name of under which
            the master flat was stored.

        """

        flat_file_list = flat_group.file.tolist()
        cleaned_flat_list = []
        master_flat_list = []
        master_flat_name = None
        self.log.info('Creating Master Flat')
        for flat_file in flat_file_list:
            # print(f_file)
            image_full_path = os.path.join(self.args.raw_path, flat_file)
            ccd = read_fits(image_full_path, technique=self.technique)
            self.log.debug('Loading flat image: ' + image_full_path)
            if master_flat_name is None:

                master_flat_name = self.name_master_flats(
                    header=ccd.header,
                    group=flat_group,
                    target_name=target_name)

            if self.technique == 'Spectroscopy':
                ccd = image_overscan(ccd, overscan_region=self.overscan_region)

                ccd = image_trim(ccd=ccd,
                                 trim_section=self.trim_section,
                                 trim_type='trimsec')

                ccd = ccdproc.subtract_bias(ccd,
                                            self.master_bias,
                                            add_keyword=False)
                ccd.header['GSP_BIAS'] = (
                    os.path.basename(self.master_bias_name),
                    'Master bias image')

            elif self.technique == 'Imaging':
                ccd = image_trim(ccd=ccd,
                                 trim_section=self.trim_section,
                                 trim_type='trimsec')

                ccd = ccdproc.subtract_bias(ccd,
                                            self.master_bias,
                                            add_keyword=False)

                ccd.header['GSP_BIAS'] = (
                    os.path.basename(self.master_bias_name),
                    'Master bias image')
            else:
                self.log.error('Unknown observation technique: ' +
                               self.technique)
            if self._is_file_saturated(ccd=ccd):
                self.log.warning('Removing saturated image {:s}. '
                                 'Use --saturation to change saturation '
                                 'level'.format(flat_file))
                continue
            else:
                cleaned_flat_list.append(flat_file)
                master_flat_list.append(ccd)

        if master_flat_list != []:
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

            self.log.info('Created Master Flat: ' + master_flat_name)
            return master_flat, master_flat_name
            # print(master_flat_name)
        else:
            self.log.error('Empty flat list. Check that they do not exceed the '
                           'saturation limit.')
            return None, None

    def name_master_flats(self, header, group, target_name='', get=False):
        """Defines the name of a master flat or what master flat is compatible
        with a given data

        Given the header of a flat image this method will look for certain
        keywords that are unique to a given instrument configuration therefore
        they are used to discriminate compatibility.
        
        It can be used to define a master flat when creating it or find a base
        name to match existing master flat files thus finding a compatible one
        for a given non-flat image.

        Args:
            header (object): Fits header. Instance of
                astropy.io.fits.header.Header
            group (object): :class:`~pandas.DataFrame` instance. Contains filenames as
                well as other important keywords that are defined in
                ccd.night_organizer.NightOrganizer.keywords
            target_name (str): Optional science target name to be added to the
                master flat name.
            get (bool): This option is used when trying to find a suitable 
                master flat for a given data.

        Returns:
            A master flat name, or basename to find a match among existing
            files.

        """
        master_flat_name = os.path.join(self.args.red_path, 'master_flat')
        sunset = datetime.datetime.strptime(self.sun_set,
                                            "%Y-%m-%dT%H:%M:%S.%f")

        sunrise = datetime.datetime.strptime(self.sun_rise,
                                             "%Y-%m-%dT%H:%M:%S.%f")

        afternoon_twilight = datetime.datetime.strptime(self.evening_twilight,
                                                        "%Y-%m-%dT%H:%M:%S.%f")

        morning_twilight = datetime.datetime.strptime(self.morning_twilight,
                                                      "%Y-%m-%dT%H:%M:%S.%f")

        date_obs = datetime.datetime.strptime(header['DATE-OBS'],
                                              "%Y-%m-%dT%H:%M:%S.%f")
        # print(sunset, date_obs, evening_twilight)
        # print(' ')
        if target_name != '':
            target_name = '_' + target_name
        if not get:
            # TODO (simon): There must be a pythonic way to do this
            if (date_obs < sunset) or (date_obs > sunrise):
                dome_sky = '_dome'
            elif (sunset < date_obs < afternoon_twilight) or\
                    (morning_twilight < date_obs < sunrise):
                dome_sky = '_sky'
            elif afternoon_twilight < date_obs < morning_twilight:
                dome_sky = '_night'
            else:
                dome_sky = '_other'
        else:
            dome_sky = '*'

        if self.technique == 'Spectroscopy':
            if group.grating.unique()[0] != '<NO GRATING>':
                flat_grating = '_' + re.sub('[A-Za-z_-]',
                                            '',
                                            group.grating.unique()[0])

                # self.spec_mode is an instance of SpectroscopicMode
                wavmode = self.spec_mode(header=header)
            else:
                flat_grating = '_no_grating'
                wavmode = ''

            flat_slit = re.sub('[A-Za-z" ]',
                               '',
                               group.slit.unique()[0])

            filter2 = group['filter2'].unique()[0]
            if filter2 == '<NO FILTER>':
                filter2 = ''
            else:
                filter2 = '_' + filter2

            master_flat_name += target_name\
                + flat_grating\
                + wavmode\
                + filter2\
                + '_'\
                + flat_slit\
                + dome_sky\
                + '.fits'

        elif self.technique == 'Imaging':
            flat_filter = re.sub('-', '_', group['filter'].unique()[0])
            master_flat_name += '_' + flat_filter + dome_sky + '.fits'
        # print(master_flat_name)
        return master_flat_name

    def process_spectroscopy_science(self, science_group, save_all=False):
        """Process Spectroscopy science images.

        This function handles the full image reduction process for science
        files. if save_all is set to True, all intermediate steps are saved.

        Args:
            science_group (object): :class:`~pandas.DataFrame` instance that contains a
                list of science images that where observed at the same pointing
                and time. It also contains a set of selected keywords from the
                image's header.
            save_all (bool): If True the pipeline will save all the intermadiate
                files such as after overscan correction or bias corrected and
                etc.

        """
        # TODO (simon): The code here is too crowded.
        # TODO cont. Create other functions as necessary

        target_name = ''
        slit_trim = None
        master_bias = None
        master_flat = None
        master_flat_name = None
        obstype = science_group.obstype.unique()
        # print(obstype)
        if 'OBJECT' in obstype or 'COMP' in obstype:
            object_comp_group = science_group[
                (science_group.obstype == 'OBJECT') |
                (science_group.obstype == 'COMP')]

            if 'OBJECT' in obstype:
                target_name = science_group.object[science_group.obstype ==
                                                   'OBJECT'].unique()[0]

                self.log.info('Processing Science Target: '
                              '{:s}'.format(target_name))
            else:
                # TODO (simon): This does not make any sense
                self.log.info('Processing Comparison Lamp: '
                              '{:s}'.format(target_name))

            if 'FLAT' in obstype and not self.args.ignore_flats:
                flat_sub_group = science_group[science_group.obstype == 'FLAT']

                master_flat, master_flat_name = \
                    self.create_master_flats(flat_group=flat_sub_group,
                                             target_name=target_name)
            elif self.args.ignore_flats:
                self.log.warning('Ignoring creation of Master Flat by request.')
                master_flat = None
                master_flat_name = None
            else:
                self.log.info('Attempting to find a suitable Master Flat')
                object_list = object_comp_group.file.tolist()

                # grab a random image from the list
                random_image = random.choice(object_list)

                # define random image full path
                random_image_full = os.path.join(self.args.raw_path,
                                                 random_image)

                # read the random chosen file
                ccd = CCDData.read(random_image_full, unit=u.adu)

                if not self.args.ignore_flats:
                    # define the master flat name
                    master_flat_name = self.name_master_flats(
                        header=ccd.header,
                        group=object_comp_group,
                        get=True)

                    # load the best flat based on the name previously defined
                    master_flat, master_flat_name = \
                        get_best_flat(flat_name=master_flat_name,
                                      path=self.args.red_path)
                    if (master_flat is None) and (master_flat_name is None):
                        self.log.critical('Failed to obtain master flat')

            if master_flat is not None and not self.args.ignore_flats:
                self.log.debug('Attempting to find slit trim section')
                slit_trim = get_slit_trim_section(master_flat=master_flat)
            elif self.args.ignore_flats:
                self.log.warning('Slit Trimming will be skipped, '
                                 '--ignore-flats is activated')
            else:
                self.log.info("Master flat inexistent, can't find slit trim "
                              "section")
            if slit_trim is not None:

                master_flat = image_trim(ccd=master_flat,
                                         trim_section=slit_trim,
                                         trim_type='slit')

                if self.master_bias is not None:
                    master_bias = image_trim(ccd=self.master_bias,
                                             trim_section=slit_trim,
                                             trim_type='slit')
            else:
                try:
                    master_bias = self.master_bias.copy()
                except AttributeError:
                    master_bias = None

            norm_master_flat = None
            norm_master_flat_name = None
            all_object_image = []
            all_comp_image = []
            for science_image in object_comp_group.file.tolist():
                self.out_prefix = ''

                # define image full path
                image_full_path = os.path.join(self.args.raw_path,
                                               science_image)

                # load image
                ccd = read_fits(image_full_path, technique=self.technique)

                # apply overscan
                ccd = image_overscan(ccd, overscan_region=self.overscan_region)
                self.out_prefix += 'o_'

                if save_all:
                    full_path = os.path.join(self.args.red_path,
                                             self.out_prefix + science_image)

                    # ccd.write(full_path, clobber=True)
                    write_fits(ccd=ccd, full_path=full_path)

                if slit_trim is not None:
                    # There is a double trimming of the image, this is to match
                    # the size of the other data
                    # TODO (simon): Potential problem here
                    ccd = image_trim(ccd=ccd,
                                     trim_section=self.trim_section,
                                     trim_type='trimsec')

                    ccd = image_trim(ccd=ccd,
                                     trim_section=slit_trim,
                                     trim_type='slit')

                    self.out_prefix = 'st' + self.out_prefix

                    if save_all:
                        full_path = os.path.join(
                            self.args.red_path,
                            self.out_prefix + science_image)

                        # ccd.write(full_path, clobber=True)
                        write_fits(ccd=ccd, full_path=full_path)

                else:
                    ccd = image_trim(ccd=ccd,
                                     trim_section=self.trim_section,
                                     trim_type='trimsec')

                    self.out_prefix = 't' + self.out_prefix

                    if save_all:
                        full_path = os.path.join(
                            self.args.red_path,
                            self.out_prefix + science_image)

                        # ccd.write(full_path, clobber=True)
                        write_fits(ccd=ccd, full_path=full_path)

                if not self.args.ignore_bias:
                    # TODO (simon): Add check that bias is compatible

                    ccd = ccdproc.subtract_bias(ccd=ccd,
                                                master=master_bias,
                                                add_keyword=False)

                    self.out_prefix = 'z' + self.out_prefix

                    ccd.header['GSP_BIAS'] = (
                        os.path.basename(self.master_bias_name),
                        'Master bias image')

                    if save_all:
                        full_path = os.path.join(
                            self.args.red_path,
                            self.out_prefix + science_image)

                        # ccd.write(full_path, clobber=True)
                        write_fits(ccd=ccd, full_path=full_path)
                else:
                    self.log.warning('Ignoring bias correction by request.')

                # Do flat correction
                if master_flat is None or master_flat_name is None:
                    self.log.warning('The file {:s} will not be '
                                     'flatfielded'.format(science_image))
                elif self.args.ignore_flats:
                    self.log.warning('Ignoring flatfielding by request.')
                else:
                    if norm_master_flat is None:

                        norm_master_flat, norm_master_flat_name = \
                            normalize_master_flat(
                                master=master_flat,
                                name=master_flat_name,
                                method=self.args.flat_normalize,
                                order=self.args.norm_order)

                    ccd = ccdproc.flat_correct(ccd=ccd,
                                               flat=norm_master_flat,
                                               add_keyword=False)

                    self.out_prefix = 'f' + self.out_prefix

                    ccd.header['GSP_FLAT'] = (
                        os.path.basename(norm_master_flat_name),
                        'Master flat image')

                    # propagate master flat normalization method
                    ccd.header['GSP_NORM'] = norm_master_flat.header['GSP_NORM']

                    if save_all:
                        full_path = os.path.join(
                            self.args.red_path,
                            self.out_prefix + science_image)

                        # ccd.write(full_path, clobber=True)
                        write_fits(ccd=ccd, full_path=full_path)

                ccd, prefix = call_cosmic_rejection(
                    ccd=ccd,
                    image_name=science_image,
                    out_prefix=self.out_prefix,
                    red_path=self.args.red_path,
                    dcr_par=self.args.dcr_par_dir,
                    keep_files=self.args.keep_cosmic_files,
                    method=self.args.clean_cosmic,
                    save=True)

                self.out_prefix = prefix

                if ccd is not None:
                    if ccd.header['OBSTYPE'] == 'OBJECT':
                        self.log.debug("Appending OBJECT image for combination")
                        all_object_image.append(ccd)
                    elif ccd.header['OBSTYPE'] == 'COMP':
                        self.log.debug("Appending COMP image for combination")
                        all_comp_image.append(ccd)
                    else:
                        self.log.error("Unknown OBSTYPE = {:s}"
                                       "".format(ccd.header['OBSTYPE']))
                else:
                    self.log.warning("Cosmic ray rejection returned a None.")

            if self.args.combine:
                self.log.warning("Combination of data is experimental.")
                if len(all_object_image) > 1:
                    print(len(all_object_image))
                    self.log.info("Combining {:d} OBJECT images"
                                  "".format(len(all_object_image)))

                    object_group = object_comp_group[
                        object_comp_group.obstype == "OBJECT"]

                    print(object_group, len(all_object_image))

                    combine_data(all_object_image,
                                 dest_path=self.args.red_path,
                                 prefix=self.out_prefix,
                                 save=True)
                elif len(all_object_image) == 1:
                    # write_fits(all_object_image[0])
                    pass

                else:
                    self.log.error("No OBJECT images to combine")

                if len(all_comp_image) > 1:
                    self.log.info("Combining {:d} COMP images"
                                  "".format(len(all_comp_image)))
                    # comp_group = object_comp_group[
                    #     object_comp_group.obstype == "COMP"]
                    combine_data(all_comp_image,
                                 dest_path=self.args.red_path,
                                 prefix=self.out_prefix,
                                 save=True)
                else:
                    self.log.error("No COMP images to combine")
            else:
                self.log.debug("Combine is disabled (default)")

        elif 'FLAT' in obstype:
            self.queue.append(science_group)
            self.log.warning('Only flats found in this group')
            flat_sub_group = science_group[science_group.obstype == 'FLAT']
            # TODO (simon): Find out if these variables are useful or not
            self.create_master_flats(flat_group=flat_sub_group)
        else:
            self.log.error('There is no valid datatype in this group')

    def process_imaging_science(self, imaging_group):
        """Does image reduction for science imaging data.

        Args:
            imaging_group (object): :class:`~pandas.DataFrame` instance that contains a
              list of science data that are compatible with a given
              instrument configuration and can be reduced together.


        """
        # pick a random image in order to get a header
        random_image = random.choice(imaging_group.file.tolist())
        path_random_image = os.path.join(self.args.raw_path, random_image)
        sample_file = CCDData.read(path_random_image, unit=u.adu)

        master_flat_name = self.name_master_flats(header=sample_file.header,
                                                  group=imaging_group,
                                                  get=True)

        self.log.debug('Got {:s} for master flat name'.format(master_flat_name))

        master_flat, master_flat_name = get_best_flat(
            flat_name=master_flat_name,
            path=self.args.red_path)

        if master_flat is not None:
            for image_file in imaging_group.file.tolist():

                # start with an empty prefix
                self.out_prefix = ''

                image_full_path = os.path.join(self.args.raw_path, image_file)
                ccd = read_fits(image_full_path, technique=self.technique)

                # Trim image
                ccd = image_trim(ccd=ccd,
                                 trim_section=self.trim_section,
                                 trim_type='trimsec')
                self.out_prefix = 't_'
                if not self.args.ignore_bias:

                    ccd = ccdproc.subtract_bias(ccd,
                                                self.master_bias,
                                                add_keyword=False)

                    self.out_prefix = 'z' + self.out_prefix

                    ccd.header['GSP_BIAS'] = (
                        os.path.basename(self.master_bias_name),
                        'Master bias image')

                # apply flat correction
                ccd = ccdproc.flat_correct(ccd,
                                           master_flat,
                                           add_keyword=False)

                self.out_prefix = 'f' + self.out_prefix

                ccd.header['GSP_FLAT'] = (
                    os.path.basename(master_flat_name),
                    'Master Flat image')

                if self.args.clean_cosmic:
                    ccd = astroscrappy_lacosmic(
                        ccd=ccd,
                        red_path=self.args.red_path,
                        save_mask=self.args.keep_cosmic_files)

                    self.out_prefix = 'c' + self.out_prefix
                else:
                    print('Clean Cosmic ' + str(self.args.clean_cosmic))

                final_name = os.path.join(self.args.red_path,
                                          self.out_prefix + image_file)
                # ccd.write(final_name, clobber=True)
                write_fits(ccd=ccd, full_path=final_name)
                self.log.info('Created science file: {:s}'.format(final_name))
        else:
            self.log.error('Can not process data without a master flat')


if __name__ == '__main__':
    pass
