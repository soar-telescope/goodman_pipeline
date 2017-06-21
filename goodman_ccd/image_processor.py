from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import ccdproc
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import logging
import random
import re
import os
import glob
import datetime
import pandas
from .wavmode_translator import SpectroscopicMode
import sys
from ccdproc import CCDData
from scipy import interpolate
from scipy import signal

from astropy.convolution import (convolve,
                                 convolve_fft,
                                 Gaussian2DKernel,
                                 Tophat2DKernel)

from astropy.modeling import (models,
                              fitting)

from .core import (image_overscan,
                   image_trim,
                   get_slit_trim_section,
                   cosmicray_rejection,
                   get_best_flat,
                   normalize_master_flat,
                   dcr_cosmicray_rejection)

log = logging.getLogger('goodmanccd.imageprocessor')


class ImageProcessor(object):

    def __init__(self, args, data_container):
        """

        Args:
            args (object): Argparse object
            data_container (object): Contains relevant information of the night
            and the data itself.
        """
        self.args = args
        self.instrument = data_container.instrument
        self.technique = data_container.technique
        self.bias = data_container.bias
        self.day_flats = data_container.day_flats
        self.dome_flats = data_container.dome_flats
        self.sky_flats = data_container.sky_flats
        self.data_groups = data_container.data_groups
        self.sun_set = data_container.sun_set_time
        self.sun_rise = data_container.sun_rise_time
        self.morning_twilight = data_container.morning_twilight
        self.evening_twilight = data_container.evening_twilight
        self.pixel_scale = 0.15 * u.arcsec
        self.queue = None
        self.trim_section = self.define_trim_section(technique=self.technique)
        self.overscan_region = self.get_overscan_region()
        self.spec_mode = SpectroscopicMode()
        self.master_bias = None
        self.out_prefix = None

    def __call__(self, *args, **kwargs):
        """

        Args:
            *args:
            **kwargs:

        Returns:

        """

        for group in [self.bias,
                      self.day_flats,
                      self.dome_flats,
                      self.sky_flats,
                      self.data_groups]:
            if group is not None:
                for sub_group in group:
                    group_obstype = sub_group.obstype.unique()
                    if len(group_obstype) == 1 and group_obstype[0] == 'BIAS':
                        log.debug('Create Master Bias')
                        self.create_master_bias(sub_group)
                    elif len(group_obstype) == 1 and group_obstype[0] == 'FLAT':
                        log.debug('Create Master FLATS')
                        self.create_master_flats(sub_group)
                    else:
                        log.debug('Process Data Group')
                        if self.technique == 'Spectroscopy':
                            self.process_spectroscopy_science(sub_group)
                        else:
                            log.info('Processing Imaging Science Data')
                            self.process_imaging_science(sub_group)
        # print('data groups ', len(self.data_groups))
        if self.queue is not None:
            if len(self.queue) > 1:
                self.queue = pandas.concat(self.queue, axis=0).reset_index()
            else:
                self.queue = self.queue[0]
            print('QUEUE')
            print(self.queue)
            for sub_group in self.queue:
                print(sub_group)

    def define_trim_section(self, technique=None):
        """

        Returns:

        """

        assert technique is not None
        trim_section = None
        # TODO (simon): Consider binning and possibly ROIs for trim section
        log.warning('Determining trim section. Assuming you have only one kind '
                    'of data in this folder')
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
                ccd = CCDData.read(sample_image, unit=u.adu)

                # serial binning - dispersion binning
                # parallel binngin - spatial binning
                serial_binning, \
                parallel_binning = [int(x) for x
                                    in ccd.header['CCDSUM'].split()]

                # Trim section is valid for Blue and Red Camera Binning 1x1 and
                # Spectroscopic ROI
                if technique == 'Spectroscopy':

                    # left
                    l = int(np.ceil(51. / serial_binning))

                    # right
                    r = int(4110 / serial_binning)

                    # bottom
                    b = 1

                    #top
                    t = int(1896 / parallel_binning)
                    # TODO (simon): Need testing
                    trim_section = '[{:d}:{:d},{:d}:{:d}]'.format(l, r, b, t)

                elif technique == 'Imaging':
                    trim_section = ccd.header['TRIMSEC']

                log.info('Trim Section: %s', trim_section)
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

        Notes
        1 based
        IRAF V3 convention oppsoite to Python (poner en returns)

        Returns:
            overscan_region (str): Region for overscan in the format
            '[min:max,:]' where min is the starting point and max is the end
            point of the overscan region.

        """
        # TODO (simon): elaborate Notes on docstrings
        log.warning('Determining Overscan Region. Assuming you have only one '
                    'kind of binning in the data.')
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
                log.debug('Overscan Sample File ' + sample_image)
                ccd = CCDData.read(sample_image, unit=u.adu)

                # Image height - spatial direction
                h = ccd.data.shape[0]

                # Image width - spectral direction
                w = ccd.data.shape[1]

                # Take the binnings
                serial_binning, parallel_binning = \
                    [int(x) for x in ccd.header['CCDSUM'].split()]

                if self.technique == 'Spectroscopy':
                    log.info('Overscan regions has been tested for ROI '
                             'Spectroscopic 1x1, 2x2 and 3x3')

                    # define l r b and t to avoid local variable might be
                    # defined before assignment warning
                    l, r, b, t = [None] * 4
                    if self.instrument == 'Red':
                        # for red camera it is necessary to eliminate the first
                        # rows/columns (depends on the point of view) because
                        # they come with an abnormal high signal. Usually the
                        # first 5 pixels. In order to find the corresponding
                        # value for the subsequent binnings divide by the
                        # binning size.
                        # The numbers 6 and 49 where obtained from visual
                        # inspection

                        # left
                        l = int(np.ceil(6. / serial_binning))
                        # right
                        r = int(49. / serial_binning)
                        # bottom
                        b = 1
                        # top
                        t = h
                    elif self.instrument == 'Blue':
                        # 16 is the length of the overscan region with no
                        # binning.

                        # left
                        l = 1
                        # right
                        r = int(16. / serial_binning)
                        # bottom
                        b = 1
                        # top
                        t = h

                    overscan_region = '[{:d}:{:d},{:d}:{:d}]'.format(l, r, b, t)

                elif self.technique == 'Imaging':
                    log.warning("Imaging mode doesn't have overscan region. "
                                "Use bias instead.")
                    if self.bias is None:
                        log.warning('Bias are needed for Imaging mode')
                    overscan_region = None
                else:
                    overscan_region = None
                log.info('Overscan Region: %s', overscan_region)
                return overscan_region

    def create_master_bias(self, bias_group):
        """

        Args:
            bias_group:

        Returns:

        """
        bias_file_list = bias_group.file.tolist()
        default_bias_name =os.path.join(self.args.red_path, 'master_bias.fits')
        search_bias_name = re.sub('.fits', '*.fits', default_bias_name)
        n_bias = len(glob.glob(search_bias_name))
        if n_bias > 0:
            new_bias_name = re.sub('.fits',
                                   '_{:d}.fits'.format(n_bias + 1),
                                   default_bias_name)

            log.info('New name for master bias: ' + new_bias_name)
        else:
            new_bias_name = default_bias_name
        # TODO (simon): Review whether it is necessary to discriminate by
        # TODO technique
        if self.technique == 'Spectroscopy':
            master_bias_list = []
            log.info('Creating master bias')
            for image_file in bias_file_list:
                # print(image_file)
                image_full_path = os.path.join(self.args.raw_path, image_file)
                log.debug('Overscan Region: {:s}'.format(self.overscan_region))
                ccd = CCDData.read(image_full_path, unit=u.adu)
                log.debug('Loading bias image: ' + image_full_path)
                ccd = image_overscan(ccd, overscan_region=self.overscan_region)
                ccd = image_trim(ccd, trim_section=self.trim_section)
                master_bias_list.append(ccd)

            # combine bias for spectroscopy
            self.master_bias = ccdproc.combine(master_bias_list,
                                               method='median',
                                               sigma_clip=True,
                                               sigma_clip_low_thresh=3.0,
                                               sigma_clip_high_thresh=3.0,
                                               add_keyword=False)

            # write master bias to file
            self.master_bias.write(new_bias_name, clobber=True)
            log.info('Created master bias: ' + new_bias_name)

        elif self.technique == 'Imaging':
            master_bias_list = []
            log.info('Creating master bias')
            for image_file in bias_file_list:
                image_full_path = os.path.join(self.args.raw_path, image_file)
                ccd = CCDData.read(image_full_path, unit=u.adu)
                log.debug('Loading bias image: {:s}'.format(image_full_path))
                ccd = image_trim(ccd, trim_section=self.trim_section)
                master_bias_list.append(ccd)

            # combine bias for imaging
            self.master_bias = ccdproc.combine(master_bias_list,
                                               method='median',
                                               sigma_clip=True,
                                               sigma_clip_low_thresh=3.0,
                                               sigma_clip_high_thresh=3.0,
                                               add_keyword=False)

            # write master bias to file
            self.master_bias.write(new_bias_name, clobber=True)
            log.info('Created master bias: ' + new_bias_name)

    def create_master_flats(self, flat_group, target_name=''):
        """

        Args:
            flat_group:
            target_name:

        Returns:

        """

        flat_file_list = flat_group.file.tolist()
        master_flat_list = []
        master_flat_name = None
        log.info('Creating Master Flat')
        for flat_file in flat_file_list:
            # print(f_file)
            image_full_path = os.path.join(self.args.raw_path, flat_file)
            ccd = CCDData.read(image_full_path, unit=u.adu)
            log.debug('Loading flat image: ' + image_full_path)
            if master_flat_name is None:

                master_flat_name = self.name_master_flats(
                    header=ccd.header,
                    group=flat_group,
                    target_name=target_name)

            if self.technique == 'Spectroscopy':
                # plt.title('Before Overscan')
                # plt.imshow(ccd.data, clim=(-100, 0))
                # plt.show()
                ccd = image_overscan(ccd, overscan_region=self.overscan_region)
                # plt.title('After Overscan')
                # plt.imshow(ccd.data, clim=(-100, 0))
                # plt.show()
                ccd = image_trim(ccd, trim_section=self.trim_section)
                # plt.title('After Trimming')
                # plt.imshow(ccd.data, clim=(-100, 0))
                # plt.show()
            elif self.technique == 'Imaging':
                ccd = image_trim(ccd, trim_section=self.trim_section)
                ccd = ccdproc.subtract_bias(ccd,
                                            self.master_bias,
                                            add_keyword=False)
            else:
                log.error('Unknown observation technique: ' + self.technique)
            if ccd.data.max() > self.args.saturation_limit:
                log.warning('Removing saturated image {:s}. Use --saturation '
                            'to change saturation level'.format(flat_file))
                # plt.plot(ccd.data[802,:])
                # plt.show()
                # print(ccd.data.max())
                continue
            else:
                master_flat_list.append(ccd)
        if master_flat_list != []:
            master_flat = ccdproc.combine(master_flat_list,
                                          method='median',
                                          sigma_clip=True,
                                          sigma_clip_low_thresh=1.0,
                                          sigma_clip_high_thresh=1.0,
                                          add_keyword=False)
            master_flat.write(master_flat_name, clobber=True)
            # plt.imshow(master_flat.data, clim=(-100,0))
            # plt.show()
            log.info('Created Master Flat: ' + master_flat_name)
            return master_flat, master_flat_name
            # print(master_flat_name)
        else:
            log.error('Empty flat list. Probably they exceed the '
                      'saturation limit.')
            return None, None


    def name_master_flats(self, header, group, target_name='', get=False):
        """

        Args:
            header:
            group:
            target_name:
            get:

        Returns:

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

            master_flat_name += target_name +\
                                flat_grating +\
                                wavmode +\
                                filter2 +\
                                '_' +\
                                flat_slit +\
                                dome_sky +\
                                '.fits'

        elif self.technique == 'Imaging':
            flat_filter = re.sub('-', '_', group['filter'].unique()[0])
            master_flat_name += '_' + flat_filter + dome_sky + '.fits'
        # print(master_flat_name)
        return master_flat_name

    def process_spectroscopy_science(self, science_group):
        """

        Args:
            science_group:

        Returns:

        """
        # TODO (simon): The code here is too crowded.
        # TODO cont. Create other functions as necessary

        target_name = ''
        slit_trim = None
        obstype = science_group.obstype.unique()
        # print(obstype)
        if 'OBJECT' in obstype or 'COMP' in obstype:
            object_group = science_group[(science_group.obstype == 'OBJECT') |
                                         (science_group.obstype == 'COMP')]
            if 'OBJECT' in obstype:
                target_name = science_group.object[science_group.obstype ==
                                                   'OBJECT'].unique()[0]

                log.info('Processing Science Target: {:s}'.format(target_name))
            else:
                log.info('Processing Comparison Lamp: {:s}'.format(target_name))

            if 'FLAT' in obstype:
                flat_sub_group = science_group[science_group.obstype == 'FLAT']

                master_flat, master_flat_name =\
                    self.create_master_flats(flat_group=flat_sub_group,
                                             target_name=target_name)
            else:
                object_list = object_group.file.tolist()

                # grab a random image from the list
                random_image = random.choice(object_list)

                # define random image full path
                random_image_full = os.path.join(self.args.raw_path,
                                                 random_image)

                # read the random chosen file
                ccd = CCDData.read(random_image_full, unit=u.adu)

                # define the master flat name
                master_flat_name = self.name_master_flats(header=ccd.header,
                                                          group=object_group,
                                                          get=True)

                # load the best flat based on the name previously defined
                master_flat, master_flat_name =\
                    get_best_flat(flat_name=master_flat_name)
                if (master_flat is None) and (master_flat_name is None):
                    # attempt to find a set of flats in all the data

                    if self.queue is None:
                        self.queue = [science_group]
                    else:
                        self.queue.append(science_group)
                        # for sgroup in self.queue:
                        #     if isinstance(sgroup, pandas.DataFrame):
                        #         if sgroup.equals(science_group):
                        #             log.info('Science group already in queue')
                        #         else:
                        #             log.info('Appending science Group')
                        #             self.queue.append(science_group)
                        #     else:
                        #         log.error('Science Group is not an instance of
                        # pandas.DataFrame')
                    log.warning('Adding image group to Queue')
            if master_flat is not None:
                log.debug('Attempting to find slit trim section')
                slit_trim = get_slit_trim_section(master_flat=master_flat)
            else:
                log.info('Master flat inexistent, cant find slit trim section')
            if slit_trim is not None:

                master_flat = image_trim(ccd=master_flat,
                                         trim_section=slit_trim)

                master_bias = image_trim(ccd=self.master_bias,
                                         trim_section=slit_trim)
            else:
                master_bias = self.master_bias.copy()

            norm_master_flat = None
            for science_image in object_group.file.tolist():
                self.out_prefix = ''

                # define image full path
                image_full_path = os.path.join(self.args.raw_path,
                                               science_image)

                # load image
                ccd = CCDData.read(image_full_path, unit=u.adu)

                # apply overscan
                ccd = image_overscan(ccd, overscan_region=self.overscan_region)
                self.out_prefix += 'o_'

                if slit_trim is not None:
                    # There is a double trimming of the image, this is to match
                    # the size of the other data
                    ccd = image_trim(ccd=ccd, trim_section=self.trim_section)
                    ccd = image_trim(ccd=ccd, trim_section=slit_trim)
                    self.out_prefix = 'st' + self.out_prefix

                else:
                    ccd = image_trim(ccd=ccd, trim_section=self.trim_section)
                    self.out_prefix = 't' + self.out_prefix

                if not self.args.ignore_bias:
                    # TODO (simon): Add check that bias is compatible

                    ccd = ccdproc.subtract_bias(ccd=ccd,
                                                master=master_bias,
                                                add_keyword=False)

                    self.out_prefix = 'z' + self.out_prefix
                    ccd.header.add_history('Bias subtracted image')
                if master_flat is None or master_flat_name is None:
                    log.warning('The file {:s} will not be'
                                'flatfielded'.format(science_image))
                else:
                    if norm_master_flat is None:

                        norm_master_flat = normalize_master_flat(
                            master=master_flat,
                            name=master_flat_name,
                            method=self.args.flat_normalize,
                            order=self.args.norm_order)

                    ccd = ccdproc.flat_correct(ccd=ccd,
                                               flat=norm_master_flat,
                                               add_keyword=False)

                    self.out_prefix = 'f' + self.out_prefix

                    ccd.header.add_history('master flat norm_'
                                           '{:s}'.format(master_flat_name))

                full_path = os.path.join(self.args.red_path,
                                         self.out_prefix + science_image)
                ccd.write(full_path, clobber=True)
                log.info('Created science image: {:s}'.format(full_path))

                if self.args.clean_cosmic:

                    in_file = self.out_prefix + science_image

                    dcr_cosmicray_rejection(data_path=self.args.red_path,
                                            in_file=in_file,
                                            prefix='c',
                                            dcr_par_dir=self.args.dcr_par_dir,
                                            delete=self.args.keep_cosmic_files)

                    # ccd = cosmicray_rejection(ccd=ccd)
                    # self.out_prefix = 'c' + self.out_prefix
                    # plt.title('Cosmic')
                    # plt.imshow(ccd.data, clim=(0, 1000))
                    # plt.show()
                else:
                    log.warning('Clean Cosmic ' + str(self.args.clean_cosmic))


                

                # print(science_group)
        elif 'FLAT' in obstype:
            self.queue.append(science_group)
            log.warning('Only flats found in this group')
            flat_sub_group = science_group[science_group.obstype == 'FLAT']
            # TODO (simon): Find out if these variables are useful or not
            master_flat, master_flat_name = \
                self.create_master_flats(flat_group=flat_sub_group)
        else:
            log.error('There is no valid datatype in this group')

    def process_imaging_science(self, imaging_group):
        """

        Args:
            imaging_group:

        Returns:

        """
        # pick a random image in order to get a header
        random_image = random.choice(imaging_group.file.tolist())
        path_random_image = os.path.join(self.args.raw_path, random_image)
        sample_file = CCDData.read(path_random_image, unit=u.adu)

        master_flat_name = self.name_master_flats(header=sample_file.header,
                                                  group=imaging_group,
                                                  get=True)

        log.debug('Got {:s} for master flat name'.format(master_flat_name))

        master_flat, master_flat_name = get_best_flat(
            flat_name=master_flat_name)

        if master_flat is not None:
            for image_file in imaging_group.file.tolist():

                # start with an empty prefix
                self.out_prefix = ''

                image_full_path = os.path.join(self.args.raw_path, image_file)
                ccd = CCDData.read(image_full_path, unit=u.adu)

                # Trim image
                ccd = image_trim(ccd, trim_section=self.trim_section)
                self.out_prefix = 't_'
                if not self.args.ignore_bias:

                    ccd = ccdproc.subtract_bias(ccd,
                                                self.master_bias,
                                                add_keyword=False)

                    self.out_prefix = 'z' + self.out_prefix
                    ccd.header.add_history('Bias subtracted image')

                # apply flat correction
                ccd = ccdproc.flat_correct(ccd,
                                           master_flat,
                                           add_keyword=False)

                self.out_prefix = 'f' + self.out_prefix

                ccd.header.add_history(
                    'Flat corrected '
                    '{:s}'.format(master_flat_name.split('/')[-1]))

                if self.args.clean_cosmic:
                    ccd = cosmicray_rejection(ccd=ccd)
                    self.out_prefix = 'c' + self.out_prefix
                else:
                    print('Clean Cosmic ' + str(self.args.clean_cosmic))

                final_name = os.path.join(self.args.red_path,
                                          self.out_prefix + image_file)
                ccd.write(final_name, clobber=True)
                log.info('Created science file: {:s}'.format(final_name))
        else:
            log.error('Can not process data without a master flat')


if __name__ == '__main__':
    pass

