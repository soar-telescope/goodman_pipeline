from __future__ import print_function
from astropy import units as u
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel, Tophat2DKernel
from astropy.modeling import models, fitting
import ccdproc
from ccdproc import CCDData
from scipy import interpolate
from scipy import signal
import numpy as np
import logging
import random
import re
import os
import glob
import datetime
import pandas
from wavmode_translator import SpectroscopicMode
import sys
import matplotlib.pyplot as plt

log = logging.getLogger('goodmanccd.imageprocessor')


class ImageProcessor(object):

    def __init__(self, args, data_container):
        """

        Args:
            args (object): Argparse object
            data_container (object): Contains relevant information of the night and the data itself.
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
        self.trim_section = self.define_trim_section()
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

        for group in [self.bias, self.day_flats, self.dome_flats, self.sky_flats, self.data_groups]:
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



    def define_trim_section(self):
        """

        Returns:

        """

        # TODO (simon): Consider binning and possibly ROIs for trim section
        log.warning('Determining trim section. Assuming you have only one kind of data in this folder')
        for group in [self.bias, self.day_flats, self.dome_flats, self.sky_flats, self.data_groups]:
            if group is not None:
                # print(self.bias[0])
                image_list = group[0]['file'].tolist()
                sample_image = os.path.join(self.args.raw_path, random.choice(image_list))
                # header = fits.getheader(sample_image)
                # trim_section = header['TRIMSEC']
                # Trim section is valid for Blue and Red Camera Binning 1x1 Spectroscopic ROI
                trim_section = '[51:4110,1:1896]'
                log.info('Trim Section: %s', trim_section)
                return trim_section

    def get_overscan_region(self):
        """

        Returns:

        """

        log.warning('Determining Overscan Region. Assuming you have only one kind of data in this folder')
        for group in [self.bias, self.day_flats, self.dome_flats, self.sky_flats, self.data_groups]:
            if group is not None:
                # print(self.bias[0])
                image_list = group[0]['file'].tolist()
                sample_image = os.path.join(self.args.raw_path, random.choice(image_list))
                log.debug('Overscan Sample File ' + sample_image)
                ccd = CCDData.read(sample_image, unit=u.adu)
                serial_binning, paralell_bining = [int(x) for x in ccd.header['CCDSUM'].split()]
                log.info('Using predefined overscan region for binning 1x1')
                if self.technique == 'Spectroscopy':
                    if self.instrument == 'Red':
                        overscan_region = '[{:d}:{:d},1:1896]'.format(int(np.ceil(6. / serial_binning)), int(49./ serial_binning))
                    elif self.instrument == 'Blue':
                        overscan_region = '[1:{:d},1:1896]'.format(int(16. / serial_binning))
                elif self.technique == 'Imaging':
                    log.warning("Imaging mode doesn't have overscan region. Use bias instead.")
                    if self.bias is None:
                        log.warning('Bias are needed for Imaging mode')
                    overscan_region = None
                else:
                    overscan_region = None
                log.info('Overscan Region: %s', overscan_region)
                return overscan_region

    @staticmethod
    def image_overscan(ccd, overscan_region):
        """

        Args:
            ccd:

        Returns:

        """

        ccd = ccdproc.subtract_overscan(ccd, median=True, fits_section=overscan_region, add_keyword=False)
        ccd.header.add_history('Applied overscan correction ' + overscan_region)
        return ccd

    @staticmethod
    def image_trim(ccd, trim_section):
        """

        Args:
            ccd (object): CCDData Object
            trim_section (str):

        Returns:

        """

        ccd = ccdproc.trim_image(ccd, fits_section=trim_section, add_keyword=False)
        ccd.header.add_history('Trimmed image to ' + trim_section)
        return ccd

    def create_master_bias(self, bias_group):
        """

        Args:
            bias_group:

        Returns:

        """
        bias_file_list = bias_group.file.tolist()
        default_bias_name =os.path.join(self.args.red_path, 'master_bias.fits')
        if len(glob.glob(os.path.join(self.args.red_path, 'master_bias*fits'))) > 0:
            new_bias_name = re.sub(".fits",
                                   '_{:d}.fits'.format(len(glob.glob(os.path.join(self.args.red_path,
                                                                                  'master_bias*fits'))) + 1),
                                   default_bias_name)
            log.info('New name for master bias: ' + new_bias_name)
        else:
            new_bias_name = default_bias_name
        # TODO (simon): Review whether it is necessary to discriminate by technique
        if self.technique == 'Spectroscopy':
            master_bias_list = []
            for image_file in bias_file_list:
                # print(image_file)
                log.debug(self.overscan_region)
                ccd = CCDData.read(os.path.join(self.args.raw_path, image_file), unit=u.adu)
                log.info('Loading bias image: ' + os.path.join(self.args.raw_path, image_file))
                o_ccd = self.image_overscan(ccd, overscan_region=self.overscan_region)
                to_ccd = self.image_trim(o_ccd, trim_section=self.trim_section)
                master_bias_list.append(to_ccd)
            self.master_bias = ccdproc.combine(master_bias_list, method='median', sigma_clip=True,
                                               sigma_clip_low_thresh=3.0, sigma_clip_high_thresh=3.0, add_keyword=False)
            self.master_bias.write(new_bias_name, clobber=True)
            log.info('Created master bias: ' + new_bias_name)
        elif self.technique == 'Imaging':
            master_bias_list = []
            for image_file in bias_file_list:
                ccd = CCDData.read(os.path.join(self.args.raw_path, image_file), unit=u.adu)
                log.info('Loading bias image: ' + os.path.join(self.args.raw_path, image_file))
                t_ccd = self.image_trim(ccd, trim_section=self.trim_section)
                master_bias_list.append(t_ccd)
            self.master_bias = ccdproc.combine(master_bias_list, method='median', sigma_clip=True,
                                               sigma_clip_low_thresh=3.0, sigma_clip_high_thresh=3.0, add_keyword=False)
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
        for flat_file in flat_file_list:
            # print(f_file)
            ccd = CCDData.read(os.path.join(self.args.raw_path, flat_file), unit=u.adu)
            log.info('Loading flat image: ' + os.path.join(self.args.raw_path, flat_file))
            if master_flat_name is None:
                master_flat_name = self.name_master_flats(header=ccd.header, group=flat_group, target_name=target_name)
            if self.technique == 'Spectroscopy':
                # plt.title('Before Overscan')
                # plt.imshow(ccd.data, clim=(-100, 0))
                # plt.show()
                ccd = self.image_overscan(ccd, overscan_region=self.overscan_region)
                # plt.title('After Overscan')
                # plt.imshow(ccd.data, clim=(-100, 0))
                # plt.show()
                ccd = self.image_trim(ccd, trim_section=self.trim_section)
                # plt.title('After Trimming')
                # plt.imshow(ccd.data, clim=(-100, 0))
                # plt.show()
            elif self.technique == 'Imaging':
                ccd = self.image_trim(ccd, trim_section=self.trim_section)
                ccd = ccdproc.subtract_bias(ccd, self.master_bias, add_keyword=False)
            else:
                log.error('Unknown observation technique: ' + self.technique)
            if ccd.data.max() > self.args.saturation_limit:
                log.info('Removing saturated image {:s}. Use --saturation to change saturation level'.format(flat_file))
                # plt.plot(ccd.data[802,:])
                # plt.show()
                print(ccd.data.max())
                continue
            else:
                master_flat_list.append(ccd)
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
        sunset = datetime.datetime.strptime(self.sun_set, "%Y-%m-%dT%H:%M:%S.%f")
        sunrise = datetime.datetime.strptime(self.sun_rise, "%Y-%m-%dT%H:%M:%S.%f")
        afternoon_twilight = datetime.datetime.strptime(self.evening_twilight, "%Y-%m-%dT%H:%M:%S.%f")
        morning_twilight = datetime.datetime.strptime(self.morning_twilight, "%Y-%m-%dT%H:%M:%S.%f")
        date_obs = datetime.datetime.strptime(header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f")
        # print(sunset, date_obs, evening_twilight)
        # print(' ')
        if target_name != '':
            target_name = '_' + target_name
        if not get:
            if (date_obs < sunset) or (date_obs > sunrise):
                dome_sky = '_dome'
            elif (sunset < date_obs < afternoon_twilight) or (morning_twilight < date_obs < sunrise):
                dome_sky = '_sky'
            elif afternoon_twilight < date_obs < morning_twilight:
                dome_sky = '_night'
            else:
                dome_sky = '_other'
        else:
            dome_sky = '*'

        if self.technique == 'Spectroscopy':
            if group.grating.unique()[0] != '<NO GRATING>':
                flat_grating = '_' + re.sub('[A-Za-z_-]', '', group.grating.unique()[0])
                wavmode = self.spec_mode(header=header)
            else:
                flat_grating = '_no_grating'
                wavmode = ''
            flat_slit = re.sub('[A-Za-z" ]', '', group.slit.unique()[0])
            filter2 = group['filter2'].unique()[0]
            if filter2 == '<NO FILTER>':
                filter2 = ''
            else:
                filter2 = '_' + filter2
            master_flat_name += target_name + flat_grating + wavmode + filter2 + '_' + flat_slit + dome_sky + '.fits'
        elif self.technique == 'Imaging':
            flat_filter = re.sub('-', '_', group['filter'].unique()[0])
            master_flat_name += '_' + flat_filter + dome_sky + '.fits'
        # print(master_flat_name)
        return master_flat_name

    @staticmethod
    def get_slit_trim_section(master_flat):
        """Find the slit edges to trim all data

        Args:
            master_flat (object): CCDData instance

        Returns:
            trim_section (str): Trim section in spatial direction in the format [:,slit_lower_limit:slit_higher_limit]

        """


        x, y = master_flat.data.shape
        middle = int(y / 2.)
        ccd_section = master_flat.data[:, middle:middle + 200]
        ccd_section_median = np.median(ccd_section, axis=1)
        # spatial_axis = range(len(ccd_section_median))

        spatial_half = len(ccd_section_median) / 2.

        pseudo_derivative = np.array(
            [abs(ccd_section_median[i + 1] - ccd_section_median[i]) for i in range(0, len(ccd_section_median) - 1)])
        filtered_data = np.where(np.abs(pseudo_derivative > 0.5 * pseudo_derivative.max()),
                                 pseudo_derivative,
                                 None)
        peaks = signal.argrelmax(filtered_data, axis=0, order=3)[0]
        # print(peaks)

        trim_section = None
        if len(peaks) > 2 or peaks == []:
            log.debug('No trim section')
        else:
            # print(peaks, flat_files.grating[flat_files.file == file], flat_files.slit[flat_files.file == file])
            if len(peaks) == 2:
                low, high = peaks
                trim_section = '[:,{:d}:{:d}]'.format(low, high)
            elif len(peaks) == 1:
                if peaks[0] <= spatial_half:
                    trim_section = '[:,{:d}:{:d}]'.format(peaks[0], len(ccd_section_median))
                else:
                    trim_section = '[:,{:d}:{:d}]'.format(0, peaks[0])
        return trim_section

    def process_spectroscopy_science(self, science_group):
        """

        Args:
            science_group:

        Returns:

        """

        target_name = ''
        slit_trim = None
        obstype = science_group.obstype.unique()
        print(obstype)
        if 'OBJECT' in obstype or 'COMP' in obstype:
            object_group = science_group[(science_group.obstype == 'OBJECT') | (science_group.obstype == 'COMP')]
            if 'OBJECT' in obstype:
                target_name = science_group.object[science_group.obstype == 'OBJECT'].unique()[0]
                log.info('Processing Science Target: ' + target_name)
            else:
                # target_name = science_group.object[science_group.obstype == 'COMP'].unique()[0]
                log.info('Processing Comparison Lamp: ' + target_name)

            if 'FLAT' in obstype:
                flat_sub_group = science_group[science_group.obstype == 'FLAT']
                master_flat, master_flat_name = self.create_master_flats(flat_group=flat_sub_group, target_name=target_name)
            else:
                ccd = CCDData.read(os.path.join(self.args.raw_path, random.choice(object_group.file.tolist())), unit=u.adu)
                master_flat_name = self.name_master_flats(header=ccd.header, group=object_group, get=True)
                master_flat, master_flat_name = self.get_best_flat(flat_name=master_flat_name)
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
                        #         log.error('Science Group is not an instance of pandas.DataFrame')
                    log.warning('Adding image group to Queue')
            if master_flat is not None:
                log.debug('Attempting to find slit trim section')
                slit_trim = self.get_slit_trim_section(master_flat=master_flat)
            else:
                log.info('Master flat inexistent, cant find slit trim section')
            if slit_trim is not None:
                    master_flat = self.image_trim(ccd=master_flat, trim_section=slit_trim)
                    master_bias = self.image_trim(ccd=self.master_bias, trim_section=slit_trim)
            else:
                master_bias = self.master_bias.copy()

            norm_master_flat = None
            for science_image in object_group.file.tolist():
                self.out_prefix = ''
                ccd = CCDData.read(os.path.join(self.args.raw_path, science_image), unit=u.adu)
                # plt.title('Raw')
                # plt.imshow(ccd.data, clim=(0, 1000))
                # plt.show()
                ccd = self.image_overscan(ccd, overscan_region=self.overscan_region)
                # plt.title('Overscanned')
                # plt.imshow(ccd.data, clim=(0, 1000))
                # plt.show()
                self.out_prefix += 'o_'
                if slit_trim is not None:
                    # There is a double trimming of the image, this is to match the size of the other data
                    ccd = self.image_trim(ccd=ccd, trim_section=self.trim_section)
                    ccd = self.image_trim(ccd=ccd, trim_section=slit_trim)
                    self.out_prefix = 'st' + self.out_prefix

                else:
                    ccd = self.image_trim(ccd=ccd, trim_section=self.trim_section)
                    self.out_prefix = 't' + self.out_prefix
                # plt.title('Trimmed')
                # plt.imshow(ccd.data)
                # plt.show()
                if not self.args.ignore_bias:
                    ccd = ccdproc.subtract_bias(ccd=ccd, master=master_bias, add_keyword=False)
                    self.out_prefix = 'z' + self.out_prefix
                    ccd.header.add_history('Bias subtracted image')
                if master_flat is None or master_flat_name is None:
                    log.warning('The file ' + science_image + ' will not be flatfielded')
                else:
                    if self.args.flat_normalize == 'mean':
                        ccd = ccdproc.flat_correct(ccd=ccd, flat=master_flat, add_keyword=False)
                        self.out_prefix = 'f' + self.out_prefix
                        ccd.header.add_history('Flat corrected ' + master_flat_name.split('/')[-1])
                        # plt.title('Flat Mean Normalized')
                        # plt.imshow(ccd.data)
                        # plt.show()
                    elif norm_master_flat is None:
                        norm_master_flat = master_flat.copy()
                        if self.args.flat_normalize == 'simple':
                            # log.warning('This part of the code was left here for experimental purposes only')
                            log.info('Normalizing flat by model')
                            model_init = models.Chebyshev1D(degree=self.args.norm_order)
                            model_fitter = fitting.LevMarLSQFitter()
                            x_size, y_size = master_flat.data.shape
                            x_axis = range(y_size)

                            profile = np.median(master_flat.data, axis=0)
                            fit = model_fitter(model_init, x_axis, profile)
                            # plt.ion()
                            for i in range(x_size):
                                # fit = model_fitter(model_init, x_axis, master_flat.data[i])
                                norm_master_flat.data[i] = master_flat.data[i] / fit(x_axis)
                            #     plt.clf()
                            #     plt.xlim(0,1500)
                            #     plt.ylim(-1000, 2000)
                            #     plt.plot(norm_master_flat.data[i], color='r')
                            #     plt.plot(fit(x_axis), color='k')
                            #     plt.plot(master_flat.data[i], color='b')
                            #     plt.draw()
                            #     plt.pause(1)
                            # plt.ioff()
                            norm_master_flat.write(os.path.join('/'.join(master_flat_name.split('/')[0:-1]), 'norm_' + master_flat_name.split('/')[-1]), clobber=True)
                            # plt.title('simple Normalized Flat ')
                            # plt.imshow(norm_master_flat.data, clim=(-100,0))
                            # plt.show()

                            ccd.data = ccd.data / norm_master_flat.data
                            self.out_prefix = 'f' + self.out_prefix
                            ccd.header.add_history('Flat normalized simple model ' + master_flat_name.split('/')[-1])
                            # plt.title('Flat simple Normalized')
                            # plt.imshow(ccd.data, clim=(-100,900))
                            # plt.show()
                        elif self.args.flat_normalize == 'full':
                            log.warning('This part of the code was left here for experimental purposes only')
                            log.info('This procedure takes a lot to process, you might want to see other method')
                            model_init = models.Chebyshev1D(degree=self.args.norm_order)
                            model_fitter = fitting.LinearLSQFitter()
                            # fit = model_fitter(model_init, x_axis, profile)
                            x_size, y_size = master_flat.data.shape
                            x_axis = range(y_size)

                            for i in range(x_size):
                                fit = model_fitter(model_init, x_axis, master_flat.data[i])
                                norm_master_flat.data[i] = master_flat.data[i] / fit(x_axis)
                            log.debug(os.path.join('/'.join(master_flat_name.split('/')[0:-1]), 'norm_' + master_flat_name.split('/')[-1]))
                            norm_master_flat.write(os.path.join('/'.join(master_flat_name.split('/')[0:-1]), 'norm_' + master_flat_name.split('/')[-1]), clobber=True)
                            ccd.data = ccd.data / norm_master_flat.data
                            self.out_prefix = 'f' + self.out_prefix
                            ccd.header.add_history('Flat normalized line by line ' + master_flat_name.split('/')[-1])
                    else:
                        # mf_min = np.min(norm_master_flat.data)
                        # mf_max = np.max(norm_master_flat.data)
                        # plt.imshow(norm_master_flat.data, clim=(1- 0.3 * mf_min, 1 + 0.3 * mf_max))
                        # plt.show()
                        ccd.data = ccd.data / norm_master_flat.data
                        self.out_prefix = 'f' + self.out_prefix
                        ccd.header.add_history('Flat normalized {:s} '.format(self.args.flat_normalize.split('-')) + master_flat_name.split('/')[-1])
                if self.args.clean_cosmic:
                    ccd = self.cosmicray_rejection(ccd=ccd)
                    self.out_prefix = 'c' + self.out_prefix
                    # plt.title('Cosmic')
                    # plt.imshow(ccd.data, clim=(0, 1000))
                    # plt.show()
                else:
                    print('Clean Cosmic ' + str(self.args.clean_cosmic))
                # ccd.data = np.nan_to_num(ccd.data)
                # print('Length CCD Data', ccd.data.shape)
                # if self.args.interp_invalid:
                #     ccd = self.convolve(ccd=ccd)
                # else:
                #     log.warning('Replacing invalid values by numbers. Use --interpolate-invalid'
                #                 ' for a more advanced method.')
                #     ccd.data = np.nan_to_num(ccd.data)
                #     print(ccd.data.shape)
                # if self.args.debug_mode:
                #     x, y = ccd.data.shape
                #     print('x ', x, ' y ', y)
                #     for i in range(x):
                #         for e in range(y):
                #             if np.isnan(ccd.data[i, e]):
                #                 # print(i, x, e, y)
                #                 plt.plot(e, i, marker='o', color='r', mec='r')
                #             elif np.isneginf(ccd.data[i, e]):
                #                 plt.plot(e, i, marker='o', color='r', mec='r')
                #                 # print(" ")
                #
                #     plt.imshow(ccd.data, clim=(5, 150), cmap='cubehelix', origin='lower', interpolation='nearest')
                #     plt.show()
                ccd.write(os.path.join(self.args.red_path, self.out_prefix + science_image))
                log.info('Created science image: ' + os.path.join(self.args.red_path, self.out_prefix + science_image))
                

                # print(science_group)
        elif 'FLAT' in obstype:
            self.queue.append(science_group)
            log.warning('Only flats found in this group')
            flat_sub_group = science_group[science_group.obstype == 'FLAT']
            master_flat, master_flat_name = self.create_master_flats(flat_group=flat_sub_group)
        else:
            log.error('There is no valid datatype in this group')

    def process_imaging_science(self, imaging_group):
        """

        Args:
            imaging_group:

        Returns:

        """

        sample_file = CCDData.read(os.path.join(self.args.raw_path, random.choice(imaging_group.file.tolist())), unit=u.adu)
        master_flat_name = self.name_master_flats(sample_file.header, imaging_group, get=True)
        print(master_flat_name)
        master_flat, master_flat_name = self.get_best_flat(flat_name=master_flat_name)
        if master_flat is not None:
            for image_file in imaging_group.file.tolist():
                self.out_prefix = ''
                ccd = CCDData.read(os.path.join(self.args.raw_path, image_file), unit=u.adu)
                ccd = self.image_trim(ccd, trim_section=self.trim_section)
                self.out_prefix = 't_'
                if not self.args.ignore_bias:
                    ccd = ccdproc.subtract_bias(ccd, self.master_bias, add_keyword=False)
                    self.out_prefix = 'z' + self.out_prefix
                    ccd.header.add_history('Bias subtracted image')
                ccd = ccdproc.flat_correct(ccd, master_flat, add_keyword=False)
                self.out_prefix = 'f' + self.out_prefix
                ccd.header.add_history('Flat corrected ' + master_flat_name.split('/')[-1])
                if self.args.clean_cosmic:
                    ccd = self.cosmicray_rejection(ccd=ccd)
                    self.out_prefix = 'c' + self.out_prefix
                else:
                    print('Clean Cosmic ' + str(self.args.clean_cosmic))
                ccd.write(self.args.red_path, self.out_prefix + image_file, clobber=True)
                log.info('Created science file: ' + os.path.join(self.args.red_path, self.out_prefix + image_file))
        else:
            log.error('Can not process data without a master flat')

    def cosmicray_rejection(self, ccd):
        """Do cosmic ray rejection

        Notes:
            OBS: cosmic ray rejection is working pretty well by defining gain = 1. It's not working when we use the real
            gain of the image. In this case the sky level changes by a factor equal the gain.
            Function to determine the sigfrac and objlim: y = 0.16 * exptime + 1.2

        Args:
            ccd (object): CCDData Object

        Returns:
            ccd (object): The modified CCDData object

        """
        # TODO (simon): Validate this method
        if ccd.header['OBSTYPE'] == 'OBJECT':
            value = 0.16 * float(ccd.header['EXPTIME']) + 1.2
            log.info('Cleaning cosmic rays... ')
            # ccd.data= ccdproc.cosmicray_lacosmic(ccd, sigclip=2.5, sigfrac=value, objlim=value,
            ccd.data, _mask = ccdproc.cosmicray_lacosmic(ccd.data, sigclip=2.5, sigfrac=value, objlim=value,

                                                         gain=float(ccd.header['GAIN']),
                                                         readnoise=float(ccd.header['RDNOISE']),
                                                         satlevel=np.inf, sepmed=True, fsmode='median',
                                                         psfmodel='gaussy', verbose=False)
            # ccd.data = np.array(ccd.data, dtype=np.double) / float(ccd.header['GAIN'])

            # difference = np.nan_to_num(ccd - nccd)
            # print(_mask)
            # plt.imshow(_mask, interpolation='none')
            # plt.show()

            # plt.imshow(difference, clim=(0,1))
            # plt.show()

            ccd.header.add_history("Cosmic rays rejected with LACosmic")
            log.info("Cosmic rays rejected with LACosmic")
        else:
            log.info('Skipping cosmic ray rejection for image of datatype: {:s}'.format(ccd.header['OBSTYPE']))
        return ccd

    def remove_nan(self, ccd):
        # do some kind of interpolation to removen NaN and -INF
        # np.nan_to_num()
        # TODO (simon): Re-write in a more comprehensive way
        log.info('Removing NaN and INF by cubic interpolation')
        x = np.arange(0, ccd.data.shape[1])
        y = np.arange(0, ccd.data.shape[0])
        # mask invalid values
        array = np.ma.masked_invalid(ccd.data)
        xx, yy = np.meshgrid(x, y)
        # get only the valid values
        x1 = xx[~array.mask]
        y1 = yy[~array.mask]
        newarr = array[~array.mask]

        ccd.data = interpolate.griddata((x1, y1), newarr.ravel(),
                                   (xx, yy),
                                   method='cubic')

        # print('Length CCD Data', ccd.data.shape)

        plt.imshow(ccd.data, clim=(5, 150), cmap='cubehelix')
        plt.show()

        return ccd

    def convolve(self, ccd):
        binning = 1
        if self.instrument == 'Red':
            binning = int(ccd.header['PG5_4'])
        elif self.instrument == 'Blue':
            binning = int(ccd.header['PARAM22'])
        else:
            log.warning('No proper camera detected.')
        seeing = float(ccd.header['SEEING']) * u.arcsec
        print('seeing '+ str(seeing))
        print('Pixel Scale '+ str(self.pixel_scale))
        print('Binning '+ str(binning))
        fwhm = seeing / (self.pixel_scale * binning)
        # to deal with the cases when there is no seeing info
        fwhm = abs(fwhm/2.)

        gaussian_kernel = Gaussian2DKernel(stddev=1)
        tophat_kernel = Tophat2DKernel(5)
        new_ccd_data = convolve_fft(ccd.data, tophat_kernel, interpolate_nan=True, allow_huge=True)

        fig = plt.figure()
        a = fig.add_subplot(2, 1, 1)
        plt.imshow(ccd.data, clim=(5, 150), cmap='cubehelix', origin='lower', interpolation='nearest')
        a = fig.add_subplot(2,1,2)
        plt.imshow(new_ccd_data, clim=(5, 150), cmap='cubehelix', origin='lower', interpolation='nearest')
        plt.show()

        print(fwhm)
        return ccd

    @staticmethod
    def get_best_flat(flat_name):
        """

        Args:
            flat_name:

        Returns:

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

if __name__ == '__main__':
    pass

