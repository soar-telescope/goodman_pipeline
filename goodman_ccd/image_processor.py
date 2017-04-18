from __future__ import print_function
from astropy.io import fits
from astropy import units as u
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel, Tophat2DKernel
import ccdproc
from ccdproc import CCDData
from scipy import interpolate
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
                header = fits.getheader(sample_image)
                # trim_section = header['TRIMSEC']
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
                header = fits.getheader(sample_image)
                if header['CCDSUM'] == '1 1':
                    log.info('Using predefined overscan region for binning 1x1')
                    if self.technique == 'Spectroscopy':
                        overscan_region = '[5:45,1:4096]'
                    elif self.technique == 'Imaging':
                        log.warning("Imaging mode doesn't have overscan region. Use bias instead.")
                        if self.bias is None:
                            log.error('Bias are needed for Imaging mode')
                        overscan_region = None
                    else:
                        overscan_region = None
                else:
                    log.warning('There is no predefined overscan region')
                    overscan_region = None

                log.info('Overscan Region: %s', overscan_region)
                return overscan_region

    def image_overscan(self, ccd):
        """

        Args:
            ccd:

        Returns:

        """

        ccd = ccdproc.subtract_overscan(ccd, median=True, fits_section=self.overscan_region, add_keyword=False)
        ccd.header.add_history('Applied overscan correction ' + self.overscan_region)
        return ccd

    def image_trim(self, ccd):
        """

        Args:
            ccd:

        Returns:

        """

        ccd = ccdproc.trim_image(ccd, fits_section=self.trim_section, add_keyword=False)
        ccd.header.add_history('Trimmed image to ' + self.trim_section)
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
                o_ccd = self.image_overscan(ccd)
                to_ccd = self.image_trim(o_ccd)
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
                t_ccd = self.image_trim(ccd)
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
                ccd = self.image_overscan(ccd)
                ccd = self.image_trim(ccd)
            elif self.technique == 'Imaging':
                ccd = self.image_trim(ccd)
                ccd = ccdproc.subtract_bias(ccd, self.master_bias, add_keyword=False)
            else:
                log.error('Unknown observation technique: ' + self.technique)
            master_flat_list.append(ccd)
        master_flat = ccdproc.combine(master_flat_list,
                                      method='median',
                                      sigma_clip=True,
                                      sigma_clip_low_thresh=1.0,
                                      sigma_clip_high_thresh=1.0,
                                      add_keyword=False)
        master_flat.write(master_flat_name, clobber=True)
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

    def process_spectroscopy_science(self, science_group):
        """

        Args:
            science_group:

        Returns:

        """

        target_name = ''
        obstype = science_group.obstype.unique()
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
                    log.warning('Adding ')

            for science_image in object_group.file.tolist():
                self.out_prefix = ''
                ccd = CCDData.read(os.path.join(self.args.raw_path, science_image), unit=u.adu)
                ccd = self.image_overscan(ccd)
                self.out_prefix += 'o_'
                ccd = self.image_trim(ccd)
                self.out_prefix = 't' + self.out_prefix
                if not self.args.ignore_bias:
                    ccd = ccdproc.subtract_bias(ccd=ccd, master=self.master_bias, add_keyword=False)
                    self.out_prefix = 'z' + self.out_prefix
                    ccd.header.add_history('Bias subtracted image')
                if master_flat is None or master_flat_name is None:
                    log.warning('The file ' + science_image + ' will not be flatfielded')
                else:
                    ccd = ccdproc.flat_correct(ccd=ccd, flat=master_flat, add_keyword=False)
                    self.out_prefix = 'f' + self.out_prefix
                    ccd.header.add_history('Flat corrected ' + master_flat_name.split('/')[-1])
                # print('Length CCD Data', ccd.data.shape)
                if self.args.clean_cosmic:
                    ccd = self.cosmicray_rejection(ccd=ccd)
                    self.out_prefix = 'c' + self.out_prefix
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
                ccd = self.image_trim(ccd)
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
        value = 0.16 * float(ccd.header['EXPTIME']) + 1.2
        log.info('Cleaning cosmic rays... ')
        ccd.data, _ = ccdproc.cosmicray_lacosmic(ccd.data, sigclip=2.5, sigfrac=value, objlim=value,
                                             gain=float(ccd.header['GAIN']),
                                             readnoise=float(ccd.header['RDNOISE']),
                                             satlevel=np.inf, sepmed=True, fsmode='median',
                                             psfmodel='gaussy', verbose=True)
        # ccd.data = np.array(ccd.data, dtype=np.double) / float(ccd.header['GAIN'])

        # difference = np.nan_to_num(ccd - nccd)

        # plt.imshow(difference, clim=(0,1))
        # plt.show()

        ccd.header.add_history("Cosmic rays rejected with LACosmic")
        log.info("Cosmic rays rejected with LACosmic")
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

