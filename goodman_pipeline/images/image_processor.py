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
from astropy.io import fits
from ccdproc import CCDData
from ..core import (astroscrappy_lacosmic,
                    call_cosmic_rejection,
                    combine_data,
                    create_master_bias,
                    create_master_flats,
                    define_trim_section,
                    get_best_flat,
                    get_slit_trim_section,
                    get_overscan_region,
                    image_overscan,
                    image_trim,
                    name_master_flats,
                    normalize_master_flat,
                    read_fits,
                    write_fits)

from ..core import SaturationValues, SpectroscopicMode

log = logging.getLogger(__name__)


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
        self.queue = []
        self.trim_section = None
        self.overscan_region = None
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
                image_list = group[0]['file'].tolist()
                sample_image = os.path.join(self.args.raw_path,
                                            random.choice(image_list))

                self.trim_section = define_trim_section(
                    sample_image=sample_image,
                    technique=self.technique)

                self.overscan_region = get_overscan_region(
                    sample_image=sample_image,
                    technique=self.technique)

                for sub_group in group:
                    group_obstype = sub_group.obstype.unique()

                    if len(group_obstype) == 1 and \
                        group_obstype[0] == 'BIAS' and \
                            not self.args.ignore_bias:

                        log.debug('Creating Master Bias')
                        bias_files = sub_group.file.tolist()
                        self.master_bias, self.master_bias_name = \
                            create_master_bias(
                                bias_files=bias_files,
                                raw_data=self.args.raw_path,
                                reduced_data=self.args.red_path,
                                technique=self.technique)

                    elif len(group_obstype) == 1 and group_obstype[0] == 'FLAT':
                        log.debug('Create Master FLATS')
                        flat_files = sub_group.file.tolist()
                        sample_header = fits.getheader(os.path.join(
                            self.args.raw_path, flat_files[0]))
                        master_flat_name = name_master_flats(
                            header=sample_header,
                            technique=self.technique,
                            reduced_data=self.args.red_path,
                            sun_set=self.sun_set,
                            sun_rise=self.sun_rise,
                            evening_twilight=self.evening_twilight,
                            morning_twilight=self.morning_twilight)

                        create_master_flats(
                            flat_files=flat_files,
                            raw_data=self.args.raw_path,
                            reduced_data=self.args.red_path,
                            technique=self.technique,
                            overscan_region=self.overscan_region,
                            trim_section=self.trim_section,
                            master_bias_name=self.master_bias_name,
                            new_master_flat_name=master_flat_name,
                            saturation_threshold=self.args.saturation_threshold,
                            ignore_bias=self.args.ignore_bias)
                    else:
                        log.debug('Process Data Group')
                        if self.technique == 'Spectroscopy':
                            self.process_spectroscopy_science(sub_group)
                        else:
                            log.info('Processing Imaging Science Data')
                            self.process_imaging_science(sub_group)

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
        valid_obstypes = ['OBJECT', 'COMP', 'SPECTRUM', 'ARC']
        # print(obstype)
        if any([value in valid_obstypes for value in obstype]):
            object_comp_group = science_group[(
                (science_group.obstype == 'OBJECT') |
                (science_group.obstype == 'COMP') |
                (science_group.obstype == 'SPECTRUM') |
                (science_group.obstype == 'ARC'))]

            if any([value in ['OBJECT', 'SPECTRUM'] for value in obstype]):
                target_name = science_group.object[(
                        (science_group.obstype == 'OBJECT') |
                        (science_group.obstype == 'SPECTRUM'))].unique()[0]

                log.info('Processing Science Target: '
                              '{:s}'.format(target_name))
            else:
                # TODO (simon): This does not make any sense
                log.info('Processing Comparison Lamp: '
                              '{:s}'.format(target_name))

            if any([value in ['FLAT', 'LAMPFLAT'] for value in obstype]) and \
                    not self.args.ignore_flats:
                flat_sub_group = science_group[(
                        (science_group.obstype == 'FLAT') |
                        (science_group.obstype == 'LAMPFLAT'))]

                flat_files = flat_sub_group.file.tolist()
                sample_header = fits.getheader(os.path.join(
                    self.args.raw_path, flat_files[0]))
                master_flat_name = name_master_flats(
                    header=sample_header,
                    technique=self.technique,
                    reduced_data=self.args.red_path,
                    sun_set=self.sun_set,
                    sun_rise=self.sun_rise,
                    evening_twilight=self.evening_twilight,
                    morning_twilight=self.morning_twilight)

                master_flat, master_flat_name = \
                    create_master_flats(
                            flat_files=flat_files,
                            raw_data=self.args.raw_path,
                            reduced_data=self.args.red_path,
                            technique=self.technique,
                            overscan_region=self.overscan_region,
                            trim_section=self.trim_section,
                            master_bias_name=self.master_bias_name,
                            new_master_flat_name=master_flat_name,
                            saturation_threshold=self.args.saturation_threshold,
                            ignore_bias=self.args.ignore_bias)
            elif self.args.ignore_flats:
                log.warning('Ignoring creation of Master Flat by request.')
                master_flat = None
                master_flat_name = None
            else:
                log.info('Attempting to find a suitable Master Flat')
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
                    master_flat_name = name_master_flats(
                        header=ccd.header,
                        technique=self.technique,
                        reduced_data=self.args.red_path,
                        sun_set=self.sun_set,
                        sun_rise=self.sun_rise,
                        evening_twilight=self.evening_twilight,
                        morning_twilight=self.morning_twilight,
                        get=True)

                    # load the best flat based on the name previously defined
                    master_flat, master_flat_name = \
                        get_best_flat(flat_name=master_flat_name,
                                      path=self.args.red_path)
                    if (master_flat is None) and (master_flat_name is None):
                        log.critical('Failed to obtain master flat')

            if master_flat is not None and not self.args.ignore_flats:
                if not self.args.skip_slit_trim:
                    log.debug('Attempting to find slit trim section')
                    slit_trim = get_slit_trim_section(master_flat=master_flat)
                    log.debug('Slit trim section found: {}'.format(slit_trim))
                else:
                    log.warning('Skipping slit trim section trimming')
            elif self.args.ignore_flats:
                log.warning('Slit Trimming will be skipped, '
                                 '--ignore-flats is activated')
            else:
                log.info("Master flat inexistent, can't find slit trim "
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
                    master_bias = image_trim(ccd=self.master_bias,
                                             trim_section=self.trim_section,
                                             trim_type='trimsec')
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
                    print(ccd.data.shape)
                    print(master_bias.data.shape)

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
                    log.warning('Ignoring bias correction by request.')

                # Do flat correction
                if master_flat is None or master_flat_name is None:
                    log.warning('The file {:s} will not be '
                                     'flatfielded'.format(science_image))
                elif self.args.ignore_flats:
                    log.warning('Ignoring flatfielding by request.')
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
                    keep_files=self.args.keep_cosmic_files,
                    method=self.args.clean_cosmic,
                    save=True)

                self.out_prefix = prefix

                if ccd is not None:
                    if ccd.header['OBSTYPE'] in ['OBJECT', 'SPECTRUM']:
                        log.debug("Appending {} image for combination".format(
                            ccd.header['OBSTYPE']))
                        all_object_image.append(ccd)
                    elif ccd.header['OBSTYPE'] in ['COMP', 'ARC']:
                        log.debug("Appending {} image for combination".format(
                            ccd.header['OBSTYPE']))
                        all_comp_image.append(ccd)
                    else:
                        log.error("Unknown OBSTYPE = {:s}"
                                  "".format(ccd.header['OBSTYPE']))
                else:
                    log.warning("Cosmic ray rejection returned a None.")
            if self.args.combine:
                log.warning("Combination of data is experimental.")
                if len(all_object_image) > 1:
                    # print(len(all_object_image))
                    log.info("Combining {:d} OBJECT images"
                             "".format(len(all_object_image)))

                    # object_group = object_comp_group[(
                    #     (object_comp_group.obstype == "OBJECT") |
                    #     (object_comp_group.obstype == "SPECTRUM"))]

                    # print(object_group, len(all_object_image))

                    combined_data = combine_data(all_object_image,
                                                 dest_path=self.args.red_path,
                                                 prefix=self.out_prefix,
                                                 save=True)

                elif len(all_object_image) == 1:
                    # write_fits(all_object_image[0])
                    log.critical("Unable to combine one single image")

                else:
                    log.error("No OBJECT images to combine")

                if len(all_comp_image) > 1:
                    log.info("Combining {:d} COMP images"
                                  "".format(len(all_comp_image)))
                    # comp_group = object_comp_group[
                    #     object_comp_group.obstype == "COMP"]
                    combine_data(all_comp_image,
                                 dest_path=self.args.red_path,
                                 prefix=self.out_prefix,
                                 save=True)
                else:
                    log.error("No COMP images to combine")
            else:
                log.debug("Combine is disabled (default)")

        elif any([value in ['FLAT', 'LAMPFLAT'] for value in obstype]):
            self.queue.append(science_group)
            log.warning('Only flats found in this group')
            flat_sub_group = science_group[(
                    (science_group.obstype == 'FLAT') |
                    (science_group.obstype == 'LAMPFLAT'))]
            # TODO (simon): Find out if these variables are useful or not
            flat_files = flat_sub_group.file.tolist()

            sample_header = fits.getheader(os.path.join(self.args.raw_path,
                                                        flat_files[0]))
            master_flat_name = name_master_flats(
                header=sample_header,
                technique=self.technique,
                reduced_data=self.args.red_path,
                sun_set=self.sun_set,
                sun_rise=self.sun_rise,
                evening_twilight=self.evening_twilight,
                morning_twilight=self.morning_twilight)

            create_master_flats(
                            flat_files=flat_files,
                            raw_data=self.args.raw_path,
                            reduced_data=self.args.red_path,
                            technique=self.technique,
                            overscan_region=self.overscan_region,
                            trim_section=self.trim_section,
                            master_bias_name=self.master_bias_name,
                            new_master_flat_name=master_flat_name,
                            saturation_threshold=self.args.saturation_threshold)
        else:
            log.error('There is no valid datatype in this group')

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

        master_flat_name = name_master_flats(
            header=sample_file.header,
            technique=self.technique,
            reduced_data=self.args.red_path,
            sun_set=self.sun_set,
            sun_rise=self.sun_rise,
            evening_twilight=self.evening_twilight,
            morning_twilight=self.morning_twilight,
            get=True)

        log.debug('Got {:s} for master flat name'.format(master_flat_name))

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
                log.info('Created science file: {:s}'.format(final_name))
        else:
            log.error('Can not process data without a master flat')


if __name__ == '__main__':
    pass
