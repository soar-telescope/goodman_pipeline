from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.modeling import (models, fitting)
from ccdproc import CCDData
from ccdproc import ImageFileCollection


from pipeline.tools.reference_lamp_factory.line_matcher import LineMatcher
from pipeline.tools.reference_lamp_factory.wcs_model_calculator import \
    GSPWcsCalculator

from ...core import (read_fits,
                     write_fits,
                     identify_targets,
                     trace,
                     extraction,
                     SpectroscopicMode,
                     combine_data)

from pipeline.spectroscopy.new_wavelength import WavelengthCalibration

import argparse
import json
import logging
import matplotlib.pyplot as plt
import os
import sys


class SettingsField(object):

    def __init__(self, field):
        for attr in field:
            if field[attr] == 'none':
                self.__setattr__(attr, None)
            else:
                self.__setattr__(attr, field[attr])


class Settings(object):
    def __init__(self, json_settings):
        for field in json_settings:
            self.__setattr__(field, SettingsField(json_settings[field]))


class ReferenceLibraryFactory(object):

    def __init__(self, arguments=None):
        self.log = logging.getLogger(__name__)
        self.args = self._get_args(arguments=arguments)
        self.settings = self._load_settings()
        self.line_identify = WavelengthCalibration(data_path=self.args.path)
        self.line_match = LineMatcher()
        self.wcs_calc = GSPWcsCalculator(save_pdf_to=self.args.save_to)
        self.image_collection = None
        self.gaussian_profile = None
        self.target_trace = None
        self.target_location = None
        self.target_fwhm = None
        self.keywords = ['date',
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
                         'roi',
                         'wavmode']

        for field in self.settings:
            self.__setattr__(field, SettingsField(self.settings[field]))

    def __call__(self, search_pattern=None, keywords=None, **kwargs):

        gaussian_amplitude = self.target.gaussian_amplitude
        gaussian_mean = self.target.gaussian_mean
        gaussian_stddev = self.target.gaussian_stddev

        trace_poly_degree = self.target.trace_poly_degree
        trace_poly_degree = int(trace_poly_degree)

        self.target_location = gaussian_mean
        self.target_fwhm = gaussian_stddev

        self.gaussian_profile = models.Gaussian1D(amplitude=gaussian_amplitude,
                                                  mean=gaussian_mean,
                                                  stddev=gaussian_stddev)

        self.target_trace = models.Polynomial1D(degree=trace_poly_degree)

        if keywords is not None and isinstance(list, keywords):
            assert all([isinstance(key, str) for key in keywords])
            self.keywords = keywords

        for ccd in self._classify_and_combine():
            print(ccd.data.shape)
            # print(combined.header['NAXIS'])
            ccd = self.line_identify(ccd=ccd, record_lines=True)
            ccd = self.line_match(ccd=ccd)

            ccd = self.wcs_calc(ccd=ccd)

            file_name = ccd.header['GSP_FNAM']
            full_path = os.path.join(self.args.save_to, file_name)
            write_fits(ccd=ccd, full_path=full_path)
            # plt.title(combined.header['GSP_FNAM'])
            # plt.plot(combined.data)
            # plt.show()

    def _classify_and_combine(self):

        search_pattern = self.prefix.search_pattern
        comb_prefix = self.prefix.combined_prefix
        extracted_prefix = self.prefix.extracted_prefix

        image_collection = ImageFileCollection(location=self.args.path,
                                               keywords=self.keywords,
                                               glob_include=search_pattern)

        if image_collection is not None:

            try:

                self.image_collection = image_collection.summary.to_pandas()

            except AttributeError as error:

                if str(error) == "'NoneType' object has no attribute " \
                                 "'to_pandas'":
                    self.log.critical("Folder {:s} does not contain "
                                      "data.".format(self.args.path))
                    raise SystemExit

            grouped_data = self._classify_images(ic=self.image_collection)

            for group in grouped_data:

                file_list = group.file.tolist()

                for _file in file_list:
                    self.log.info(_file)

                output_name = self._get_combined_name(file_list=file_list,
                                                      prefix=comb_prefix)

                combined = self._combine_data(file_list=file_list,
                                              output_name=output_name,
                                              data_path=self.args.path,
                                              out_path=self.args.save_to)

                extracted_name = self._get_extracted_name(
                    combined_name=output_name)

                if combined.header['OBSTYPE'] == 'COMP':
                    extracted = self._extract_lamp(ccd=combined,
                                                   output_name=extracted_name)
                    yield extracted

    @staticmethod
    def _create_plot(ccd, x_label='', y_label=''):
        plt.title("Lamp: {:s}\nGrating: {:s}\nSlit: {:s}".format(
            ccd.header['OBJECT'],
            ccd.header['WAVMODE'],
            ccd.header['SLIT']))
        plt.plot(ccd.data)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.tight_layout()
        plt.show()

    def _extract_lamp(self, ccd, output_name):
        assert isinstance(ccd, CCDData)
        assert ccd.header['OBSTYPE'] == 'COMP'

        spatial, dispersion = ccd.data.shape

        # print(spatial, dispersion)

        if self.target_location is None or self.target_location == 'none':
            self.target_location = 0.5 * spatial

        if self.target_fwhm is None or self.target_location == 'none':
            self.target_fwhm = 0.2 * spatial

        self.gaussian_profile.mean.value = self.target_location
        self.gaussian_profile.stddev.vale = self.target_fwhm

        self.target_trace.c0.value = self.target_location

        extracted = extraction(ccd=ccd,
                               target_trace=self.target_trace,
                               spatial_profile=self.gaussian_profile,
                               extraction_name='fractional')

        print("Saving File: {:s}".format(output_name))
        write_fits(ccd=extracted, full_path=os.path.join(self.args.save_to,
                                                         output_name))
        return extracted

    @staticmethod
    def _combine_data(file_list, output_name, data_path, out_path):
        image_list = []
        for image in file_list:
            ccd = read_fits(os.path.join(data_path, image))
            image_list.append(ccd)
        try:
            combined = combine_data(image_list=image_list, dest_path=out_path,
                                    output_name=output_name, save=True)

            return combined
        except AssertionError as error:
            print(error)
            return image_list[0]

    @staticmethod
    def _get_combined_name(file_list, prefix=''):
        number_range = "{:s}-{:s}".format(file_list[0].split('_')[1],
                                          file_list[-1].split('_')[1])
        split_name = file_list[0].split('_')
        split_name[1] = number_range

        output_name = prefix + '_'.join(split_name)

        return output_name

    @staticmethod
    def _get_extracted_name(combined_name, prefix=''):
        base_name = os.path.basename(combined_name)

        split_base_name = base_name.split('_')

        split_base_name[0] = 'ext'

        extracted_name = prefix + '_'.join(split_base_name)

        return extracted_name

    @staticmethod
    def _classify_images(ic):
        classification_keywords = ['gain',
                                   'rdnoise',
                                   'roi',
                                   'grating',
                                   'filter2',
                                   'wavmode',
                                   'obstype',
                                   'object',
                                   'exptime']
        grouped_data = []
        groups = ic.groupby(
            classification_keywords
            ).size().reset_index().rename(columns={0: 'count'})

        for i in groups.index:
            sub_group = ic[
                ((ic['grating'] == groups.iloc[i]['grating']) &
                 (ic['filter2'] == groups.iloc[i]['filter2']) &
                 (ic['object'] == groups.iloc[i]['object']) &
                 (ic['wavmode'] == groups.iloc[i]['wavmode']) &
                 (ic['exptime'] == groups.iloc[i]['exptime']) &
                 (ic['obstype'] == groups.iloc[i]['obstype']))]
            grouped_data.append(sub_group)

        return grouped_data

    @staticmethod
    def _get_args(arguments=None):
        parser = argparse.ArgumentParser(description='Tool to create reference '
                                                     'lamps for using with the '
                                                     'Goodman Spectroscopic '
                                                     'Pipeline.')

        parser.add_argument('--from',
                            action='store',
                            default=os.getcwd(),
                            dest='path',
                            help="Path to raw data.")

        parser.add_argument('--save-to',
                            action='store',
                            default=os.getcwd(),
                            dest='save_to',
                            help='Where to save new data.')

        parser.add_argument('--settings',
                            action='store',
                            default='settings.json',
                            dest='json_settings',
                            help='Name of the json files containing the'
                                 'settings')

        args = parser.parse_args(args=arguments)
        if not os.path.exists(args.path) or not os.path.isdir(args.path):
            parser.exit("Input folder {:s} does not exist".format(args.path))
        if not os.path.exists(args.save_to) or not os.path.isdir(args.save_to):
            parser.exit("Output folder {:s} does not "
                        "exist".format(args.save_to))
        return args

    def _load_settings(self):
        """Load json file with settings

        First search for it in the data folder and then search for the default
        one

        Raises:
            AttributeError
        """
        json_settings_full_path = os.path.join(self.args.path,
                                               self.args.json_settings)

        default_settings_file = os.path.join(os.path.dirname(__file__),
                                             self.args.json_settings)

        if os.path.isfile(json_settings_full_path):
            file_full_path = json_settings_full_path
        elif os.path.isfile(default_settings_file):
            file_full_path = default_settings_file
        else:
            raise AttributeError('Can not find settings file.')

        with open(file_full_path) as json_settings:
            settings = json.load(json_settings)
            return settings




if __name__ == '__main__':
    factory = ReferenceLibraryFactory()
    factory()
