import sys
import matplotlib.pyplot as plt
import glob
import os
import re
from ccdproc import ImageFileCollection
from ccdproc import CCDData
from astropy import units as u
import matplotlib.pyplot as plt

sys.path.append('../..')

from pipeline.core import (read_fits,
                           write_fits,
                           identify_targets,
                           trace,
                           extraction,
                           SpectroscopicMode,
                           combine_data)

from pipeline.spectroscopy import (get_args, WavelengthCalibration)

from pipeline.wcs import WCS
from astropy.modeling import (models, fitting)


class CombineAndExtract(object):

    def __init__(self):
        self.data_path = None
        self.output_path = None
        self.image_collection = None
        self.spec_mode = SpectroscopicMode()
        self.gaussian_profile = models.Gaussian1D(amplitude=10000,
                                                  mean=1,
                                                  stddev=30)
        self.target_trace = models.Polynomial1D(degree=2)

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
        self.wcs = WCS()

    def __call__(self, data_path, output_path=None, glob_include='*.fits'):
        self.data_path = data_path
        if output_path is None:
            self.output_path = data_path
        else:
            self.output_path = output_path
        image_collection = ImageFileCollection(location=self.data_path,
                                               keywords=self.keywords,
                                               glob_include=glob_include)
        if image_collection is not None:

            self.image_collection = image_collection.summary.to_pandas()

            grouped_data = self._classify_images(ic=self.image_collection)

            for group in grouped_data:
                file_list = group.file.tolist()
                for _file in file_list:
                    print(_file)
                output_name = self._get_combined_name(file_list=file_list,
                                                      prefix='comb_')

                combined = self._combine_data(file_list=file_list,
                                              output_name=output_name,
                                              data_path=self.data_path,
                                              out_path=self.output_path)

                extracted_name = self._get_extracted_name(combined_name=output_name)

                if combined.header['OBSTYPE'] == 'COMP':
                    extracted = self._extract_lamp(ccd=combined,
                                                   output_name=extracted_name)
                # ext_copy = extracted.copy()

                # ADD COPY AND THEN ADD THE GSP WAY OF WCS TO THE HEADER

                # self._create_plot(ccd=extracted, x_label="Dispersion (pixels)",
                #                   y_label='Intensity (ADU)')

                # ws_model = self._wavelength_calibration(ccd=extracted)
                #
                # if ws_model is not None:
                #     nccd = self.wcs.write_gsp_wcs(ccd=ext_copy, model=ws_model)
                #
                #     gsp_name = os.path.join(self.output_path,
                #                             re.sub('ext_', 'gsp_', extracted_name))
                #     print(gsp_name)
                #     nccd.write(gsp_name, overwrite=True)


    def _create_plot(self, ccd, x_label='', y_label=''):
        plt.title("Lamp: {:s}\nGrating: {:s}\nSlit: {:s}".format(
            ccd.header['OBJECT'],
            ccd.header['WAVMODE'],
            ccd.header['SLIT']))
        plt.plot(ccd.data)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.tight_layout()
        plt.show()

    def _wavelength_calibration(self, ccd):
        args = get_args()
        wavelength_calibration = WavelengthCalibration(args=args)
        wavelength_calibration(ccd, [ccd])
        return wavelength_calibration.wsolution

    def _extract_lamp(self, ccd, output_name):
        assert isinstance(ccd, CCDData)
        assert ccd.header['OBSTYPE'] == 'COMP'

        spatial, dispersion = ccd.data.shape

        # print(spatial, dispersion)

        center = 0.5 * spatial

        fifteen_percent = 0.15 * spatial
        fake_fwhm = 0.2 * spatial

        self.gaussian_profile.mean.value = center
        self.gaussian_profile.stddev.vale = fake_fwhm

        self.target_trace.c0.value = center

        extracted = extraction(ccd=ccd,
                               target_trace=self.target_trace,
                               spatial_profile=self.gaussian_profile,
                               extraction_name='fractional')


        # print(center - fifteen_percent, center + fifteen_percent)

        # plt.imshow(ccd.data, clim=(10, 80))
        # plt.axhspan(center - fifteen_percent, center + fifteen_percent, alpha=0.3, color='r')
        # plt.plot(self.target_trace(range(dispersion)), color='g')
        # plt.show()
        # return 0
        print("Saving File: {:s}".format(output_name))
        write_fits(ccd=extracted, full_path=os.path.join(self.output_path,
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
        grouped_data = []
        groups = ic.groupby(
            ['gain',
             'rdnoise',
             'roi',
             'grating',
             'filter2',
             'wavmode',
             'obstype',
             'object',
             'exptime']).size().reset_index().rename(columns={0: 'count'})

        for i in groups.index:
            # print(pandas_ic.iloc[i])
            # print("\nGROUP {:d}".format(i))
            # print(groups.iloc[i]['grating'])
            # print(groups.iloc[i]['filter2'])
            # print(groups.iloc[i]['object'] )
            # print(groups.iloc[i]['wavmode'])
            # print(groups.iloc[i]['exptime'])
            # print(groups.iloc[i]['obstype'])
            sub_group = ic[
                ((ic['grating'] == groups.iloc[i]['grating']) &
                 (ic['filter2'] == groups.iloc[i]['filter2']) &
                 (ic['object'] == groups.iloc[i]['object']) &
                 (ic['wavmode'] == groups.iloc[i]['wavmode']) &
                 (ic['exptime'] == groups.iloc[i]['exptime']) &
                 (ic['obstype'] == groups.iloc[i]['obstype']))]
            grouped_data.append(sub_group)

        return grouped_data


if __name__ == '__main__':
    # data_path = '/user/simon/data/soar/comp_lamp_lib/work/'
    # destination = '/user/simon/data/soar/comp_lamp_lib/work/goodman_comp'
    #
    # directories = ["2017-12-19_GSP/RED"]
    data_path = '/user/simon/data/soar/work/2018-02-08_David/'
    destination = '/user/simon/data/soar/work/2018-02-08_David/comp'

    directories = ["RED"]
    combine_and_extract = CombineAndExtract()

    for folder in directories:
        path = os.path.join(data_path, folder)
        print("Processing folder: {:s}".format(path))
        combine_and_extract(data_path=path,
                            output_path=destination,
                            glob_include='cfzsto*fits')