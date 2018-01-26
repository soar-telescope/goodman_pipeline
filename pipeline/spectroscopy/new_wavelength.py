from ccdproc import (CCDData, ImageFileCollection)
from scipy import (signal)
import os
import sys
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.io.fits.header import Header
from astropy.modeling import (models, fitting, Model)
import re

import sys

# import goodman
from pipeline.core import (read_fits, write_fits, interpolate)

plt.rcParams["figure.figsize"] = [16, 9]


class DataValidationError(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)


def validate_input(ccd_data):
    assert isinstance(ccd_data, str) or isinstance(ccd_data, CCDData)

    if isinstance(ccd_data, CCDData):
        ccd = ccd_data.copy()
    elif not isinstance(ccd_data, CCDData):
        ccd = read_fits(ccd_data)
    else:
        raise NotImplementedError(
            "Can't process {:s}".format(str(type(ccd_data))))

    if not isinstance(ccd, CCDData):
        raise DataValidationError(
            'Invalid Input: Data is not an instance of CCDData.')
    elif ccd.header['NAXIS'] != 1:
        raise DataValidationError("Wrong data dimensions")

    return ccd


def identify_lines(lamp):
    print("identifying liines in lamp ", lamp.header['GSP_FNAM'])

    interpolation_size = 200

    x_axis, data = interpolate(lamp.data, interpolation_size)

    filtered_data = np.where(
        np.abs(data > data.min() + 0.03 * data.max()),
        data,
        np.zeros(data.shape))

    peaks = signal.argrelmax(filtered_data, axis=0, order=7 * 200)[0]

    slit_size = interpolation_size * float(re.sub('[a-zA-Z" ]', '', lamp.header['slit']))
    print(slit_size)
    for peak in peaks:
        recenter_line(data=data, center=peak, slit_size=slit_size)
        print(x_axis[peak], peak)
        # plt.axvline(x_axis[peak], color='r')
    plt.plot(lamp.data)
    plt.show()


def get_wavelength_solution(comp_lamps):
    assert isinstance(comp_lamps, list)
    assert all([isinstance(lamp, CCDData) for lamp in comp_lamps])

    for lamp in comp_lamps:
        lamp_lines = identify_lines(lamp=lamp)


def validate_emission_line(data, line_center, left_limit, right_limit,slit_size):

    data_sample = data[left_limit:right_limit]
    data_x_axis = np.linspace(left_limit, right_limit, len(data_sample))
    model_fitter = fitting.LevMarLSQFitter()
    box_width = slit_size / 0.15
    slit_model = models.Box1D(amplitude=data[int(line_center)],
                              x_0=line_center,
                              width=box_width)
    gaussian_model = models.Gaussian1D(amplitude=data[int(line_center)],
                                       mean=line_center,
                                       stddev=right_limit - left_limit)

    fitted_gaussian = model_fitter(gaussian_model * slit_model,
                                   data_x_axis,
                                   data_sample)

    final_model = fitted_gaussian
    print(final_model)

    plt.title("Data and Model")
    plt.plot(data_x_axis, data_sample, label="Data")
    plt.plot(data_x_axis, final_model(data_x_axis), label="Final Model")
    # plt.xlim(left_limit - .75 * (right_limit - left_limit),
    #          right_limit + .75 * (right_limit - left_limit))
    # plt.ylim(- .3 * data.max(), data.max() + .3 * data.max())
    plt.legend(loc='best')
    plt.show()


def recenter_line(data, center, slit_size):
    int_center = int(round(center))
    peak_value = data[int_center]

    estimated_fwhm = 0.5 * peak_value

    filtered = np.where(np.abs(data < 0.5 * peak_value), data,
                        np.zeros(data.shape))
    # plt.plot(data)
    # plt.plot(filtered)
    # plt.show()

    left_side_pix = int_center
    right_side_pix = int_center
    print("Center ", int_center)
    i = 1
    while left_side_pix == int_center or right_side_pix == int_center:
        left = int_center - i
        right = int_center + i

        # print(filtered[left], filtered[right])
        if filtered[left] != 0 and left > 0 and left_side_pix == int_center:
            # print("Left ", left)
            left_side_pix = left
        elif left <= 0 and left_side_pix == int_center:
            left_side_pix = 0
        if filtered[right] != 0 and right < len(
                data) - 1 and right_side_pix == int_center:
            # print("Right ", right)
            right_side_pix = right
        else:
            pass
            # print(left, right)
        i += 1
        # print(i)
        # if left_side_pix != int_center and right_side_pix != int_center:
        #     break


    print(filtered[left_side_pix], filtered[right_side_pix])
    left_weight = filtered[left_side_pix] / (
    filtered[left_side_pix] + filtered[right_side_pix])
    right_weight = filtered[right_side_pix] / (
    filtered[left_side_pix] + filtered[right_side_pix])
    print(left_weight, right_weight)
    # print(1/(estimated_fwhm - filtered[left_side_pix]), 1/(estimated_fwhm - filtered[right_side_pix]))

    estimated_center = np.mean([left_side_pix, right_side_pix])
    weighted_center = np.average([left_side_pix, right_side_pix],
                                 weights=[left_weight, right_weight])

    validate_emission_line(data=data,
                           line_center=estimated_center,
                           left_limit=left_side_pix,
                           right_limit=right_side_pix,
                           slit_size=slit_size)

    print(estimated_center)

    plt.title("Initial Center: {:.3f}\nCorrected Center: {:.3f}".format(center,
                                                                        estimated_center))
    if left_side_pix != right_side_pix:
        plt.axvline(left_side_pix, color='c')
        plt.axvline(right_side_pix, color='c')

    plt.plot(data, color='b')
    plt.axvline(center, color='r', alpha=.4)
    plt.axhline(peak_value, color='g')
    plt.axhline(estimated_fwhm, color='y')
    plt.plot(filtered)
    plt.xlim(center - .75 * (right_side_pix - left_side_pix),
             center + .75 * (right_side_pix - left_side_pix))
    plt.ylim(- .3 * peak_value, peak_value + .3 * peak_value)
    plt.axvline(weighted_center, color='m', alpha=.4,
                label='Weighted Center {:.3f}'.format(weighted_center))
    plt.axvline(estimated_center, color='k',
                label="Estimated Center {:.3f}".format(estimated_center))
    plt.legend(loc='best')
    plt.show()

    return weighted_center


class WavelengthCalibration(object):
    def __init__(self, data_path=None):
        assert os.path.isdir(data_path)
        self.path = data_path
        self.image_collection = ImageFileCollection(self.path)
        self.comparison_lamps = None

    def __call__(self, file_name, lamps=None):
        print(file_name)
        if os.path.isabs(file_name) and os.path.isfile(file_name):
            print("is abs path")
            self.wavelength_solution(ccd_data=file_name, lamps=lamps)
        elif os.path.isfile(os.path.join(self.path, file_name)):
            print("file exist")
            full_path = os.path.join(self.path, file_name)
            self.wavelength_solution(ccd_data=full_path, lamps=lamps)
        else:
            print("Can't locate file: {:s}".format(os.path.join(self.path,
                                                                file_name)))

    def wavelength_solution(self, ccd_data, lamps=None):

        assert lamps is None or isinstance(lamps, list)
        lamps_list = None

        try:
            ccd = validate_input(ccd_data=ccd_data)
            if lamps is not None:
                lamps_list = []
                for lamp in lamps:
                    print(lamp)
                    ccd_lamp = validate_input(ccd_data=lamp)
                    lamps_list.append(ccd_lamp)
        except DataValidationError as error:
            print("DataValidationError:", error)
            return

        if ccd.header['OBSTYPE'] == 'OBJECT':
            print('OBJECT')
            if lamps_list is not None:
                print("Lamp is a list or CCData or a filename")
                get_wavelength_solution(comp_lamps=lamps_list)
            else:
                print("No lamps provided")
                print("save non-calibrated object")

                # lamps_list = self.get_comparison_lamp(data=ccd)

            print("Call wavelength calibration finder...")

        elif ccd.header['OBSTYPE'] == 'COMP':
            print('Input is a comparison Lamp COMP')
            get_wavelength_solution(comp_lamps=[ccd])
        # if lamps_list is not None:
        #                 print("Warning: lamps will be treated as independent")
        #                 for ccd_lamp in lamps_list:
        #                     plt.plot(ccd_lamp.data)
        #                     plt.title("Lamp: {:s}".format(ccd_lamp.header['GSP_FNAM']))
        #                     plt.show()
        #             print("Call wavelength calibration finder...")
        else:
            print(ccd.header['OBSTYPE'])


if __name__ == '__main__':
    path = '/user/simon/data/soar/work/aller/2017-06-11/RED'
    file_1 = 'gcfzsto_0131_CuHeAr_G1200M2_slit103_1.fits'
    file_2 = 'cfzsto_0131_CuHeAr_G1200M2_slit103.fits'
    file_3 = 'gcfzsto_0144_Abell36_G1200M2_slit103.fits'

    lamps_list = ['gcfzsto_0131_CuHeAr_G1200M2_slit103_1.fits',
                  'gcfzsto_0170_CuHeAr_G1200M2_slit103.fits',
                  'gcfzsto_0174_CuHeAr_G1200M2_slit103.fits']

    lamps = [os.path.join(path, lamp_file) for lamp_file in lamps_list]

    wavelength_calibration = WavelengthCalibration(data_path=path)
    full_path = os.path.join(path, file_1)
    wavelength_calibration(full_path)
    wavelength_calibration(file_3)
