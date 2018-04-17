from ccdproc import (CCDData, ImageFileCollection)
from scipy import (signal)
import os
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import (models, fitting)
from astropy.convolution import (convolve, Gaussian1DKernel)
import re
import glob
import logging


# import goodman
from pipeline.core import (read_fits,
                           write_fits,
                           SpectroscopicMode)

from pipeline.wcs import WCS

plt.rcParams["figure.figsize"] = [16, 9]


class DataValidationError(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)


class WavelengthCalibration(object):

    def __init__(self, data_path=None):
        """The class is initialized only with the path to the data
        then instance is called for each image."""
        assert os.path.isdir(data_path)
        self.log = logging.getLogger(__name__)
        self.path = data_path
        self.comparison_lamps = None
        self.image_collection = ImageFileCollection(self.path)
        self.spectroscopic_mode = SpectroscopicMode()
        self.wcs = WCS()
        self.record_lines = False

    def __call__(self, ccd, lamps=None, **kwargs):
        self.record_lines = kwargs.get('record_lines', False)
        # print(file_name)
        # if os.path.isabs(file_name) and os.path.isfile(file_name):
        #     print("is abs path")
        #     self._wavelength_solution(ccd_data=file_name, lamps=lamps)
        # full_path = os.path.join(self.path, file_name)
        # if os.path.isfile(full_path):
        #     print("file exist")
        ccd = self._wavelength_solution(ccd_data=ccd, lamps=lamps)
        return ccd
        # else:
        #     print("Can't locate file: {:s}".format(os.path.join(self.path,
        #                                                        file_name)))

    @staticmethod
    def _add_angstrom_lines_slot(ccd, line_list):
        """Adds keyords for registering Angstrom values"""
        assert isinstance(ccd, CCDData)
        assert isinstance(line_list, list)
        for i in range(len(line_list)):
            # self.log.info("GSP_A{:03d}".format(i + 1), line_list[i])
            ccd.header.set("GSP_A{:03d}".format(i + 1),
                           value=0,
                           comment="Line location in angstrom value")
        return ccd

    @staticmethod
    def _add_pixel_lines(ccd, line_list):
        assert isinstance(ccd, CCDData)
        assert isinstance(line_list, list)
        for i in range(len(line_list)):
            print("GSP_P{:03d}".format(i + 1), line_list[i])
            ccd.header.set("GSP_P{:03d}".format(i + 1),
                           value=line_list[i],
                           comment="Line location in pixel value")
        return ccd

    def _get_wavelength_solution(self, comp_lamps, record_lines=False):
        # MAKE SURE THIS WILL RECORD THE LINES ONLY ON COMPARISON LAMPS
        assert isinstance(comp_lamps, list)
        assert all([isinstance(lamp, CCDData) for lamp in comp_lamps])
        all_lamps = []
        for lamp in comp_lamps:
            lamp_lines = self._identify_lines(ccd=lamp)
            if record_lines:
                print("Recording Lines")
                lamp = self._add_pixel_lines(ccd=lamp, line_list=lamp_lines)
                lamp = self._add_angstrom_lines_slot(ccd=lamp,
                                                     line_list=lamp_lines)

                new_name = 'll_' + lamp.header['GSP_FNAM']

                write_fits(ccd=lamp,
                           full_path=os.path.join(self.path, new_name),
                           parent_file=lamp.header['GSP_FNAM'])
            all_lamps.append(lamp)
        return all_lamps

    def _identify_lines(self, ccd, show_plots=False):
        assert isinstance(ccd, CCDData)

        raw_pixel_axis = range(ccd.data.size)
        nccd = ccd.copy()
        nccd.data = np.nan_to_num(nccd.data)

        filtered_data = np.where(
            np.abs(nccd.data > nccd.data.min() +
                   0.03 * nccd.data.max()),
            nccd.data,
            None)

        # replace None to zero and convert it to an array
        none_to_zero = [0 if it is None else it for it in filtered_data]
        filtered_data = np.array(none_to_zero)

        _upper_limit = nccd.data.min() + 0.03 * nccd.data.max()
        slit_size = np.float(re.sub('[a-zA-Z" ]', '', nccd.header['SLIT']))

        serial_binning, parallel_binning = [
            int(x) for x in nccd.header['CCDSUM'].split()]

        new_order = int(round(float(slit_size) / (0.15 * serial_binning)))
        # self.log.debug('New Order:  {:d}'.format(new_order))

        print("New Order: ", round(new_order))
        peaks = signal.argrelmax(filtered_data, axis=0, order=new_order)[0]
        # ines_center = peaks

        if slit_size >= 5.:

            lines_center = self._recenter_broad_lines(lamp_data=nccd.data,
                                                      lines=peaks,
                                                      order=new_order)
        else:
            lines_center = self._recenter_lines(ccd=nccd, lines=peaks)

        # print(len(peaks), len(lines_center))
        if show_plots:
            fig = plt.figure(1)
            fig.canvas.set_window_title('Lines Detected')
            plt.title('Lines detected in Lamp\n'
                      '{:s}'.format(ccd.header['OBJECT']))
            plt.xlabel('Pixel Axis')
            plt.ylabel('Intensity (counts)')

            # Build legends without data to avoid repetitions
            plt.plot([], color='k', label='Comparison Lamp Data')

            plt.plot([],
                     color='k',
                     linestyle=':',
                     label='Spectral Line Detected')

            plt.axhline(_upper_limit, color='r')

            for line in peaks:
                plt.axvline(line, color='k', linestyle=':')
            for rline in lines_center:
                plt.axvline(rline, color='r', alpha=0.5, linestyle='-.')

            # plt.axhline(median + stddev, color='g')
            # for rc_line in lines_center:
            #     plt.axvline(rc_line, color='r')

            plt.plot(raw_pixel_axis, nccd.data, color='k')  # , marker='o')
            # plt.plot(filtered_data, color='b')
            plt.legend(loc='best')
            plt.show()

        return lines_center

    @staticmethod
    def _recenter_broad_lines(lamp_data, lines, order):
        """Recenter broad lines

        Notes:
            This method is used to recenter broad lines only, there is a special
            method for dealing with narrower lines.

        Args:
            lamp_data (array): numpy.ndarray instance. It contains the lamp
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

    @staticmethod
    def _recenter_lines(ccd, lines, plots=False):
        """Finds the centroid of an emission line

        For every line center (pixel value) it will scan left first until the
        data stops decreasing, it assumes it is an emission line and then will
        scan right until it stops decreasing too. Defined those limits it will
        use the line data in between and calculate the centroid.

        Notes:
            This method is used to recenter relatively narrow lines only, there
            is a special method for dealing with broad lines.

        Args:
            ccd (object): a ccdproc.CCDData instance.
            lines (list): A line list in pixel values.
            plots (bool): If True will plot spectral line as well as the input
              center and the recentered value.

        Returns:
            A list containing the recentered line positions.

        """
        new_center = []
        x_size = ccd.data.shape[0]
        median = np.median(ccd.data)
        for line in lines:
            # TODO (simon): Check if this definition is valid, so far is not
            # TODO (cont..): critical
            left_limit = 0
            right_limit = 1
            condition = True
            left_index = int(line)

            while condition and left_index - 2 > 0:

                if (ccd.data[left_index - 1] > ccd.data[left_index]) and \
                        (ccd.data[left_index - 2] > ccd.data[left_index - 1]):

                    condition = False
                    left_limit = left_index

                elif ccd.data[left_index] < median:
                    condition = False
                    left_limit = left_index

                else:
                    left_limit = left_index

                left_index -= 1

            # id right limit
            condition = True
            right_index = int(line)
            while condition and right_index + 2 < x_size - 1:

                if (ccd.data[right_index + 1] > ccd.data[right_index]) and \
                        (ccd.data[right_index + 2] > ccd.data[right_index + 1]):

                    condition = False
                    right_limit = right_index

                elif ccd.data[right_index] < median:
                    condition = False
                    right_limit = right_index

                else:
                    right_limit = right_index
                right_index += 1
            index_diff = [abs(line - left_index), abs(line - right_index)]

            sub_x_axis = range(line - min(index_diff),
                               (line + min(index_diff)) + 1)

            sub_data = ccd.data[line - min(index_diff):(line + min(index_diff))
                                                       + 1]
            centroid = np.sum(sub_x_axis * sub_data) / np.sum(sub_data)

            # checks for asymmetries
            differences = [abs(ccd.data[line] - ccd.data[left_limit]),
                           abs(ccd.data[line] - ccd.data[right_limit])]

            if max(differences) / min(differences) >= 2.:
                if plots:
                    plt.axvspan(line - 1, line + 1, color='g', alpha=0.3)
                new_center.append(line)
            else:
                new_center.append(centroid)
        if plots:
            fig = plt.figure(1)
            fig.canvas.set_window_title('Lines Detected in Lamp')
            plt.axhline(median, color='b')

            plt.plot(range(len(ccd.data)),
                     ccd.data,
                     color='k',
                     marker='o',
                     label='Lamp Data')

            for line in lines:

                plt.axvline(line + 1,
                            color='k',
                            linestyle=':',
                            label='First Detected Center')

            for center in new_center:

                plt.axvline(center,
                            color='k',
                            linestyle='.-',
                            label='New Center')

            plt.show()
        return new_center

    def _to_header(self, ccd):
        """Add wavelength solution to the new FITS header

        Defines FITS header keyword values that will represent the wavelength
        solution in the header so that the image can be read in any other
        astronomical tool. (e.g. IRAF)

        Notes:
            This method also saves the data to a new FITS file, This should be
            in separated methods to have more control on either process.

        Args:
            ccd (CCDData): CCDData instance with 1D data.

        Returns:
            ccd (object): An instance of `ccdproc.CCDData` with an updated
              header.

        """
        rms_error, n_points, n_rejections = self.evaluate_solution()

        ccd.header.set('GSP_WRMS', value=rms_error)
        ccd.header.set('GSP_WPOI', value=n_points)
        ccd.header.set('GSP_WREJ', value=n_rejections)

        new_crpix = 1
        new_crval = ccd.data[0][new_crpix - 1]
        new_cdelt = ccd.data[0][new_crpix] - ccd.data[0][new_crpix - 1]

        ccd.header.set('BANDID1',
                       value='spectrum - background none, weights none, '
                             'clean no', comment='')
        # ccd.header.set(['APNUM1'] = '1 1 1452.06 1454.87'
        ccd.header.set('WCSDIM',
                       value=1,
                       comment='')
        ccd.header.set('CTYPE1',
                       value='LINEAR  ',
                       comment='')
        ccd.header.set('CRVAL1',
                       value=new_crval,
                       comment='')
        ccd.header.set('CRPIX1',
                       value=new_crpix,
                       comment='')
        ccd.header.set('CDELT1',
                       value=new_cdelt,
                       comment='')
        ccd.header.set('CD1_1',
                       value=new_cdelt,
                       comment='')
        ccd.header.set('LTM1_1', value=1.,
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
        print(self.calibration_lamp)
        ccd.header.set('DCLOG1',
                       value='REFSPEC1 = {:s}'.format(self.calibration_lamp),
                       comment='')

        return ccd

    @staticmethod
    def _validate_input(ccd_data):
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

    def _wavelength_solution(self, ccd_data, lamps=None):

        assert lamps is None or isinstance(lamps, list)
        lamps_list = None

        try:
            ccd = self._validate_input(ccd_data=ccd_data)
            if lamps is not None:
                lamps_list = []
                for lamp in lamps:
                    print(lamp)
                    try:
                        ccd_lamp = self._validate_input(ccd_data=lamp)
                        lamps_list.append(ccd_lamp)
                    except DataValidationError as error:
                        print("DataValidationError: Lamp ", error)
        except DataValidationError as error:
            print("DataValidationError:", error)
            return

        if ccd.header['OBSTYPE'] == 'OBJECT':
            print('OBJECT')
            if lamps_list is not None:
                print("Lamp is a list or CCData or a filename")
                self._get_wavelength_solution(comp_lamps=lamps_list,
                                              record_lines=self.record_lines)
            else:
                print("No lamps provided")
                print("save non-calibrated object")

                # lamps_list = self.get_comparison_lamp(data=ccd)

            print("Call wavelength calibration finder...")

        elif ccd.header['OBSTYPE'] == 'COMP':
            print('Input is a comparison Lamp COMP')
            ccd_list = self._get_wavelength_solution(
                comp_lamps=[ccd],
                record_lines=self.record_lines)
            return ccd_list[0]
        else:
            print(ccd.header['OBSTYPE'])


if __name__ == '__main__':
    # path = '/user/simon/data/soar/work/aller/2017-06-11/RED'
    path = '/user/simon/data/soar/comp_lamp_lib/work/goodman_comp'
    # path = '/user/simon/data/soar/work/2018-02-08_David/comp'
    # file_1 = 'gcfzsto_0131_CuHeAr_G1200M2_slit103.fits'
    # file_2 = 'cfzsto_0131_CuHeAr_G1200M2_slit103.fits'
    # file_3 = 'gcfzsto_0144_Abell36_G1200M2_slit103.fits'

    # files = ["ext_cfzsto_0029-0033_goodman_comp_400M1_CuHeAr.fits",
    #          "ext_cfzsto_0041-0045_goodman_comp_930M3_GG385_CuHeAr.fits",
    #          "ext_cfzsto_0046-0050_goodman_comp_930M6_OG570_CuHeAr.fits",
    #          "ext_cfzsto_0050-0054_goodman_comp_600Red_GG495_CuHeAr.fits",
    #          "ext_cfzsto_0051-0056_goodman_comp_930M6_OG570_CuHeAr.fits",
    #          "ext_cfzsto_0063-0067_goodman_comp_1200M6_GG495_CuHeAr.fits",
    #          "ext_cfzsto_0068-0072_goodman_comp_930M4_GG495_CuHeAr.fits",
    #          "ext_cfzsto_0077-0081_goodman_comp_400M2_GG455_CuHeAr.fits",
    #          "ext_cfzsto_0108-0112_goodman_comp_600UV_CuHeAr.fits",
    #          "ext_cfzsto_0119-0123_goodman_comp_1200M7_OG570_CuHeAr.fits",
    #          "ext_cfzsto_0149-0153_goodman_comp_930M5_GG495_CuHeAr.fits"]

    # lamps_list = ['gcfzsto_0137_CuHeAr_G1200M2_slit103.fits',
    #               'gcfzsto_0157_CuHeAr_G1200M2_slit103.fits',
    #               'gcfzsto_0174_CuHeAr_G1200M2_slit103.fits',
    #               'gcfzsto_0142_CuHeAr_G1200M3_slit103.fits',
    #               'gcfzsto_0161_CuHeAr_G1200M2_slit103.fits',
    #               'gcfzsto_0177_CuHeAr_G1200M2_slit103.fits']
    #
    # lamps = [os.path.join(path, lamp_file) for lamp_file in lamps_list]

    files = glob.glob(path + "/ext_*fits")

    wavelength_calibration = WavelengthCalibration(data_path=path)
    for file in files:
        # full_path = os.path.join(path, file)
        # wavelength_calibration(full_path)
        print(file)
        wavelength_calibration(file, record_lines=True)
