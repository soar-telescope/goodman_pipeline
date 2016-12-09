"""Contains the tools to produce a wavelength solution

This module gets the extracted data to produce a wavelength solution, linerize the spectrum and write the solution
to the image's header following the FITS standard.
"""
from __future__ import print_function
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip
import scipy.interpolate
import logging as log
import matplotlib.image as mpimg
import matplotlib.ticker as mtick
from scipy import signal
import time
from linelist import ReferenceData
import wsbuilder

FORMAT = '%(levelname)s:%(filename)s:%(module)s: 	%(message)s'
log.basicConfig(level=log.INFO, format=FORMAT)


class WavelengthCalibration(object):

    def __init__(self, sci_pack, science_object, args):
        self.args = args
        self.wsolution = None
        self.rms_error = None
        self.reference_data = ReferenceData(self.args)
        self.science_object = science_object
        self.slit_offset = 0
        self.interpolation_size = 200
        self.line_search_method = 'derivative'
        """Instrument configuration and spectral characteristics"""
        self.gratings_dict = {'SYZY_400': 400,
                              'KOSI_600': 600,
                              '930': 930,
                              'RALC_1200-BLUE': 1200}
        self.grating_frequency = 0
        self.grating_angle = float(0)
        self.camera_angle = float(0)
        self.binning = 1
        self.pixel_count = 0
        self.alpha = 0
        self.beta = 0
        self.center_wavelength = 0
        self.blue_limit = 0
        self.red_limit = 0
        """Interactive wavelength finding"""
        self.reference_clicks_x = []
        self.reference_clicks_y = []
        self.raw_data_clicks_x = []
        self.raw_data_clicks_y = []
        self.click_input_enabled = True
        self.reference_bb = None
        self.raw_data_bb = None
        self.contextual_bb = None
        self.i_fig = None
        self.ax1 = None
        self.ax2 = None
        self.ax3 = None
        self.ax4 = None
        self.ax4_plots = None
        self.ax4_com = None
        self.ax4_rlv = None
        self.points_ref = None
        self.points_raw = None
        self.line_raw = None
        self.filling_value = 1000
        self.events = True
        self.first = True
        self.evaluation_comment = None
        # self.binning = self.lamp_header[]
        self.pixelcenter = []
        """this data must come parsed"""
        self.path = self.args.source
        self.science_pack = sci_pack
        self.sci_filename = self.science_object.file_name
        # self.history_of_lamps_solutions = {}
        self.reference_solution = None

    def __call__(self, wsolution_obj=None):
        log.info('Processing Science Target: %s', self.science_pack.headers[0]['OBJECT'])
        if wsolution_obj is None:
            if self.science_object.lamp_count > 0:
                for lamp_index in range(self.science_object.lamp_count):
                    self.calibration_lamp = self.science_object.lamp_file[lamp_index - 1]
                    self.lamp_data = self.science_pack.lamps_data[lamp_index]
                    self.raw_pixel_axis = range(1, len(self.lamp_data) + 1, 1)
                    # self.raw_pixel_axis = range(len(self.lamp_data))
                    self.lamp_header = self.science_pack.lamps_headers[lamp_index]
                    self.lamp_name = self.lamp_header['OBJECT']
                    log.info('Processing Comparison Lamp: %s', self.lamp_name)
                    self.data1 = self.interpolate(self.lamp_data)
                    # self.lines_limits = self.get_line_limits()
                    # self.lines_center = self.get_line_centers(self.lines_limits)
                    self.lines_center = self.get_lines_in_lamp()
                    self.spectral = self.get_spectral_characteristics()
                    if self.args.interactive_ws:
                        self.interactive_wavelength_solution()
                    else:
                        log.warning('Automatic Wavelength Solution is not fully implemented yet')
                        self.automatic_wavelength_solution()
                        # self.wsolution = self.wavelength_solution()
                    if self.wsolution is not None:
                        self.linear_lamp = self.linearize_spectrum(self.lamp_data)
                        self.lamp_header = self.add_wavelength_solution(self.lamp_header,
                                                                        self.linear_lamp,
                                                                        self.science_object.lamp_file[lamp_index - 1])
                        for target_index in range(self.science_object.no_targets):
                            log.debug('Processing target %s', target_index + 1)
                            new_data = self.science_pack.data[target_index]
                            new_header = self.science_pack.headers[target_index]
                            if self.science_object.no_targets > 1:
                                new_index = target_index + 1
                            else:
                                new_index = None
                            self.linearized_sci = self.linearize_spectrum(new_data)
                            self.header = self.add_wavelength_solution(new_header,
                                                                       self.linearized_sci,
                                                                       self.sci_filename,
                                                                       index=new_index)
                        wavelength_solution = WavelengthSolution(solution_type='non_linear',
                                                                 model_name='chebyshev',
                                                                 model_order=3,
                                                                 model=self.wsolution,
                                                                 ref_lamp=self.calibration_lamp,
                                                                 eval_comment=self.evaluation_comment,
                                                                 header=self.header)

                        return wavelength_solution
                    else:
                        log.error('It was not possible to get a wavelength solution from this lamp.')
                        return None

            else:
                log.error('There are no lamps to process')
        else:
            self.wsolution = wsolution_obj.wsolution
            self.calibration_lamp = wsolution_obj.reference_lamp
            self.evaluation_comment = wsolution_obj.evaluation_comment
            # print('wavelengthSolution ', self.wsolution)
            # print('Evaluation Comment', self.evaluation_comment)
            # repeat for all sci
            for target_index in range(self.science_object.no_targets):
                log.debug('Processing target %s', target_index + 1)
                new_data = self.science_pack.data[target_index]
                new_header = self.science_pack.headers[target_index]
                if self.science_object.no_targets > 1:
                    new_index = target_index + 1
                else:
                    new_index = None
                self.linearized_sci = self.linearize_spectrum(new_data)
                self.header = self.add_wavelength_solution(new_header,
                                                           self.linearized_sci,
                                                           self.sci_filename,
                                                           self.evaluation_comment,
                                                           index=new_index)

    def get_wsolution(self):
        """Get the mathematical model of the wavelength solution

        The wavelength solution is a callable mathematical function from astropy.modeling.models
        By obtaining this solution it can be applied to a pixel axis.

        Returns:
            wsolution (callable): A callable mathematical function

        """
        if self.wsolution is not None:
            return self.wsolution
        else:
            log.error("Wavelength Solution doesn't exist!")
            return None

    def get_calibration_lamp(self):
        """Get the name of the calibration lamp used for obtain the solution

        The filename of the lamp used to obtain must go to the header for documentation

        Returns:
            calibration_lamp (str): Filename of calibration lamp used to obtain wavelength solution

        """
        if self.wsolution is not None and self.calibration_lamp is not None:
            return self.calibration_lamp
        else:
            log.error('Wavelength solution has not been calculated yet.')

    def get_lines_in_lamp(self):
        """Identify peaks in a lamp spectrum

        This method goes through all the data points in a 1D spectrum recording whether there is a trend
        to increase or to decrease. It also records the length of the trend. If there was rising trend
        it will check if trend length is 3 or larger and if the point is above the median of the spectrum.
        If these criteria are met the last rising pixel value is recorded as a line candidate.

        Notes:
            This method will only work for comparison lamps.

        Returns:
            lines_candidates (list): A common list containing pixel values at approximate location of lines.
        """
        # lines_candidates = []
        # prev_point = None
        # status = None
        # trend_length = 0
        # median = np.median(self.lamp_data)
        # stddev = np.std(self.lamp_data)
        #
        # for pixel_index in range(len(self.lamp_data)):
        #     point = self.lamp_data[pixel_index]
        #     # print(point)
        #     if prev_point is None:
        #         prev_point = point
        #         status = 0
        #     else:
        #         if point > prev_point:
        #             if status == 1:
        #                 trend_length += 1
        #             elif status == -1:
        #                 trend_length = 1
        #             status = 1
        #         elif point < prev_point:
        #             if status == 1 and trend_length > 2 and point > median:
        #                 lines_candidates.append(pixel_index - 1)
        #             status = -1
        #             trend_length = 1
        #         else:
        #             pass
        #     prev_point = point

        filtered_data = np.where(np.abs(self.lamp_data > self.lamp_data.min() + 0.05 * self.lamp_data.max()),
                                        self.lamp_data,
                                        0)
        peaks = signal.argrelmax(filtered_data, axis=0, order=6)[0]
        for value in peaks:
            plt.axvline(value, color='r')

        lines_center = self.recenter_lines(self.lamp_data, peaks)

        if self.args.plots_enabled:
            fig = plt.figure(1)
            fig.canvas.set_window_title('Lines Detected')
            plt.title('Lines detected in Lamp\n%s' % self.lamp_header['OBJECT'])
            plt.xlabel('Pixel Axis')
            plt.ylabel('Intensity (counts)')
            # Build legends without data
            plt.plot([], color='k', label='Comparison Lamp Data')
            plt.plot([], color='k', linestyle=':', label='Spectral Line Detected')
            for line in peaks:
                plt.axvline(line + 1, color='k', linestyle=':')
            # plt.axhline(median + stddev, color='g')
            plt.plot(self.raw_pixel_axis, self.lamp_data, color='k')
            plt.legend(loc='best')
            plt.show()
        return lines_center

    def recenter_lines(self, data, lines, plots=False):
        new_center = []
        x_size = data.shape[0]
        median = np.median(data)
        for line in lines:
            # id left limit
            # print(line)
            condition = True
            left_index = int(line)
            while condition and left_index - 2 > 0:
                if (data[left_index - 1] > data[left_index]) and (data[left_index - 2] > data[left_index - 1]):
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
                if (data[right_index + 1] > data[right_index]) and (data[right_index + 2] > data[right_index + 1]):
                    condition = False
                    right_limit = right_index
                elif data[right_index] < median:
                    condition = False
                    right_limit = right_index
                else:
                    right_limit = right_index
                right_index += 1
            index_diff = [abs(line - left_index), abs(line - right_index)]

            sub_x_axis = range(line - min(index_diff), (line + min(index_diff)) + 1)
            sub_data = data[line - min(index_diff):(line + min(index_diff)) + 1]
            centroid = np.sum(sub_x_axis * sub_data) / np.sum(sub_data)
            # checks for asymmetries
            differences = [abs(data[line] - data[left_limit]), abs(data[line] - data[right_limit])]
            if max(differences) / min(differences) >= 2.:
                if plots:
                    plt.axvspan(line - 1, line + 1, color='g', alpha=0.3)
                new_center.append(line + 1)
            else:
                new_center.append(centroid + 1)
            # if plots:
                # plt.axvline(centroid, color='m')
                # plt.plot(sub_x_axis, sub_data)
                # plt.plot(data, color='r', linestyle='--')
                # plt.xlim([left_limit, right_limit])
                # plt.ylim([min(sub_data), max(sub_data)])
                # plt.show()
        if plots:
            fig = plt.figure(1)
            fig.canvas.set_window_title('Lines Detected in Lamp')
            plt.axhline(median, color='b')
            plt.plot(self.raw_pixel_axis, data, color='k', label='Lamp Data')
            for line in lines:
                plt.axvline(line + 1, color='k', linestyle=':', label='First Detected Center')
            for center in new_center:
                plt.axvline(center, color='k', linestyle='.-', label='New Center')
            plt.show()
        return new_center

    # def get_line_limits(self):
    #     """Method for identifying lines in a spectrum
    #
    #     This is the spectral line identifying method. It calculates a pseudo-derivative of the spectrum thus determining
    #     peaks. For a typical spectrum a pseudo-derivative, from now on just "derivative", will produce a series of
    #     positive and negative peaks, for emission lines. Since we will be using it for comparison lamps only we don't
    #     have to worry about absorption lines and therefore this method only detects emission lines.
    #     A threshold is defined by calculating the 75 percent of the standard deviation of the derivative. There
    #     is no particular reason for choosing 75 percent is just the value that worked best for all the test subjects.
    #     Then any search for lines must be done for values above and below of the threshold and its negative
    #     respectively.  There are a few control mechanisms that ensures that an actual good line is being detected. For
    #     instance: For each line both positive and negative derivative's peaks are stored and are called limits. In order
    #     to add a (first) positive limits the number of previously stored limits must be even, meaning that we are
    #     "opening" a new line. Correspondingly if we are adding  a (last) negative limit we have to have added previously
    #     a odd amount of limits meaning that we are "closing" a line. At this point there is another constraint, it
    #     cannot be separated by an amount larger than a defined spacing, or we would be having a very broad line and lamp
    #     lines are expected to be very narrow.
    #
    #     Notes:
    #         The lines detected by this method are usually the best lines in the spectrum. There is a certain amount of
    #         lines missed but in very crowded regions.
    #
    #         Very important to note that the line centers found using the data produced here are very precise.
    #
    #     Returns:
    #         limits (list): List of line limits in consecutive pairs. This is considered by the method that uses this
    #         list for further processing. The values are index values not pixel values.
    #
    #     """
    #     if self.line_search_method == 'derivative':
    #         # in fact is a pseudo-derivative
    #         derivative = []
    #         # derivative2 = []
    #         faux_x = range(0, len(self.data1[1]) - 1)
    #         # faux_x2 = range(0, len(self.data1[1]) - 2)
    #         # print faux_x
    #         for i in faux_x:
    #             derivative.append(self.data1[1][i + 1] - self.data1[1][i])
    #         # for e in faux_x2:
    #         #    derivative2.append(derivative[e] - derivative[e+1])
    #         threshold = np.std(derivative) * .75
    #         new_range = 0
    #         spacing = 1500
    #         limits = []
    #         for i in range(len(derivative) - 1):
    #             if i > new_range:
    #                 if derivative[i] > threshold and derivative[i + 1] - derivative[i] >= 0:
    #                     partial_max = i + np.argmax(derivative[i:i + spacing])
    #                     # print i, partial_min
    #                     new_range = partial_max + (partial_max - i)
    #                     if limits == [] or len(limits) % 2 == 0:
    #                         limits.append(partial_max)
    #                     else:
    #                         plt.axvline(partial_max, color='k')
    #                 elif derivative[i] < -threshold and derivative[i + 1] - derivative[i] <= 0:
    #                     partial_min = i + np.argmin(derivative[i:i + spacing])
    #                     new_range = partial_min + (partial_min - i)
    #                     if len(limits) % 2 == 1 and partial_min - limits[-1] < spacing:
    #                         limits.append(partial_min)
    #                         # plt.axvline(partial_max, color='g')
    #                     elif limits != []:
    #                         if partial_min - limits[-1] > spacing:
    #                             plt.axvline(partial_min, color='m')
    #                             limits = limits[:-1]
    #         if len(limits) % 2 == 1:
    #             limits = limits[:-1]
    #
    #         # Produce Plots
    #         if self.args.plots_enabled:
    #             for i in range(len(limits)):
    #                 if i % 2 == 0:
    #                     plt.axvline(limits[i], color='r')
    #                 elif i % 2 == 1:
    #                     plt.axvline(limits[i], color='g')
    #
    #             plt.title('Line Identification')
    #             # plt.plot(self.data1[1], label='Spectrum')
    #             plt.plot(faux_x, derivative, label='1st Derivative')
    #             plt.axhline(0, color='m')
    #             # plt.plot(faux_x2, derivative2, label='2nd')
    #             plt.axhline(threshold)
    #             plt.axhline(-threshold)
    #             plt.legend(loc='best')
    #             # plt.plot(pixel_axis, self.data1[1])
    #             plt.savefig(self.path + 'line-identification-201.png', dpi=300)
    #             plt.show()
    #             # """
    #         return limits
    #
    # def get_line_centers(self, limits):
    #     """Finds the center of the lines using limits previously found
    #
    #     This method is very simple and could be integrated in the get_line_limits method but I'd rather have the limit
    #     information available for later. Basically finds the mean value of the line limits and then finds the
    #     correspoing pixel value, adds it up to the "centers" list.
    #
    #     Args:
    #         limits (list): Line limits in the list's index domain.
    #
    #     Returns:
    #         centers (list): Line centers in pixel values as floats.
    #
    #     """
    #     centers = []
    #     for i in range(0, len(limits), 2):
    #         center = (float(limits[i]) + float(limits[i + 1])) / 2.
    #         width = limits[i + 1] - limits[i]
    #         pixel_width = self.data1[0][limits[i + 1]] - self.data1[0][limits[i]]
    #         log.debug('Approximate FWHM: %s pix %s Angstrom (pix * 0.65)', pixel_width, pixel_width * 0.65)
    #         i_min = int(center - 2 * width)
    #         i_max = int(center + 2 * width)
    #         pixel_axis = self.data1[0][i_min:i_max]
    #         data_axis = self.data1[1][i_min:i_max]
    #         pixel_center = self.data1[0][int(round(center))]
    #         center_val = self.data1[1][int(round(center))]
    #         # new_center = self.recenter_line_by_model(pixel_axis, data_axis, center_val, pixel_center, 'gauss')
    #         if self.args.plots_enabled:
    #             plt.plot(pixel_axis, data_axis)
    #             plt.axvline(pixel_center)
    #             # plt.axvline(new_center, color='m')
    #             plt.show()
    #         self.pixelcenter.append([pixel_center, center_val])
    #         centers.append(pixel_center)
    #         # print(center, width)
    #     return centers

    def get_spectral_characteristics(self):
        """Calculates some Goodman's specific spectroscopic values.

        From the Header value for Grating, Grating Angle and Camera Angle it is possible to estimate what are the limits
        wavelength values and central wavelength. It was necessary to add offsets though, since the formulas provided
        are slightly off. The values are only an estimate.

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
        blue_correction_factor = -90
        red_correction_factor = -60
        self.grating_frequency = self.gratings_dict[self.lamp_header['GRATING']]
        self.grating_angle = float(self.lamp_header['GRT_ANG'])
        self.camera_angle = float(self.lamp_header['CAM_ANG'])
        # binning = self.lamp_header[]
        # TODO(simon): Make sure which binning is the important, parallel or serial
        # self.binning = 1
        # PG5_4 parallel
        # PG5_9 serial
        # PARAM18 serial
        # PARAM22 parallel
        try:
            self.binning = self.lamp_header['PG5_4']
            # serial_binning = self.lamp_header['PG5_9']
        except KeyError:
            self.binning = self.lamp_header['PARAM22']
            # serial_binning = self.lamp_header['PARAM18']


        # self.pixel_count = len(self.lamp_data)
        # Calculations
        self.alpha = self.grating_angle + self.slit_offset
        self.beta = self.camera_angle - self.grating_angle
        self.center_wavelength = 10 * (1e6 / self.grating_frequency) * (
            np.sin(self.alpha * np.pi / 180.) + np.sin(self.beta * np.pi / 180.))
        self.blue_limit = 10 * (1e6 / self.grating_frequency) * (
            np.sin(self.alpha * np.pi / 180.) + np.sin((self.beta - 4.656) * np.pi / 180.)) + blue_correction_factor
        self.red_limit = 10 * (1e6 / self.grating_frequency) * (
            np.sin(self.alpha * np.pi / 180.) + np.sin((self.beta + 4.656) * np.pi / 180.)) + red_correction_factor
        pixel_one = self.predicted_wavelength(1)
        pixel_two = self.predicted_wavelength(2)
        log.debug('Center Wavelength : %s Blue Limit : %s Red Limit : %s',
                  self.center_wavelength,
                  self.blue_limit,
                  self.red_limit)
        spectral_characteristics = {'center': self.center_wavelength,
                                    'blue': self.blue_limit,
                                    'red': self.red_limit,
                                    'alpha': self.alpha,
                                    'beta': self.beta,
                                    'pix1': pixel_one,
                                    'pix2': pixel_two}
        return spectral_characteristics

    def interpolate(self, spectrum):
        """Creates an interpolated version of the input spectrum

        This method creates an interpolated version of the input array, it is used mainly for a spectrum but it can
        also be used with any unidimensional array, assuming you are happy with the interpolation_size attribute
        defined for this class. The reason for doing interpolation is that it allows to find the lines and its
        respective center more precisely. The default interpolation size is 200 (two hundred) points.

        Args:
            spectrum (array): an uncalibrated spectrum or any unidimensional array.

        Returns:
            Two dimensional array containing x-axis and interpolated array. The x-axis preserves original pixel values.

        """
        x_axis = range(1, spectrum.size + 1)
        first_x = x_axis[0]
        last_x = x_axis[-1]
        new_x_axis = np.linspace(first_x, last_x, spectrum.size * self.interpolation_size)

        tck = scipy.interpolate.splrep(x_axis, spectrum, s=0)
        new_spectrum = scipy.interpolate.splev(new_x_axis, tck, der=0)
        return [new_x_axis, new_spectrum]

    def recenter_line_by_data(self, data_name, x_data):
        if data_name == 'reference':
            pseudo_center = np.argmin(abs(self.reference_solution[0] - x_data))
            reference_line_index = np.argmin(abs(self.reference_data.get_line_list_by_name(self.lamp_name) - x_data))
            reference_line_value = self.reference_data.get_line_list_by_name(self.lamp_name)[reference_line_index]
            sub_x = self.reference_solution[0][pseudo_center - 10: pseudo_center + 10]
            sub_y = self.reference_solution[1][pseudo_center - 10: pseudo_center + 10]
            center_of_mass = np.sum(sub_x * sub_y) / np.sum(sub_y)
            # print 'centroid ', center_of_mass
            # plt.figure(3)
            if self.ax4_plots is not None or self.ax4_com is not None or self.ax4_rlv is not None:
                try:
                    self.ax4.cla()
                    self.ax4.relim()
                except:
                    pass

            self.ax4.set_title('Reference Data Clicked Line')
            self.ax4.set_xlabel('Wavelength (Angstrom)')
            self.ax4.set_ylabel('Intensity (Counts)')
            self.ax4_plots = self.ax4.plot(sub_x, sub_y, color='k', label='Data')
            self.ax4_rlv = self.ax4.axvline(reference_line_value, linestyle='-', color='r', label='Reference Line Value')
            self.ax4_com = self.ax4.axvline(center_of_mass, linestyle='--', color='b', label='Centroid')
            self.ax4.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
            ### self.ax4.legend(loc=4)
            self.i_fig.canvas.draw()
            # return center_of_mass
            return reference_line_value
        elif data_name == 'raw-data':
            pseudo_center = np.argmin(abs(self.raw_pixel_axis - x_data))
            raw_line_index = np.argmin(abs(self.lines_center - x_data))
            raw_line_value = self.lines_center[raw_line_index]
            # print(raw_line_value)
            sub_x = self.raw_pixel_axis[pseudo_center - 10: pseudo_center + 10]
            sub_y = self.lamp_data[pseudo_center - 10: pseudo_center + 10]
            center_of_mass = np.sum(sub_x * sub_y) / np.sum(sub_y)
            # print 'centroid ', center_of_mass
            # plt.figure(3)
            if self.ax4_plots is not None or self.ax4_com is not None or self.ax4_rlv is not None:
                self.ax4.cla()
                self.ax4.relim()

                # except:
                #     pass
            self.ax4.set_title('Raw Data Clicked Line')
            self.ax4.set_xlabel('Pixel Axis')
            self.ax4.set_ylabel('Intensity (Counts)')
            self.ax4_plots = self.ax4.plot(sub_x, sub_y, color='k', label='Data')
            self.ax4_rlv = self.ax4.axvline(raw_line_value, linestyle='-', color='r', label='Line Center')
            self.ax4_com = self.ax4.axvline(center_of_mass, linestyle='--',color='b', label='Centroid')
            self.ax4.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
            ### self.ax4.legend(loc=4)
            self.i_fig.canvas.draw()
            # return center_of_mass
            return raw_line_value
        else:
            log.error('Unrecognized data name')

    def predicted_wavelength(self, pixel):
        alpha = self.alpha
        beta = self.beta
        # pixel_count = self.pixel_count
        binning = self.binning
        grating_frequency = self.grating_frequency
        wavelength = 10 * (1e6 / grating_frequency) * (np.sin(alpha * np.pi / 180.)
                                                       + np.sin((beta * np.pi / 180.)
                                                                + np.arctan((pixel * binning - 2048) * 0.015 / 377.2)))
        return wavelength

    def automatic_wavelength_solution(self):
        # needs:
        #   - self.sci, self.header
        #   - self.lines_center
        raise NotImplemented
        return

    def interactive_wavelength_solution(self):
        """Find the wavelength solution interactively



        """
        start = time.time()
        plt.switch_backend('GTK3Agg')
        reference_file = self.reference_data.get_reference_lamps_by_name(self.lamp_name)
        if reference_file is not None:
            log.info('Using reference file: %s', reference_file)
            reference_plots_enabled = True
            ref_data = fits.getdata(reference_file)
            ref_header = fits.getheader(reference_file)
            fits_ws_reader = wsbuilder.ReadWavelengthSolution(ref_header, ref_data)
            self.reference_solution = fits_ws_reader()
        else:
            reference_plots_enabled = False
            log.error('Please Check the OBJECT Keyword of your reference data')

        # ------- Plots -------
        plot_creation = time.time() - start
        # print('plot creation ', plot_creation)
        self.i_fig, ((self.ax1, self.ax2), (self.ax3, self.ax4)) = plt.subplots(2,
                                                                               2,
                                                                               gridspec_kw={'width_ratios': [4, 1]})
        self.i_fig.canvas.set_window_title('Science Target: %s' % self.science_object.name)
        manager = plt.get_current_fig_manager()
        manager.window.maximize()
        # manager.window.attributes('-topmost', 0)
        # self.ax1 = plt.subplot(211)
        self.ax1.set_title('Raw Data - %s' % self.lamp_name)
        self.ax1.set_xlabel('Pixels')
        self.ax1.set_ylabel('Intensity (counts)')
        self.ax1.plot([], linestyle='--', color='r', label='Detected Lines')
        for idline in self.lines_center:
            self.ax1.axvline(idline, linestyle='--', color='r')
        self.ax1.plot(self.raw_pixel_axis, self.lamp_data, color='k', label='Raw Data')
        # self.ax1.plot(self.lamp_data, color='b')
        self.ax1.set_xlim((0, len(self.lamp_data)))
        # ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
        self.ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
        ### self.ax1.legend(loc='best')
        # self.ax1.set_yscale('log')
        first_plot = time.time() - start
        # print('first plot ', first_plot)

        # self.ax3 = plt.subplot(212)
        self.ax3.set_title('Reference Data')
        self.ax3.set_xlabel('Wavelength (Angstrom)')
        self.ax3.set_ylabel('Intensity (counts)')
        # self.ax3.axvline(self.blue_limit, color='k')
        # self.ax3.axvline(self.center_wavelength, color='k')
        # self.ax3.axvline(self.red_limit, color='k')
        self.ax3.set_xlim((self.blue_limit, self.red_limit))
        self.ax3.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
        self.ax3.plot([], linestyle=':', color='r', label='Reference Line Values')
        segundo_plot = time.time() - start
        # print('segundo plot', segundo_plot)
        for rline in self.reference_data.get_line_list_by_name(self.lamp_name):
            self.ax3.axvline(rline, linestyle=':', color='r')
        if reference_plots_enabled:
            l1, = self.ax3.plot(self.reference_solution[0], self.reference_solution[1], color='k', label='Reference Lamp Data')
            # self.ax3.set_xlim((self.reference_solution[0][0], self.reference_solution[0][-1]))
            # self.ax3.set_xlim((self.blue_limit, self.red_limit))
        # self.ax3.set_yscale('log')
        ### self.ax3.legend(loc='best')

        # help plot
        self.ax2.set_title('Help')
        self.ax2.set_xticks([])
        self.ax2.set_yticks([])
        self.leg = self.ax2.legend([l1], ['xxxx'], mode='expand', frameon=False)


        # self.ax2.set_axis_off()
        # ToDo (simon): figure out what to do here.

        # self.ax2.text(1, 11, 'F1:', fontsize=15)
        # self.ax2.text(1.3, 11, 'Prints Help (remove this one)', fontsize=13)
        # self.ax2.text(1, 10, 'F2:', fontsize=15)
        # self.ax2.text(1.3, 10, 'Fit Wavelength Solution to points collected', fontsize=13)
        # self.ax2.text(1, 9, 'F3: Find new lines, use with caution!.', fontsize=15)
        # self.ax2.set_ylim((0, 12))
        # self.ax2.set_xlim((0.95, 3.5))


        # zoomed plot
        self.ax4.set_title('Contextual Information')
        self.ax4_plots = self.ax4.text(0.25, 0.25, "Goodman\nSpectrograph\nSoar Telescope", fontsize=20)
        self.ax4_com, = self.ax4.plot([])
        self.ax4_rlv, = self.ax4.plot([])
        # image = mpimg.imread('/user/simon/development/soar/goodman/soar_outside.png')
        # self.ax4_plots = self.ax4.imshow(image)

        plt.subplots_adjust(left=0.05, right=0.99, top=0.97, bottom=0.04, hspace=0.17, wspace=0.11)
        self.raw_data_bb = self.ax1.get_position()
        self.reference_bb = self.ax3.get_position()
        self.contextual_bb = self.ax4.get_position()

        if self.click_input_enabled:
            self.i_fig.canvas.mpl_connect('button_press_event', self.on_click)
            self.i_fig.canvas.mpl_connect('key_press_event', self.key_pressed)
            # print self.wsolution
            plt.show()
        return True

    def on_click(self, event):
        # print event.button
        self.events = True
        if event.button == 2:
            if event.xdata is not None and event.ydata is not None:
                figure_x, figure_y = self.i_fig.transFigure.inverted().transform((event.x, event.y))
                if self.reference_bb.contains(figure_x, figure_y):
                    # self.reference_clicks.append([event.xdata, event.ydata])
                    self.reference_clicks_x.append(self.recenter_line_by_data('reference', event.xdata))
                    self.reference_clicks_y.append(event.ydata)
                    self.update_clicks_plot('reference')
                elif self.raw_data_bb.contains(figure_x, figure_y):
                    # self.raw_data_clicks.append([event.xdata, event.ydata])
                    self.raw_data_clicks_x.append(self.recenter_line_by_data('raw-data', event.xdata))
                    self.raw_data_clicks_y.append(event.ydata)
                    self.update_clicks_plot('raw_data')
                # self.ref_click_plot.set_xdata(np.array(self.reference_clicks[:][0]))
                # self.ref_click_plot.set_ydata(np.array(self.reference_clicks[:][1]))
                # self.ref_click_plot.draw()
                else:
                    log.debug(figure_x, figure_y, 'Are not contained')
                # print 'click ', event.xdata, ' ', event.ydata, ' ', event.button
                # print event.x, event.y
            else:
                log.error('Clicked Region is out of boundaries')
        elif event.button == 3:
            if len(self.reference_clicks_x) == len(self.raw_data_clicks_x):
                self.click_input_enabled = False
                log.info('Leaving interactive mode')
            else:
                if len(self.reference_clicks_x) < len(self.raw_data_clicks_x):
                    log.info('There is %s click missing in the Reference plot',
                             len(self.raw_data_clicks_x) - len(self.reference_clicks_x))
                else:
                    log.info('There is %s click missing in the New Data plot',
                             len(self.reference_clicks_x) - len(self.raw_data_clicks_x))

    def key_pressed(self, event):
        self.events = True
        if event.key == 'f1':
            log.info('Print help regarding interactive mode')
            print("F1 : Prints Help.")
            print("F2 : Fit wavelength solution model.")
            print("F3 : Find new lines.")
            print("F4 : Evaluate solution")
            print("F6 : Linearize data (for testing not definitive)")
            print("d : deletes closest point")
            # print("l : resample spectrum to a linear dispersion axis")
            print("ctrl+d : deletes all recorded clicks")
            print("ctrl+b : Go back to previous solution (deletes automatic added points")
            print('Middle Button Click: records data location.')
            print("Right Button Click: Leaves interactive mode.")
        elif event.key == 'f2':
            log.debug('Calling function to fit wavelength Solution')
            self.fit_pixel_to_wavelength()
            self.plot_raw_over_reference()
        elif event.key == 'f3':
            if self.wsolution is not None:
                self.find_more_lines()
                self.update_clicks_plot('reference')
                self.update_clicks_plot('raw_data')
        elif event.key == 'f4':
            if self.wsolution is not None and len(self.raw_data_clicks_x) > 0:
                self.evaluate_solution(plots=True)
        elif event.key == 'd':
            figure_x, figure_y = self.i_fig.transFigure.inverted().transform((event.x, event.y))
            if self.raw_data_bb.contains(figure_x, figure_y):
                log.debug('Deleting point')
                # print abs(self.raw_data_clicks_x - event.xdata)
                closer_index = int(np.argmin(abs(self.raw_data_clicks_x - event.xdata)))
                # print 'Index ', closer_index
                if len(self.raw_data_clicks_x) == len(self.reference_clicks_x):
                    self.raw_data_clicks_x.pop(closer_index)
                    self.raw_data_clicks_y.pop(closer_index)
                    self.reference_clicks_x.pop(closer_index)
                    self.reference_clicks_y.pop(closer_index)
                    self.update_clicks_plot('reference')
                    self.update_clicks_plot('raw_data')
                else:
                    self.raw_data_clicks_x.pop(closer_index)
                    self.raw_data_clicks_y.pop(closer_index)
                    self.update_clicks_plot('raw_data')
            elif self.reference_bb.contains(figure_x, figure_y):
                log.debug('Deleting point')
                # print 'reference ', self.reference_clicks_x, self.re
                # print self.reference_clicks_x
                # print abs(self.reference_clicks_x - event.xdata)
                closer_index = int(np.argmin(abs(self.reference_clicks_x - event.xdata)))
                if len(self.raw_data_clicks_x) == len(self.reference_clicks_x):
                    self.raw_data_clicks_x.pop(closer_index)
                    self.raw_data_clicks_y.pop(closer_index)
                    self.reference_clicks_x.pop(closer_index)
                    self.reference_clicks_y.pop(closer_index)
                    self.update_clicks_plot('reference')
                    self.update_clicks_plot('raw_data')
                else:
                    self.reference_clicks_x.pop(closer_index)
                    self.reference_clicks_y.pop(closer_index)
                    self.update_clicks_plot('reference')
            elif self.contextual_bb.contains(figure_x, figure_y):
                closer_index_ref = int(np.argmin(abs(self.reference_clicks_x - event.ydata)))
                closer_index_raw = int(np.argmin(abs(self.raw_data_clicks_x - event.xdata)))
                # print(closer_index_raw, closer_index_ref)
                self.raw_data_clicks_x.pop(closer_index_raw)
                self.raw_data_clicks_y.pop(closer_index_raw)
                self.reference_clicks_x.pop(closer_index_ref)
                self.reference_clicks_y.pop(closer_index_ref)
                self.update_clicks_plot('reference')
                self.update_clicks_plot('raw_data')
                self.update_clicks_plot('eval_plots')

        elif event.key == 'f6':
            log.info('Linearize spectrum')
            if self.wsolution is not None:
                self.linearize_spectrum(self.lamp_data, plots=True)
        elif event.key == 'ctrl+b':
            log.info('Deleting automatic added points. If exist.')
            if self.raw_data_clicks_x is not [] and self.reference_clicks_x is not []:
                to_remove = []
                for i in range(len(self.raw_data_clicks_x)):
                    # print self.raw_data_clicks[i], self.filling_value
                    if self.raw_data_clicks_y[i] == self.filling_value:
                        to_remove.append(i)
                        # print to_remove
                to_remove = np.array(sorted(to_remove, reverse=True))
                if len(to_remove) > 0:
                    for index in to_remove:
                        self.raw_data_clicks_x.pop(index)
                        self.raw_data_clicks_y.pop(index)
                        self.reference_clicks_x.pop(index)
                        self.reference_clicks_y.pop(index)
                    self.update_clicks_plot('reference')
                    self.update_clicks_plot('raw_data')
                    # else:
                    # print self.raw_click_plot, self.ref_click_plot, 'mmm'
        elif event.key == 'ctrl+d':
            log.info('Deleting all recording Clicks')
            answer = raw_input('Are you sure you want to delete all clicks? only typing "No" will stop it! : ')
            if answer.lower() != 'no':
                self.reference_clicks_x = []
                self.reference_clicks_y = []
                self.raw_data_clicks_x = []
                self.raw_data_clicks_y = []
                self.update_clicks_plot('delete')
                self.plot_raw_over_reference(remove=True)
            else:
                log.info('No click was deleted this time!.')
        else:
            # print event.key
            pass

    def find_more_lines(self):
        """Method to add more lines given that a wavelength solution already exists

        This method is part of the interactive wavelength solution mechanism. If a wavelength solution exist it uses the
        line centers in pixels to estimate their respective wavelength and then search for the closest value in the list
        of reference lines for the elements in the comparison lamp. Then it filters the worst of them by doing sigma
        clipping. Finally it adds them to the class' variables that contains the list of reference points.

        Better results are obtained if the solution is already good. Visual inspection also improves final result.
        """
        new_physical = []
        new_wavelength = []
        square_differences = []
        if self.wsolution is not None:
            wlines = self.wsolution(self.lines_center)
            for i in range(len(wlines)):
                closer_index = np.argmin(abs(self.reference_data.get_line_list_by_name(self.lamp_name) - wlines[i]))
                rline = self.reference_data.get_line_list_by_name(self.lamp_name)[closer_index]
                rw_difference = wlines[i] - rline
                # print('Difference w - r ', rw_difference, rline)
                square_differences.append(rw_difference ** 2)
                new_physical.append(self.lines_center[i])
                new_wavelength.append(rline)
            clipped_differences = sigma_clip(square_differences, sigma=2, iters=3)
            if len(new_wavelength) == len(new_physical) == len(clipped_differences):
                for i in range(len(new_wavelength)):
                    if clipped_differences[i] is not np.ma.masked and new_wavelength[i] not in self.reference_clicks_x:
                        self.reference_clicks_x.append(new_wavelength[i])
                        self.reference_clicks_y.append(self.filling_value)
                        self.raw_data_clicks_x.append(new_physical[i])
                        self.raw_data_clicks_y.append(self.filling_value)
        return True

    def update_clicks_plot(self, action=None, pixel_axis=None, differences=None, **kwargs):
        # print(type(action), type(pixel_axis), type(differences))
        if action == 'reference':
            if self.points_ref is not None:
                try:
                    self.points_ref.remove()
                    self.ax3.relim()
                except:
                    pass
            self.points_ref, = self.ax3.plot(self.reference_clicks_x,
                                             self.reference_clicks_y,
                                             linestyle='None',
                                             marker='o',
                                             color='r')
            self.i_fig.canvas.draw()
        elif action == 'raw_data':
            # print self.points_raw
            # print dir(self.points_raw)
            if self.points_raw is not None:
                try:
                    self.points_raw.remove()
                    self.ax1.relim()
                except:
                    pass
            self.points_raw, = self.ax1.plot(self.raw_data_clicks_x,
                                             self.raw_data_clicks_y,
                                             linestyle='None',
                                             marker='o',
                                             color='r')
            self.i_fig.canvas.draw()
        elif action == 'delete':
            if self.points_raw is not None and self.points_ref is not None:
                self.points_raw.remove()
                self.ax1.relim()
                self.points_ref.remove()
                self.ax3.relim()
                self.i_fig.canvas.draw()

    def plot_raw_over_reference(self, remove=False):
        if self.wsolution is not None:
            if self.line_raw is not None:
                try:
                    self.line_raw.remove()
                    self.ax3.relim()
                except:
                    pass
            if not remove:
                # TODO(simon): catch TypeError Exception and correct what is causing it
                self.line_raw, = self.ax3.plot(self.wsolution(self.raw_pixel_axis),
                                               self.lamp_data,
                                               linestyle='-',
                                               color='r',
                                               alpha=0.4,
                                               label='New Wavelength Solution')
            # changed yesterday
            ### self.ax3.legend(loc='best')
            self.i_fig.canvas.draw()

    def evaluate_solution(self, plots=False):
        if self.wsolution is not None:
            differences = np.array([])
            wavelength_line_centers = self.wsolution(self.lines_center)

            for wline in wavelength_line_centers:
                closer_index = np.argmin(abs(self.reference_data.get_line_list_by_name(self.lamp_name) - wline))
                rline = self.reference_data.get_line_list_by_name(self.lamp_name)[closer_index]
                rw_difference = wline - rline
                # print 'Difference w - r ', rw_difference, rline
                differences = np.append(differences, rw_difference)
                # differences.append(rw_difference)
            # differences = np.array(differences)
            clipping_sigma = 1.5
            # print(differences)
            clipped_differences = sigma_clip(differences, sigma=clipping_sigma, iters=5, cenfunc=np.ma.median)
            once_clipped_differences = sigma_clip(differences, sigma=clipping_sigma, iters=1, cenfunc=np.ma.median)

            # sigma_value = clipped_differences.std()
            # data_mean = clipped_differences.mean()
            npoints = len(clipped_differences)
            n_rejections = np.ma.count_masked(clipped_differences)
            square_differences = []
            for i in range(len(clipped_differences)):
                if clipped_differences[i] is not np.ma.masked:
                    square_differences.append(clipped_differences[i] ** 2)
            old_rms_error = None
            if self.rms_error is not None:
                old_rms_error = float(self.rms_error)
            self.rms_error = np.sqrt(np.sum(square_differences) / len(square_differences))
            log.info('RMS Error : %s', self.rms_error)
            if plots:
                if self.ax4_plots is not None or self.ax4_com is not None or self.ax4_rlv is not None:
                    try:
                        self.ax4.cla()
                        self.ax4.relim()
                    except:
                        pass
                self.ax4.set_title('RMSE %.3f \n %s points and %s rejections' % (self.rms_error, npoints, n_rejections))
                self.ax4.set_ylim(once_clipped_differences.min(), once_clipped_differences.max())
                # self.ax4.set_ylim(- rms_error, 2 * rms_error)
                self.ax4.set_xlim(np.min(self.lines_center), np.max(self.lines_center))
                self.ax4_rlv = self.ax4.scatter(self.lines_center, differences, marker='x', color='k', label='Removed Points')
                # print(data_mean - clipping_sigma / 2. * sigma_value)
                # print(data_mean + clipping_sigma / 2. * sigma_value)
                # print(clipped_differences.min(), )
                self.ax4_com = self.ax4.axhspan(clipped_differences.min(),
                                                clipped_differences.max(),
                                                color='k',
                                                alpha=0.4,
                                                label='%s Sigma' % clipping_sigma)
                self.ax4_plots = self.ax4.scatter(self.lines_center, clipped_differences, label='Differences')
                if self.rms_error is not None and old_rms_error is not None:
                    # increment_color = 'white'
                    rms_error_difference = self.rms_error - old_rms_error
                    if rms_error_difference > 0.001:
                        increment_color = 'red'
                    elif rms_error_difference < -0.001:
                        increment_color = 'green'
                    else:
                        increment_color = 'white'
                    self.ax4.text(0.05, 0.95,
                                  '%+.3f' % rms_error_difference,
                                  verticalalignment='top',
                                  horizontalalignment='left',
                                  transform=self.ax4.transAxes,
                                  color=increment_color,
                                  fontsize=15)

                #         ax.text(0.95, 0.01, 'colored text in axes coords',
                # verticalalignment='bottom', horizontalalignment='right',
                # transform=ax.transAxes,
                # color='green', fontsize=15)
                # self.ax4.axhline(rms_error, color='r')
                self.ax4.set_xlabel('Pixel Axis (Pixels)')
                self.ax4.set_ylabel('Difference (Angstroms)')

                ### self.ax4.legend(loc=4)
                self.i_fig.canvas.draw()

            return [self.rms_error, npoints, n_rejections]
        else:
            log.error('Solution is still non-existent!')

    def fit_pixel_to_wavelength(self):
        if len(self.reference_clicks_x) and len(self.raw_data_clicks_x) > 0:
            pixel = []
            angstrom = []
            for i in range(len(self.reference_clicks_x)):
                pixel.append(self.raw_data_clicks_x[i])
                angstrom.append(self.reference_clicks_x[i])
            wavelength_solution = wsbuilder.WavelengthFitter(model='chebyshev', degree=3)
            self.wsolution = wavelength_solution.ws_fit(pixel, angstrom)
            # changed yesterday
            self.evaluate_solution(plots=True)
        else:
            log.error('Clicks record is empty')
            if self.wsolution is not None:
                self.wsolution = None

    def linearize_spectrum(self, data, plots=False):
        pixel_axis = range(1, len(data) + 1, 1)
        if self.wsolution is not None:
            x_axis = self.wsolution(pixel_axis)
            new_x_axis = np.linspace(x_axis[0], x_axis[-1], len(data))
            tck = scipy.interpolate.splrep(x_axis, data, s=0)
            linearized_data = scipy.interpolate.splev(new_x_axis, tck, der=0)
            if plots:
                fig6 = plt.figure(6)
                fig6.canvas.set_window_title('Linearized Data')
                plt.plot(x_axis, data, color='b')
                plt.plot(new_x_axis, linearized_data, color='r', linestyle=':')
                plt.tight_layout()
                plt.show()
                fig7 = plt.figure(7)
                fig7.canvas.set_window_title('Wavelength Solution')
                plt.plot(x_axis, color='b')
                plt.plot(new_x_axis, color='r')
                plt.tight_layout()
                plt.show()

            return [new_x_axis, linearized_data]

    def add_wavelength_solution(self, new_header, spectrum, original_filename, evaluation_comment=None, index=None):
        if evaluation_comment is None:
            rms_error, n_points, n_rejections = self.evaluate_solution()
            self.evaluation_comment = 'Lamp Solution RMSE = %s Npoints = %s, NRej = %s' % (rms_error,
                                                                                           n_points,
                                                                                           n_rejections)
            new_header['HISTORY'] = self.evaluation_comment
        else:
            new_header['HISTORY'] = evaluation_comment

        new_crpix = 1
        new_crval = spectrum[0][new_crpix - 1]
        new_cdelt = spectrum[0][new_crpix] - spectrum[0][new_crpix - 1]

        new_header['BANDID1'] = 'spectrum - background none, weights none, clean no'
        # new_header['APNUM1'] = '1 1 1452.06 1454.87'
        new_header['WCSDIM'] = 1
        new_header['CTYPE1'] = 'LINEAR  '
        new_header['CRVAL1'] = new_crval
        new_header['CRPIX1'] = new_crpix
        new_header['CDELT1'] = new_cdelt
        new_header['CD1_1'] = new_cdelt
        new_header['LTM1_1'] = 1.
        new_header['WAT0_001'] = 'system=equispec'
        new_header['WAT1_001'] = 'wtype=linear label=Wavelength units=angstroms'
        new_header['DC-FLAG'] = 0
        new_header['DCLOG1'] = 'REFSPEC1 = %s' % self.calibration_lamp

        # print(new_header['APNUM*'])
        if index is None:
            f_end = '.fits'
        else:
            f_end = '_%s.fits' % index
        # idea
        #  remove .fits from original_filename
        # define a base original name
        # modify in to _1, _2 etc in case there are multitargets
        # add .fits

        new_filename = self.args.destiny + self.args.output_prefix + original_filename.replace('.fits', '') + f_end

        fits.writeto(new_filename, spectrum[1], new_header, clobber=True)
        # print new_header
        return new_header


class WavelengthSolution(object):
    def __init__(self,
                 solution_type=None,
                 model_name=None,
                 model_order=0,
                 model=None,
                 ref_lamp=None,
                 eval_comment='',
                 header=None):
        self.dtype_dict = {None: -1, 'linear': 0, 'log_linear': 1, 'non_linear': 2}
        # if solution_type == 'non_linear' and model_name is not None:
        self.ftype_dict = {'chebyshev': 1,
                           'legendre': 2,
                           'cubic_spline': 3,
                           'linear_spline': 4,
                           'pixel_coords': 5,
                           'samples_coords': 6,
                           None: None}
        self.solution_type = solution_type
        self.model_name = model_name
        self.model_order = model_order
        self.wsolution = model
        self.reference_lamp = ref_lamp
        self.evaluation_comment = eval_comment
        self.spectral_dict = self.set_spectral_features(header)
        # self.aperture = 1  # aperture number
        # self.beam = 1  # beam
        # self.dtype = self.dtype_dict[solution_type]  # data type
        # self.dispersion_start = 0  # dispersion at start
        # self.dispersion_delta = 0  # dispersion delta average
        # self.pixel_number = 0  # pixel number
        # self.doppler_factor = 0  # doppler factor
        # self.aperture_low = 0  # aperture low (pix)
        # self.aperture_high = 0  # aperture high
        # # funtions parameters
        # self.weight = 1
        # self.zeropoint = 0
        # self.ftype = self.ftype_dict[model_name]  # function type
        # self.forder = model_order  # function order
        # self.pmin = 0  # minimum pixel value
        # self.pmax = 0  # maximum pixel value
        # self.fpar = []  # function parameters

    @staticmethod
    def set_spectral_features(header):
        if header is None:
            log.error('Header has not been parsed')
        else:
            try:
                log.debug('Red Camera')
                dict = {'camera': 'red',
                        'grating': header['GRATING'],
                        'roi': header['ROI'],
                        'filter1': header['FILTER'],
                        'filter2': header['FILTER2'],
                        'slit': header['SLIT'],
                        'instconf': header['INSTCONF'],
                        'wavmode': header['WAVMODE'],
                        'cam_ang': header['CAM_ANG'],
                        'grt_ang': header['GRT_ANG']}
                # for key in dict.keys():
                    # print(key, dict[key])
                return dict
            except KeyError:
                log.debug('Blue Camera')
                dict = {'camera': 'blue',
                        'grating': header['GRATING'],
                        'ccdsum': header['CCDSUM'],
                        'filter1': header['FILTER'],
                        'filter2': header['FILTER2'],
                        'slit': header['SLIT'],
                        'serial_bin': header['PARAM18'],
                        'parallel_bin': header['PARAM22'],
                        'cam_ang': header['CAM_ANG'],
                        'grt_ang': header['GRT_ANG']}
                # for key in dict.keys():
                    # print(key, dict[key])
                return dict

    def check_compatibility(self, header=None):
        if header is not None:
            new_dict = self.set_spectral_features(header)
            for key in new_dict.keys():
                if self.spectral_dict['camera'] == 'red':
                    if key in ['grating', 'roi', 'instconf', 'wavmode'] and new_dict[key] != self.spectral_dict[key]:
                        log.debug('Keyword: %s does not Match', key.upper())
                        return False
                    elif key in ['cam_ang',  'grt_ang'] and abs(new_dict[key] - self.spectral_dict[key]) > 1:
                        log.debug('Keyword: %s Lamp: %s Data: %s',
                                  key,
                                  self.spectral_dict[key],
                                  new_dict[key])
                        return False
                    else:
                        return True
                elif self.spectral_dict['camera'] == 'blue':
                    if key in ['grating', 'ccdsum', 'serial_bin', 'parallel_bin']and new_dict[key] != self.spectral_dict[key]:
                        log.debug('Keyword: %s does not Match', key.upper())
                        return False
                    elif key in ['cam_ang',  'grt_ang'] and abs(float(new_dict[key]) - float(self.spectral_dict[key])) > 1:
                        log.debug('Keyword: %s Lamp: %s Data: %s',
                                  key,
                                  self.spectral_dict[key],
                                  new_dict[key])
                        return False
                    else:
                        return True
            return True
        else:
            log.error('Header has not been parsed')
            return False

    def linear_solution_string(self, header):
        pass


if __name__ == '__main__':
    log.error('This can not be run on its own.')
