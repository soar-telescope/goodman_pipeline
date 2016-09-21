from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from specutils.wcs import specwcs
import scipy.interpolate
import linelist
import logging as log
import glob # remove later
import sys


class WavelengthCalibration:
    def __init__(self, path, sci_pack, science_object):
        self.science_object = science_object
        self.slit_offset = 0
        self.interpolation_size = 200
        self.line_search_method = 'derivative'
        self.gratings_dict = {'SYZY_400': 400, 'KOSI_600': 600, '930': 930, '1200': 1200}
        self.grating_frequency = 0
        self.grating_angle = float(0)
        self.camera_angle = float(0)
        # self.binning = self.header0[]
        self.binning = 1
        self.pixel_count = 0
        self.alpha = 0
        self.beta = 0
        self.center_wavelength = 0
        self.blue_limit = 0
        self.red_limit = 0
        """this data must come parsed"""
        # self.path = '/data/simon/data/SOAR/work/goodman/test/extraction-tests-2/'
        # self.image_name = 'lamp-HgAr.fits'
        # self.path = '/data/simon/data/SOAR/work/goodman/test/extraction-tests/'
        self.path = path
        # self.image_name = 'eefc_0046.SO2016A-019_0320.fits'
        # self.image_name = 'eefc_0047.SO2016A-019_0320.fits'
        self.all_data = sci_pack[0]
        self.all_headers = sci_pack[1]
        self.sci = self.all_data[0]
        self.header = self.all_headers[0]
        log.info('Processing Science Target: %s' % self.header['OBJECT'])
        if self.science_object.lamp_count > 0:
            for l in range(1,self.science_object.lamp_count + 1):
                self.data0 = self.all_data[l]
                self.header0 = self.all_headers[l]
                log.info('Processing Comparison Lamp: %s' % self.header0['OBJECT'])
                self.data1 = self.interpolate(self.data0)
                self.lines_limits = self.get_line_limits()
                self.lines_center = self.get_line_centers(self.lines_limits)
                self.spectral = self.get_spectral_characteristics()
                self.wsolution = self.wavelength_solution()
            self.header = self.add_wavelength_solution(self.header)
        else:
            log.error('There are no lamps to process')

    def interpolate(self, spectrum):
        """Creates an interpolated version of the input spectrum

        This method creates an interpolated version of the input array, it is used mainly for a spectrum but it can
        also be used with any unidimensional array, assuming you are happy with the interpolation_size attribute
        defined for this class. The reason for doing interpolation is that allows to find the lines and its respective
        center more precisely. The default interpolation size is 200 (two hundred) points.

        Args:
            spectrum (array): an uncalibrated spectrum or any unidimensional array.

        Returns:
            Two dimensional array containing x-axis and interpolated array. The x-axis preserves original pixel values.

        """
        x_axis = range(1, len(spectrum) + 1)
        x0 = x_axis[0]
        x1 = x_axis[-1]
        new_x_axis = np.linspace(x0, x1, len(spectrum) * self.interpolation_size)

        tck = scipy.interpolate.splrep(x_axis, spectrum, s=0)
        new_spectrum = scipy.interpolate.splev(new_x_axis, tck, der=0)
        """
        plt.title('Interpolated Spectrum')
        plt.plot(new_x_axis, new_spectrum)
        plt.show()
        # """
        return [new_x_axis, new_spectrum]

    @staticmethod
    def recenter_line(x, y, max_val, max_pos, model='gauss'):
        """Recalculates the center of a line by fitting a model

        This method fits one of the three line profiles to the input data. It works well in lines with good signal to
        noise ratios and with little contaminants (please see Notes). This method will be most likely removed since
        the combination of the line detection method plus the interpolation of the spectrum gives a very good result.

        Notes:
            This method was kept for development purposes. According to the test performed the line centers that were
            fed to this method were always better than the output since the fits were confused by contaminants.

        Args:
            x (array): x-axis of data in pixel values.
            y (array): the data to be fitted. A section around an spectral line.
            max_val (float): Value at line's peak.
            max_pos (int): Index of line's peak.
            model (str): Name of the line profile to be fitted. One of this three: "gauss", "lorentz" or "voigt".

        Returns:
            new_center (float): Value of the re-calculated center of the line.

        """
        if len(x) == len(y):
            if model == 'gauss':
                gauss_init = models.Gaussian1D(amplitude=max_val, mean=max_pos, stddev=4)
                fit_gaussian = fitting.LevMarLSQFitter()
                gauss = fit_gaussian(gauss_init, x, y)
                fitted_model = gauss
                new_center = gauss.mean.value
            elif model == 'lorentz':
                lorentz_init = models.Lorentz1D(amplitude=max_val, x_0=max_pos, fwhm=8)
                # log.debug('Amplitude: %s X_0: %s FWHM: %s', max_val, max_pos, 8)
                fit_lorentz = fitting.LevMarLSQFitter()
                lorentz = fit_lorentz(lorentz_init, x, y)
                # print(lorentz)
                fitted_model = lorentz
                new_center = lorentz.x_0.value
            elif model == 'voigt':
                voigt_init = models.Voigt1D(x_0=max_pos, amplitude_L=max_val, fwhm_L=8, fwhm_G=8)
                fit_voigt = fitting.LevMarLSQFitter()
                voigt = fit_voigt(voigt_init, x, y)
                print(voigt)
                fitted_model = voigt
                new_center = voigt.x_0.value
            # return new_center
            return new_center
        else:
            return False

    def get_line_limits(self):
        """Method for identifying lines in a spectrum

        This is the spectral line identifying method. It calculates a pseudo-derivative of the spectrum thus determining
        peaks. For a typical spectrum a pseudo-derivative, from now on just "derivative", will produce a series of
        positive and negative peaks, for emission lines. Since we will be using it for comparison lamps we don't have
        to worry about absorption lines and therefore this method only detects emission lines.
        A threshold is defined by calculating the 75 percent of the standard deviation of the derivative. There
        is no particular reason for choosing 75 percent is just the value that worked best for all the test subjects.
        Then any search for lines must be done for values above and below of the threshold and its negative
        respectively.  There are a few control mechanisms that ensures that an actual good line is being detected. For
        instance: For each line both positive and negative derivative's peaks are stored and are called limits. In order
        to add a (first) positive limits the number of previously stored limits must be even, meaning that we are
        "opening" a new line. Correspondingly if we are adding  a (last) negative limit we have to have added previously
        a odd amount of limits meaning that we are "closing" a line. At this point there is another constraint, it
        cannot be separated by an amount larger than a defined spacing, or we would be having a very broad line and lamp
        lines are expected to be very narrow.

        Notes:
            The lines detected by this method are usually the best lines in the spectrum. There is a certain amount of
            lines missed but in very crowded regions.

            Very important to note that the line centers found using the data produced here are very precise.

        Returns:
            limits (list): List of line limits in consecutive pairs. This is considered by the method that uses this
            list for further processing. The values are index values not pixel values.

        """
        """
        pixel_axis = self.data1[0]
        if self.line_search_method == 'threshold':
            median_value = np.median(self.data1[1])
            mean_value = np.mean(self.data1[1])
            threshold = median_value + mean_value
            keep_searching = True
            features = []
            while keep_searching:
                feature = []
                subx = []
                for i in range(len(self.data1[1])):
                    if self.data1[1][i] > threshold:
                        feature.append(self.data1[1][i])
                        subx.append(pixel_axis[i])
                    elif feature != [] and len(feature) >= 3:
                        features.append([subx, feature])
                        # print(len(feature))
                        feature = []
                        subx = []
                    else:
                        feature = []
                        subx = []
                    if i == len(self.data1[1]) - 1:
                        keep_searching = False
            print(len(features))
            return False
        """
        if self.line_search_method == 'derivative':
            """in fact is a pseudo-derivative"""
            derivative = []
            # derivative2 = []
            faux_x = range(0, len(self.data1[1]) - 1)
            # faux_x2 = range(0, len(self.data1[1]) - 2)
            # print faux_x
            for i in faux_x:
                derivative.append(self.data1[1][i + 1] - self.data1[1][i])
            # for e in faux_x2:
            #    derivative2.append(derivative[e] - derivative[e+1])
            threshold = np.std(derivative) * .75
            new_range = 0
            spacing = 1500
            limits = []
            for i in range(len(derivative) - 1):
                if i > new_range:
                    if derivative[i] > threshold and derivative[i + 1] - derivative[i] >= 0:
                        partial_max = i + np.argmax(derivative[i:i + spacing])
                        # print i, partial_min
                        new_range = partial_max + (partial_max - i)
                        if limits == [] or len(limits) % 2 == 0:
                            limits.append(partial_max)
                        else:
                            plt.axvline(partial_max, color='k')
                    elif derivative[i] <  -threshold and derivative[i + 1] - derivative[i] <= 0:
                        partial_min = i + np.argmin(derivative[i:i + spacing])
                        new_range = partial_min + (partial_min - i)
                        if len(limits) % 2 == 1 and partial_min - limits[-1] < spacing:
                            limits.append(partial_min)
                            # plt.axvline(partial_max, color='g')
                        elif limits != []:
                            if partial_min - limits[-1] > spacing:
                                plt.axvline(partial_min, color='m')
                                limits = limits[:-1]
            if len(limits) % 2 == 1:
                limits = limits[:-1]
            """
            # Produce Plots
            for i in range(len(limits)):
                if i % 2 == 0:
                    plt.axvline(limits[i], color='r')
                elif i % 2 == 1:
                    plt.axvline(limits[i], color='g')

            plt.title('Line Identification')
            # plt.plot(self.data1[1], label='Spectrum')
            plt.plot(faux_x, derivative, label='1st Derivative')
            plt.axhline(0, color='m')
            # plt.plot(faux_x2, derivative2, label='2nd')
            plt.axhline(threshold)
            plt.axhline(-threshold)
            plt.legend(loc='best')
            # plt.plot(pixel_axis, self.data1[1])
            plt.savefig(self.path + 'line-identification-201.png', dpi=300)
            plt.show()
            # """
            return limits

    def get_line_centers(self, limits):
        """Finds the center of the lines using limits previously found

        This method is very simple and could be integrated in the get_line_limits method but I'd rather have the limit
        information available for later. Basically finds the mean value of the line limits and then finds the
        correspoing pixel value, adds it up to the "centers" list.

        Args:
            limits (list): Line limits in the list's index domain.

        Returns:
            centers (list): Line centers in pixel values as floats.

        """
        centers = []
        for i in range(0, len(limits), 2):
            center = (limits[i] + limits[i + 1]) / 2.
            width = limits[i + 1] - limits[i]
            pixel_width = self.data1[0][limits[i + 1]] - self.data1[0][limits[i]]
            log.debug('Approximate FWHM: %s pix %s Angstrom (pix * 0.65)', pixel_width, pixel_width * 0.65)
            i_min = int(center - 2 * width)
            i_max = int(center + 2 * width)
            pixel_axis = self.data1[0][i_min:i_max]
            data_axis = self.data1[1][i_min:i_max]
            pixel_center = self.data1[0][int(round(center))]
            center_val = self.data1[1][int(round(center))]
            new_center = self.recenter_line(pixel_axis, data_axis, center_val, pixel_center, 'gauss')
            """
            plt.plot(pixel_axis, data_axis)
            plt.axvline(pixel_center)
            plt.axvline(new_center, color='m')
            plt.show()
            """
            centers.append(pixel_center)
            # print(center, width)
        return centers

    def get_spectral_characteristics(self):
        """Calculates some Goodman's specific spectroscopic values.

        Returns:
            spectral_characteristics (dict):

        """
        self.grating_frequency = self.gratings_dict[self.header0['GRATING']]
        self.grating_angle = float(self.header0['GRT_ANG'])
        self.camera_angle = float(self.header0['CAM_ANG'])
        # binning = self.header0[]
        # TODO GET BINNING FROM THE RIGHT SOURCE
        self.binning = 1
        # self.pixel_count = len(self.data0)
        """Calculations"""
        self.alpha = self.grating_angle + self.slit_offset
        self.beta = self.camera_angle - self.grating_angle
        self.center_wavelength = 10 * (1e6 / self.grating_frequency) * (
        np.sin(self.alpha * np.pi / 180.) + np.sin(self.beta * np.pi / 180.))
        self.blue_limit = 10 * (1e6 / self.grating_frequency) * (
        np.sin(self.alpha * np.pi / 180.) + np.sin((self.beta - 4.656) * np.pi / 180.))
        self.red_limit = 10 * (1e6 / self.grating_frequency) * (
        np.sin(self.alpha * np.pi / 180.) + np.sin((self.beta + 4.656) * np.pi / 180.))
        pixel_one = self.predicted_wavelength(1)
        pixel_two = self.predicted_wavelength(2)
        # print pixel_one, pixel_two, pixel_two-pixel_one
        spectral_characteristics = {'center': self.center_wavelength,
                                    'blue': self.blue_limit,
                                    'red': self.red_limit,
                                    'alpha': self.alpha,
                                    'beta': self.beta,
                                    'pix1': pixel_one,
                                    'pix2': pixel_two}
        print spectral_characteristics
        return spectral_characteristics

    def predicted_wavelength(self, pixel):
        a = self.alpha
        b = self.beta
        # c = self.pixel_count
        d = self.binning
        e = self.grating_frequency
        wavelength = 10 * (1e6 / e) * (np.sin(a * np.pi / 180.) +
                                  np.sin((b * np.pi / 180.) +
                                         np.arctan((pixel * d - 2048) * 0.015 / 377.2)))
        return wavelength

    def wavelength_solution(self):
        """Preliminary Solution"""
        pixel_axis = range(1, 4097)
        # pixel_axis = pixel_axis - pixel_axis.size // 2
        slope = self.spectral['pix2'] - self.spectral['pix1']
        intercept = self.spectral['pix1'] - slope
        first_solution = models.Linear1D(slope, intercept)

        fit_linear = fitting.LinearLSQFitter()
        """Pixel to wavelength parameters"""
        # manual input for testing
        # pix0 = [205.39, 805.99, 1280.07, 2967.24, 3445.09, 3477.87]
        # wav0 = [3650.153, 4046.563, 4358.328, 5460.735, 5769.598, 5790.663]
        """
        pix1 = [1563.51,
                1617.94,
                1633.06,
                1645.16,
                1663.31,
                1853.83,
                1902.22,
                1962.7,
                2026.21,
                3021.17,
                3042.34,
                3117.94,
                3139.11,
                3195.55,
                3257.06]
        # """
        pix1 = [1564.375,
                1616.652,
                1632.724,
                1662.721,
                1841.678,
                1855.135,
                1899.701,
                1962.554,
                2026.459,
                3021.701,
                3037.413,
                3118.79,
                3139.954,
                3192.895,
                3260.514]
        wav1 = [4545.0519,
                4579.3495,
                4589.8978,
                4609.5673,
                4726.8683,
                4735.9058,
                4764.8646,
                4806.0205,
                4847.8095,
                5495.8738,
                5506.1128,
                5558.702,
                5572.5413,
                5606.733,
                5650.7043]

        cuhear = [[3888.646, 'HeI'],
                  [3948.979, 'ArI'],
                  [3964.727, 'HeI'],
                  [4026.189, 'HeI'],
                  [4044.418, 'ArI'],
                  [4158.591, 'ArI'],
                  [4181.884, 'ArI'],
                  [4200.675, 'ArI'],
                  [4259.362, 'ArI'],
                  [4333.561, 'ArI'],
                  [4471.477, 'HeI'],
                  [4545.052, 'ArII'],
                  [4657.901, 'ArII'],
                  [4713.143, 'HeI'],
                  [4764.865, 'ArII'],
                  [4806.021, 'ArII'],
                  [4879.864, 'ArII'],
                  [4921.929, 'HeI'],
                  [4965.080, 'ArII'],
                  [5015.675, 'HeI'],
                  [5162.285, 'ArI'],
                  [5187.746, 'ArI'],
                  [5495.874, 'ArI'],
                  [5739.520, 'ArI'],
                  [5875.618, 'HeI'],
                  [5912.085, 'ArI'],
                  [6416.307, 'ArI'],
                  [6752.834, 'ArI'],
                  [6871.289, 'ArI'],
                  [6937.664, 'ArI'],
                  [6965.431, 'ArI'],
                  [7147.042, 'ArI'],
                  [7206.980, 'ArI'],
                  [7272.936, 'ArI'],
                  [7281.349, 'HeI'],
                  [7353.293, 'ArI'],
                  [7372.118, 'ArI'],
                  [7383.981, 'ArI'],
                  [7503.869, 'ArI'],
                  [7514.652, 'ArI'],
                  [7635.106, 'ArI'],
                  [7948.176, 'ArI'],
                  [8006.157, 'ArI'],
                  [8014.786, 'ArI'],
                  [8103.693, 'ArI'],
                  [8115.311, 'ArI'],
                  [8264.523, 'ArI'],
                  [8408.210, 'ArI'],
                  [8424.648, 'ArI'],
                  [8521.442, 'ArI'],
                  [8667.944, 'ArI'],
                  [8799.088, 'ArI']]

        cuhear_colors = {'HeI': 'c', 'ArI': 'r', 'ArII': 'm',}


        lines_in_range = self.get_lines_in_range(self.spectral['blue'],
                                                 self.spectral['red'],
                                                 self.header0['OBJECT'])
        fos = fit_linear(first_solution, pix1, wav1)
        print fos
        inverse_solution = models.Linear1D(1/fos.slope.value, - fos.intercept.value / fos.slope.value)

        lines_to_pix = []
        for rline in lines_in_range:
            lines_to_pix.append(inverse_solution(rline))
        pixel_offset = self.pixel_axis_cross_correlate(self.lines_center, lines_in_range, lines_to_pix)
        new_lines_estimated = []
        differences = []
        for pixel in self.lines_center:
            new_pixel = pixel + pixel_offset
            new_first_solution = first_solution(new_pixel)
            dif = abs(np.array(lines_in_range) - new_first_solution)
            i_min = np.argmin(dif)
            differences.append(dif[i_min])
            # print(dif[i_min],new_first_solution,lines_in_range[i_min])
            new_lines_estimated.append(new_first_solution)
        wavel_center = []
        pixel_center = []
        for i in range(len(new_lines_estimated)):
            if np.min(abs(np.array(lines_in_range) - new_lines_estimated[i])) > 20:
                print('Line Rejected ', new_lines_estimated[i], self.lines_center[i])
            else:
                arg_min = np.argmin(abs(np.array(lines_in_range) - new_lines_estimated[i]))
                ref_value = lines_in_range[arg_min]
                wavel_center.append(ref_value)
                pixel_center.append(self.lines_center[i])
        slope_array = []
        for e in range(0,len(wavel_center) - 1):
            new_slope = (wavel_center[e + 1] - wavel_center[e])/(pixel_center[e + 1] - pixel_center[e])
            slope_array.append(new_slope)
            print('Slope : ', new_slope)
        if slope_array != []:
            first_solution.slope = np.median(slope_array)
            print('Median New Slope : ', np.median(slope_array))

        # self.lines_center.pop(4)
        # new_lines_estimated.pop(4)
        print pixel_center
        print wavel_center

        """experiment"""
        forced_sol = models.Linear1D(slope=0.65, intercept=first_solution.intercept.value)

        for line in lines_in_range:
            plt.axhline(line, color='r')
        for line in self.lines_center:
            plt.axvline(line, color='g')
        plt.title('f(x) = %s x + %s' % (first_solution.slope.value, first_solution.intercept.value))
        plt.xlabel('Pixels')
        plt.ylabel('Angstrom')
        plt.plot(pixel_axis, first_solution(pixel_axis), color='b')
        plt.plot(pixel_axis, forced_sol(pixel_axis), color='m')
        plt.plot(pix1, wav1, marker='o', color='c')


        first_solution = models.Linear1D(slope, intercept)
        fit_wavelength = fitting.LinearLSQFitter()
        second_solution = fit_wavelength(first_solution, pix1, wav1)
        plt.plot(pixel_axis, second_solution(pixel_axis), linestyle='--', color='k')
        plt.show()

        reference_file = '/data/simon/data/SOAR/work/goodman/test/extraction-tests/CuHeAr_600.fits'
        reference_data = fits.getdata(reference_file)
        reference_header = fits.getheader(reference_file)
        # reference_value = float(reference_header['CRVAL1'])
        # reference_pixel = int(reference_header['CRPIX1'])
        # delta = float(reference_header['CDELT1'])
        # text = reference_header['WAT1_001'].split()
        # text_dict = {}
        # for value in text:
        #    key, val = value.split('=')
        #     text_dict[key] = val
        # length = len(reference_data)

        # start = reference_value - (reference_pixel - 1) * delta
        # stop = start + (length - 1) * delta
        # print 'Start ', start, 'Stop ', stop, 'Length ', length
        # wavelength_axis = np.linspace(start, stop, length)
        wavelength_axis = self.get_wavelength_solution(reference_header)
        rf = '/data/simon/data/SOAR/work/goodman/test/extraction-tests/CuHeAr_600_nonlinear.fits'
        rd = fits.getdata(rf)
        rh = fits.getheader(rf)
        rwa = self.get_wavelength_solution(rh)




        plt.figure(1)
        plt.subplot(211)
        plt.xlim((0, 4096))
        plt.plot(pixel_axis, self.data0)
        for line in self.lines_center:
            plt.axvline(line, color='g')
        plt.xlabel('Pixel')



        """This solution works"""
        poli_init = models.Polynomial1D(degree=2) # , domain=[0, len(self.data0)])
        poli_init.c0 = 1
        poli_init.c1 = 1
        poli_init.c2 = 1
        fit_poli = fitting.LinearLSQFitter()
        poli_solution = fit_poli(poli_init, pixel_center, wavel_center)

        """--------------"""
        plt.subplot(212)
        # plt.xlim((start, stop))

        for culine in cuhear:
            plt.axvline(culine[0], color=cuhear_colors[culine[1]])
        # for w in wav1:
            # plt.axvline(w, color='c')
        # for line in lines_in_range:
            # plt.axvline(line, color='r')
        # plt.plot(poli_solution(pixel_axis), self.data0)
        plt.plot(wavelength_axis, reference_data, color='g')
        # plt.plot(rwa, rd, color='b')
        # for p in linelist.line_list['CuAr']:
            # plt.axvline(p, color='k', linestyle='-.')
        plt.xlabel('Angstrom - poli')
        plt.tight_layout()
        plt.show()
        """--------------"""


        plt.plot(pixel_axis, poli_solution(pixel_axis), color='m', label='Poli')

        # plt.plot(pixel_axis, solution(pixel_axis), label='Solution')

        plt.plot(pixel_center, wavel_center, color='g', linestyle='--', marker='o',label='New Last')
        plt.plot(pix1, wav1, color='r', linestyle='--', marker='o', label='Real')
        plt.plot(pixel_axis, fos(pixel_axis), label='First Order')
        plt.legend(loc='best')
        plt.show()

        """Continue work from here...."""

        pixel_to_angstrom_ii = poli_solution(self.lines_center)

        for line in pixel_to_angstrom_ii:
            plt.axvline(line, color='c', linestyle='-.')
        """
        for nline in wavel_center:
            plt.axvline(nline, color='r')
            print(nline)
        for rline in lines_in_range:
            plt.axvline(rline, linestyle='--', color='g')
        # """
        plt.plot(second_solution(pixel_axis), self.data0, color='r')
        plt.plot(fos(pixel_axis), self.data0, color='c')

        plt.plot(poli_solution(pixel_axis), self.data0)
        plt.xlabel('Wavelength (angstrom)')
        plt.ylabel('Intensity (counts)')
        plt.legend()
        plt.show()




        plt.plot(pixel_axis, self.data0)
        for pixel in self.lines_center:
            plt.axvline(pixel, color='r')
        plt.show()
        print first_solution
        print second_solution
        print poli_solution
        # """
        return poli_solution

    @staticmethod
    def get_lines_in_range(blue, red, lamp):
        lines = []
        if len(lamp) % 2 == 0:
            for e in range(0,len(lamp),2):
                element = lamp[e:e+2]
                all_lines = linelist.line_list[element]
                for line in all_lines:
                    if blue <= line <= red:
                        lines.append(line)
                        # print(" YES %s <= %s <= %s" % (blue, line, red))
                    # else:
                        # print(" NO %s <= %s <= %s" % (blue, line, red))
        # print lines
        return sorted(lines)

    def pixel_axis_cross_correlate(self, reference, lines_in_range, pixel_lines):

        root_pixel_axis = np.linspace(0, 4096, 4096)
        reference_axis = np.zeros(len(root_pixel_axis))
        pixel_lines_axis = np.zeros(len(root_pixel_axis))
        for ref_line in reference:
            gaussian = models.Gaussian1D(amplitude=1, mean=ref_line, stddev=5)
            reference_axis += gaussian(root_pixel_axis)
        for pix_line in pixel_lines:
            gaussian = models.Gaussian1D(amplitude=1, mean=pix_line, stddev=5)
            pixel_lines_axis += gaussian(root_pixel_axis)

        """Cross correlate"""
        lag_position = range(-400, 400, 1)
        correlation = []
        for lag in lag_position:
            correlation_value = 0
            for i in range(len(reference_axis)):
                i_ref = i
                i_com = i + lag
                if 0 < i_com < len(reference_axis):
                    correlation_value += reference_axis[i_ref] * pixel_lines_axis[i_com]
            correlation.append(correlation_value)
        i_max = np.argmax(correlation)

        # """
        plt.title('Lag Position: %s' % lag_position[i_max])
        plt.axvline(lag_position[i_max])
        plt.plot(lag_position, correlation)
        plt.show()


        plt.plot(root_pixel_axis, reference_axis, color='g', label='Reference')
        plt.plot(root_pixel_axis, pixel_lines_axis, color='r', label='Lines')
        plt.legend(loc='best')
        plt.show()
        # """
        return lag_position[i_max]


    def add_wavelength_solution(self, new_header):

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
        new_header['DCLOG1'] = 'REFSPEC1 = %s' % c_lamp

        # fits.writeto(out_lamp, new_data, new_header, clobber=True)
        return new_header

    @staticmethod
    def get_wavelength_solution(header):
        """Reproduces wavelength solution from the image's header



        Args:
            header:

        Returns:

        """
        ctype1 = header['CTYPE1']
        if ctype1 == 'LINEAR':
            reference_value = float(header['CRVAL1'])
            reference_pixel = int(header['CRPIX1'])
            delta = float(header['CDELT1'])
            text = header['WAT1_001'].split()
            text_dict = {}
            for value in text:
                key, val = value.split('=')
                text_dict[key] = val
            if int(header['NAXIS']) == 1:
                length = int(header['NAXIS1'])
                start = reference_value - (reference_pixel - 1) * delta
                stop = start + (length - 1) * delta
                log.debug('Start %s Stop %s Length %s', start, stop, length)
                wavelength_axis = np.linspace(start, stop, length)
                return wavelength_axis
            else:
                log.error('Can not work with multi-axis files.')
                return False
        elif ctype1 == 'MULTISPE':
            if int(header['WCSDIM']) == 2:
                try:
                    ctype2 = header['CTYPE2']
                    if ctype2 == 'MULTISPE':
                        cdelt2 = header['CDELT2']
                        cd22 = header['CD2_2']
                        ltm22 = header['LTM2_2']
                        waxmap01 = header['WAXMAP01']
                        wat2_keys = header['WAT2_*']
                        wat2 = ''
                        for key in wat2_keys:
                            wat2 += header[key]

                        wat = wat2.split('"')
                        """
                        taken from ftp://iraf.noao.edu/iraf/web/projects/fitswcs/specwcs.html
                        and http://shaileshahuja.blogspot.cl/2014/06/analysing-iraf-multispec-format-fits.html
                        wat index and meaning
                        0 : Aperture
                        1 : Beam
                        2 : dtype
                            -1 Spectrum not calibrated
                            0 Linear dispersion sampling
                            1 Log-linear dispersion
                            2 Nonlinear dispersion.
                        3 : Dispersion Value at start
                        4 : Average Dispersion delta
                        5 : Number of Pixels
                        6 : Doppler Factor
                        7 : Aperture Low
                        8 : Aperture High
                        ---
                        9 : Weight
                        10 : Zero-point offset
                        11 : Function type
                            1 Chebyshev polynomial
                            2 Legendre polynomial
                            3 Cubic spline
                            4 Linear spline
                            5 Pixel coordinate array
                            6 Sampled coordinate array
                        12 : Order of function
                        13 : Minimum Pixel Value
                        14 : Maximum Pixel Value
                        15+ : Coefficients of functions i

                        """

                        if 'spec' not in wat[1]:
                            print(True)

                        print(wat2)
                        print(wat[1])
                        aperture = wat[1][0]
                        beam = wat[1][1]
                        dtype = wat[1][2]
                        
                    else:
                        log.debug('Nothing to do!')
                except KeyError:
                    log.error('KeyError')
        elif ctype1 == 'PIXEL':
            log.error('This header does not contain a wavelength solution')
            return False






"""
# path = '/data/simon/data/SOAR/work/goodman/test/extraction-tests/'
# image_list = ['eefc_0046.SO2016A-019_0320.fits']
for image in image_list:
    data0 = fits.getdata(path + image)
    header1 = fits.getheader(path + image)

    data1 = interpolate(data0)
    x_data = interpolate(np.linspace(0,len(data0),len(data0)))

    median = np.median(data1)
    mean = np.mean(data1)
    threshold = median + mean
    keepsearching = True
    features = []
    while keepsearching:
        feature = []
        subx = []
        for i in range(len(data1)):
            if data1[i] > threshold:
                feature.append(data1[i])
                subx.append(x_data[i])
            elif feature != [] and len(feature) >= 3:
                features.append([subx, feature])
                print(len(feature))
                feature = []
                subx = []
            else:
                feature = []
                subx = []
            if i == len(data1) - 1:
                keepsearching = False

    print(len(features))

    # crval1 = float(header1['CRVAL1'])
    # crpix1 = int(header1['CRPIX1'])
    # delta1 = float(header1['CDELT1'])

    # start1 = crval1 - (crpix1 - 1) * delta1
    # stop1 = start1 + (len(data1) - 1) * delta1
    # xa1 = np.linspace(start1, stop1, len(data1))


    # plt.plot(xa1, data1, label='New Calibrated')
    # plt.axhline(median)
    print('Please Write down feature number and wavelength value, it will be requested later')
    for i in range(len(features)):
        line = features[i]
        plt.plot(line[0], line[1])
        imax = np.max(line[1])
        amax = line[0][np.argmax(line[1])]
        plt.annotate(str(i), xy=(amax, imax + 5), xycoords='data')
    plt.title(header1['GRATING'])
    plt.axhline(threshold, color='r')
    plt.legend(loc='best')
    # plt.savefig(path + 'features-id.png', dpi=300)
    plt.show()
    plt.clf()
    print('Please enter wavelength value for feature')
    print('Enter 0 to ignore')
    line_list = {0: 3650.153,
                 1: 0,
                 2: 4046.563,
                 3: 0,
                 4: 0,
                 5: 4358.328,
                 6: 5460.735,
                 7: 0,
                 8: 5769.598,
                 9: 5790.663,
                 10: 0,
                 11: 0,
                 12: 0,
                 13: 0,
                 14: 0,
                 15: 0,
                 16: 0,
                 17: 0,
                 18: 0}
    input_data = []
    for i in range(len(features)):
        line = features[i]
        imax = line[0][np.argmax(line[1])]
        # wavelength = raw_input('Wavelength value for feature %s : (%s)' % (i, line_list[i]))
        wavelength = raw_input('Wavelength value for feature %s : ' % (i)) or '0'
        # if line_list[i] != 0:
        if wavelength != '0':
            max_val = np.max(line[1])
            gauss = recenter_line(line[0], line[1], max_val, imax, model='gauss')
            # nmax = voigt.x_0.value
            nmax = gauss.mean.value
            input_data.append([nmax, float(wavelength)])
            plt.title('Feature %s'%i)
            plt.plot(line[0],line[1], label='Data')
            plt.plot(line[0],gauss(line[0]), label='Gauss')
            plt.axvline(nmax, color='r', label='New Center')
            plt.axvline(imax, label='Old Center')
            plt.legend(loc='best')
            filename = path + 'recentering-task-feature-%s.png'%i
            # plt.savefig(filename,dpi=300)
            # plt.show()
            plt.clf()
    print input_data
    p = []
    w = []
    for point in input_data:
        p.append(point[0])
        w.append(point[1])
    wavelength_fit = 'linear'
    if wavelength_fit == 'linear':
        linear_init = models.Linear1D(1, 1)
        fit_wavl = fitting.LinearLSQFitter()
        wav = fit_wavl(linear_init, p, w)
    elif wavelength_fit == 'poli':
        poli_init = models.Polynomial1D(degree=2, domain=[p[0], p[-1]])
        poli_init.c0 = 1
        poli_init.c1 = 1
        poli_init.c2 = 1
        fit_poli = fitting.LinearLSQFitter()
        wav = fit_poli(poli_init, p, w)

    print(wav)
    crpix = 1
    crval = wav(crpix)
    cdelt = wav(crpix + 1) - wav(crpix)
    print('CRPIX1 %s CRVAL1 %s CDELT1 %s' % (crpix, crval, cdelt))

    header1['CRVAL1'] = crval
    header1['CRPIX1'] = crpix
    header1['CDELT1'] = cdelt
    header1['CD1_1'] = cdelt
    header1['CTYPE1'] = 'LINEAR'
    header1['WAT1_001'] = 'wtype=linear label=Wavelength units=angstroms'

    fits.writeto(path + 'python-calibrated.fits', data0, header1, clobber=True)
    print(wav(1), wav(2) - wav(1))
    plt.plot(p, w, marker='o', linestyle='--', label='points')
    plt.plot(p, wav(p), label='Model')
    plt.legend(loc='best')
    plt.show()
"""

if __name__ == '__main__':
    wav_cal = WavelengthCalibration()
