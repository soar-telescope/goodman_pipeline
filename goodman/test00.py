import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import numpy as np
import sys
import scipy.interpolate
import pandas
import logging as log
import wsbuilder
from astropy.io import fits



class LineDetector(object):

    def __init__(self):
        self.args_reference_dir = '/user/simon/data/soar/work/eng/2016-10-14/RED/'
        # self.file_name = '/data/simon/data/soar/work/20161114_eng/reduced_data/fzh.0298_CVSO166_400m2_gg455.fits'
        self.file_name = '/user/simon/data/soar/work/goodman/test/etest2/gfzh.0047.SO2016A-019_0320.fits'
        self.data = fits.getdata(self.file_name)
        self.header = fits.getheader(self.file_name)
        self.region = None

    def __call__(self, *args, **kwargs):
        # lines_candidates = self.get_lines_in_lamp()
        # self.identify_spectra()
        solu = self.linear_solution()
        plt.plot(solu[0], solu[1])
        plt.show()

    def linear_solution(self):
        """Linear solution reader

        This method read the apropriate keywords and defines a linear wavelength solution

        Returns:
            solution (class): Callable wavelength solution model
        """
        crval = float(self.header['CRVAL1'])
        crpix = int(self.header['CRPIX1'])
        cdelt = float(self.header['CDELT1'])
        self.wcs_dict = {'crval': crval,
                         'crpix': crpix,
                         'cdelt': cdelt,
                         'dtype': 0}
        function = ReadMathFunctions(self.wcs_dict)
        solution = function.get_solution()
        # data = fits.getdata('/data/simon/data/soar/work/goodman/test/extraction-tests/CuHeAr_600.fits')
        x_axis = range(len(self.data))

        # plt.xlabel("%s (%s)" % (self.wat_wcs_dict['label'], self.wat_wcs_dict['units']))
        wave_intens = [solution(x_axis), self.data]
        # return solution
        print(crval)
        print(crpix)
        print(cdelt)
        return wave_intens

    def identify_spectra(self):
        # define data shape and samples

        y_size, x_size  = self.data.shape
        self.region = np.ones(y_size)
        half_width = int(0.03 * x_size)
        sample_loc = int(x_size / 2.)
        # print x_size, y_size, int(0.05 * x_size)
        sample = np.median(self.data[:, sample_loc:sample_loc + 2 * half_width], axis=1)
        sample_median = np.median(sample)
        sample_std = np.std(sample)

        # plt.figure()
        # plt.hist(sample, bins=100, log=True)
        # plt.show()

        # search for spectra
        all_found = []
        rising = False
        falling = False
        trend_length = 0
        for i in range(len(sample) - 1):
            if sample[i] > max(1.1 * sample_median, sample_std):
                if sample[i] > sample[i + 1]:
                    if rising and trend_length > 4:
                        all_found.append(i)
                        plt.axvline(i, color='r')
                        trend_length = 0
                    rising = False
                    falling = True
                    trend_length += 1
                else:
                    if falling and trend_length > 4:
                        plt.axvline(i, color='g')
                        trend_length = 0
                        pass
                    rising = True
                    falling = False
                    trend_length += 1

        # validate targets
        identified_targets = []
        for spectrum in all_found:
            # print spectrum
            for i in range(0, min(spectrum, abs(y_size - spectrum)) - 100):
                # print spectrum - i, spectrum + i, y_size - spectrum, x_size
                if sample[spectrum - i] < sample[(spectrum - i) - 1] or sample[spectrum + i] < sample[(spectrum + i) + 1]:
                    # plt.plot([x for x in range(spectrum - i, spectrum + i)], sample[spectrum - i:spectrum + i], color='g')
                    sample_width = abs(2 * i)

                    try:
                        sub_sample_x_axis = [x for x in range(spectrum - i, spectrum + i)]
                        sub_sample = sample[spectrum - i:spectrum + i]
                        # centroid = np.sum(sub_sample_x_axis * sub_sample) / np.sum(sub_sample)
                        # old procedure
                        voigt_init = models.Voigt1D(x_0=spectrum, amplitude_L=sample[spectrum], fwhm_L=8, fwhm_G=8)
                        fit_voigt = fitting.LevMarLSQFitter()
                        voigt = fit_voigt(voigt_init, sub_sample_x_axis, sub_sample)
                        print(voigt)
                        if abs(voigt.x_0.value - spectrum) < 5:
                            identified_targets.append(IdentifiedTarget(
                                voigt.amplitude_L.value,
                                voigt.x_0.value,
                                voigt.fwhm_L.value,
                                voigt.fwhm_G.value,
                                sample_width,
                                sample_loc))
                            # fitted_model = voigt
                            if True:
                                plt.axvspan(spectrum - i, spectrum + i, color='r', alpha=0.3, label='Spectrum width')
                                plt.plot(sub_sample_x_axis, voigt(sub_sample_x_axis), label='Voigt Fit')
                                # plt.axvline(centroid, color='k')
                                plt.axvline(voigt.x_0.value, color='m', label='Voigt x_0')
                                plt.axvline(spectrum, color='c', label='Spectrum location')
                        else:
                            log.info('Spectrum found at pixel %s is discarded', spectrum)
                    except TypeError:
                        pass
                    break

        if True:
            plt.axhline(sample_median, color='m', label='Median')
            plt.axhline(1.1 * sample_median, color='c', label='110% Median')
            plt.plot(sample, label='Data')
            plt.xlabel('Pixel (Spatial Direction)')
            plt.ylabel('Intensity')
            plt.legend(loc='best')
            plt.show()
        for target in identified_targets:
            print target
        return identified_targets

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
        lines_candidates = []
        prev_point = None
        status = None
        trend_length = 0
        median = np.median(self.data0)
        for pixel in range(len(self.data0)):
            point = self.data0[pixel]
            # print(point)
            if prev_point is None:
                prev_point = point
                status = 0
            else:
                if point > prev_point:
                    if status == 1:
                        trend_length += 1
                    elif status == -1:
                        trend_length = 1
                    status = 1
                elif point < prev_point:
                    if status == 1 and trend_length > 2 and point > median:
                        lines_candidates.append(pixel - 1)
                    status = -1
                    trend_length = 1
                else:
                    pass
            prev_point = point
        lines_center = self.recenter_lines(self.data0, lines_candidates)
        if True:
            for line in lines_candidates:
                plt.axvline(line, color='m')
            plt.axhline(median, color='g')
            plt.plot(self.data0, color='b')
            plt.show()
        return lines_center

    @staticmethod
    def recenter_lines(data, lines):
        new_center = []
        x_size = data.shape[0]
        median = np.median(data)
        for line in lines:
            #id left limit
            print line
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
                if True:
                    plt.axvspan(line - 1, line + 1, color='g', alpha=0.3)
                new_center.append(line)
            else:
                new_center.append(centroid)
            if True:
                # plt.axvline(centroid, color='m')
                plt.plot(sub_x_axis, sub_data)
                # plt.plot(data, color='r', linestyle='--')
                # plt.xlim([left_limit, right_limit])
                # plt.ylim([min(sub_data), max(sub_data)])
                # plt.show()
        if True:
            plt.axhline(median, color='b')
            plt.plot(data, color='r', linestyle='--', alpha=0.5)
            for line in lines:
                plt.axvline(line, color='m')
            for center in new_center:
                plt.axvline(center, color='c')
            plt.show()
        return new_center


class IdentifiedTarget(object):
    """Allows for easy storage and manipulation of the targets found.

    """

    def __init__(self, amplitude, mean, stddev, fwhmg=0, raw_width = None, sample_loc = None):
        self.amplitude = amplitude
        self.mean = mean
        self.stddev = stddev
        self.fwhm_g = fwhmg
        self.raw_width = raw_width
        self.sample_loc = sample_loc

class ReadMathFunctions(object):

    def __init__(self, wcs_dict):
        self.wcs = wcs_dict
        self.solution = None
        if self.wcs['dtype'] == -1:
            self.none()
        elif self.wcs['dtype'] == 0:
            self.solution = self.linear_solution()
        elif self.wcs['dtype'] == 1:
            self.log_linear()
        elif self.wcs['dtype'] == 2:
            if self.wcs['ftype'] == '1':
                self.solution = self.chebyshev()
            elif self.wcs['ftype'] == '2':
                self.non_linear_legendre()
            elif self.wcs['ftype'] == '3':
                self.non_linear_cspline()
            elif self.wcs['ftype'] == '4':
                self.non_linear_lspline()
            elif self.wcs['ftype'] == '5':
                # pixel coordinates
                pass
            elif self.wcs['ftype'] == '6':
                # sampled coordinate array
                pass
            else:
                log.error('Not Implemented')
        else:
            log.error('Not Implemented')

    @staticmethod
    def none():
        return 0

    def linear_solution(self):
        intercept = self.wcs['crval'] - (self.wcs['crpix'] - 1) * self.wcs['cdelt']
        # intercept = self.wcs['crval'] - self.wcs['crpix'] * self.wcs['cdelt']
        linear = models.Linear1D(slope=self.wcs['cdelt'], intercept=intercept)
        return linear

    @staticmethod
    def log_linear():
        return False

    def chebyshev(self):
        cheb = models.Chebyshev1D(degree=self.wcs['order'], domain=[self.wcs['pmin'], self.wcs['pmax']], )
        for param_index in range(self.wcs['order']):
            cheb.parameters[param_index] = self.wcs['fpar'][param_index]
        return cheb

    def non_linear_legendre(self):
        raise NotImplementedError

    def non_linear_lspline(self):
        raise NotImplementedError

    def non_linear_cspline(self):
        """Cubic Spline"""
        # cubic_spline = models.
        raise NotImplementedError

    def get_solution(self):
        if self.solution is not None:
            return self.solution
        else:
            log.error("The solution hasn't been found")


if __name__ == '__main__':
    lde = LineDetector()
    lde()
