from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting
import logging as log

log.basicConfig(level=log.DEBUG)


# import sys


class Process:
    """hrllo

    hello

    """

    def __init__(self, source_path, sci_obj):
        self.science_object = sci_obj
        self.path = source_path
        self.data = fits.getdata(self.path + self.science_object.file_name)
        self.header = fits.getheader(self.path + self.science_object.file_name)
        self.lamps_data = []
        self.lamps_header = []
        if self.science_object.lamp_count > 0:
            for lamp_index in range(self.science_object.lamp_count):
                lamp_data = fits.getdata(self.path + self.science_object.lamp_file[lamp_index])
                # lamp_type = self.science_object.lamp_type[lamp_index]
                self.lamps_data.append(lamp_data)
                lamp_header = fits.getheader(self.path + self.science_object.lamp_file[lamp_index])
                self.lamps_header.append(lamp_header)
        self.targets = self.identify_spectra()
        self.traces = self.trace(self.targets)
        self.extracted_data = self.extract(self.traces)
        self.wavelength_calibration(self.extracted_data)

    def identify_spectra(self):
        """Finds the location of the spectrum or spectra in case there are more

        Performs an analysis along the spatial direction to determine the presence of important peaks which will be
        identified as science targets (or objects). A sample of the image of 200 pixels is taken more or less in the
        middle of the spectrum and then flattened along the dispersion direction, thus having enough signal to identify
        and differentiate a target or the targets from background signal. The method to find the spectrum or spectra
        is by defining a threshold which will be the Median + Standard Deviation from this unidimentional profile
        (spatial direction). For every group of points that are above the threshold a python list will be created and
        then its maximum value's position is found. Then for every position it tries to fit a gaussian profile in the
        whole subsample (not just the list with the points above the threshold). [THIS PART IS NOT IMPLEMENTED] If the
        gaussian fit yields something rare the candidate is discarded.

        I doesn't take arguments because it will read the <object> **self.science_object** that contains all the
        necessary data.

        Returns:
            A list of IdentifiedTargets class objects

        """

        x, y = self.data.shape
        # print(x,y)
        sample_data = np.median(self.data[:, int(y / 2.) - 100:int(y / 2.) + 100], axis=1)
        x_axis = np.linspace(0, len(sample_data), x)
        std = np.std(sample_data)
        median = np.median(sample_data)
        threshold = median + std
        keep_searching = True
        all_candidates = []
        while keep_searching:
            candidate = []
            candidate_index = []
            found = False
            for i in range(len(sample_data)):
                if sample_data[i] > threshold:
                    found = True
                    candidate.append(sample_data[i])
                    candidate_index.append(i)
                elif found:
                    all_candidates.append([candidate, candidate_index])
                    candidate = []
                    candidate_index = []
                    found = False
                if i == len(sample_data) - 1:
                    keep_searching = False
        identified_targets = []
        for single_candidate in all_candidates:
            cmax = np.argmax(single_candidate[0])
            max_val = np.max(single_candidate[0]) / 2.
            max_pos = single_candidate[1][cmax]
            # print("max-val", max_val, "max-pos", max_pos)
            gauss_init = models.Gaussian1D(amplitude=max_val, mean=max_pos, stddev=1)
            fit_gaussian = fitting.LevMarLSQFitter()
            gauss = fit_gaussian(gauss_init, x_axis, sample_data)
            identified_targets.append(IdentifiedTarget(gauss.amplitude.value, gauss.mean.value, gauss.stddev.value))
        return identified_targets

    def trace(self, targets):
        """Finds the trace of a spectrum given an initial location

        From the input data it finds the location of the spectrum in the image and an approximate width of six sigmas,
        with that information it takes a sample of the data just covering those limits to speed up the process.
        It takes about fifty sub-samples along the dispersion direction with a width of fifty pixels each.
        All the parameters can be changed with the variables half_n_sigma, ends_pix_spacing, n_samples and
        sample_width variables. The sub-samples are flattened along the dispersion direction using numpy.median and
        the location of the maximum value is recorded. Alternatively is possible to do a gaussian fit but my tests
        suggest is not much what you earn with this now is disabled by default. Enable it by changing do_gaussian_fit to
        True. Once there is a list of maximum locations a Chebyshev 1D of second order is fitted to define the trace.
        The Chebyshev function is defined for all the range of pixels in the dispersion direction of the image.

        Notes:
            Is not worth doing gaussian fits to the sub samples since it will not improve the result. The only case this
            might be useful is the case when someone does an interpolation along the spatial direction.

        Args:
            targets (list): Each element is a class that stores the parameters of the gaussian fitted to the data in
            previous steps. This data tells the location in the image of the target or targets.

        Returns:
            traces (list): Every element is a list with two elements. The first is the fitted Chebyshev class and the
            second and last is the width to be extracted.

        """

        """True enables gaussian fits and false disables it. Gaussian fit was left for experimental purposes"""
        do_gaussian_fit = False
        "half of number of sigmas to be sub_sampled"
        half_n_sigma = 3
        """space allowance in pixels for extreme of images"""
        ends_pix_spacing = 50
        """Number of samples to be done in the dispersion direction"""
        n_samples = 50
        """sample width must be smaller or equal than ends_pix_spacing"""
        sample_width = 50

        traces = []

        for target in targets:
            x_min = int(target.mean - half_n_sigma * target.stddev)
            x_max = int(target.mean + half_n_sigma * target.stddev)
            width = x_max - x_min

            sample_data = self.data[x_min:x_max, :]
            sx, sy = sample_data.shape

            max_positions = []
            max_index = []
            for y in np.linspace(0, sy - ends_pix_spacing, n_samples, dtype=int):

                sub_sample = sample_data[:, y:y + n_samples]

                sub_median = np.median(sub_sample, axis=1)
                sub_x_axis = range(len(sub_median))

                sub_argmax = np.argmax(sub_median)
                sub_max = np.max(sub_median)

                if do_gaussian_fit:
                    """Will leave this here just in case someone wants to experiment"""
                    gauss_init = models.Gaussian1D(amplitude=sub_max, mean=sub_argmax, stddev=1.)
                    fit_gaussian = fitting.LevMarLSQFitter()
                    gauss = fit_gaussian(gauss_init, sub_x_axis, sub_median)
                    max_positions.append(gauss.mean.value + x_min)
                    max_index.append(y + int(sample_width / 2.))
                else:
                    max_positions.append(sub_argmax + x_min)
                    max_index.append(y + int(sample_width / 2.))
            """chebyshev fitting for defining the trace"""
            chebyshev_init = models.Chebyshev1D(2, domain=[0, sy])
            fit_cheb = fitting.LinearLSQFitter()
            cheb = fit_cheb(chebyshev_init, max_index, max_positions)
            # current_trace
            traces.append([cheb, width])
        return traces

    def extract(self, traces):
        """Extracts the spectra and returns an stack of science data and lamps

        Needs to add background subtraction

        Args:
            traces (list):

        Returns:

        """
        if len(traces) > 0:
            '''loop through traces'''
            sci_pack = []
            for trace in traces:
                '''Initial checks'''
                extracted_object = []
                if self.science_object.lamp_count > 0:
                    all_lamps = []
                    for l in range(self.science_object.lamp_count):
                        lamp = np.array([])
                        all_lamps.append(lamp)
                else:
                    log.warning('There are no lamps available for this Target.')
                chebyshev, width = trace
                if width % 2 == 1:
                    half_width = int((width - 1) / 2)
                else:
                    half_width = int(width / 2)
                '''Getting data shape'''
                x, y = self.data.shape
                sci = []
                background = []
                for i in range(y):
                    x_min = int(round(chebyshev(i))) - half_width
                    x_max = int(round(chebyshev(i))) + half_width
                    data_point = np.sum(self.data[x_min:x_max, i])
                    sci.append(data_point)
                    # print(x_min, x_max,data_point)
                    for e in range(self.science_object.lamp_count):
                        lamp_point = np.sum(self.lamps_data[e][x_min:x_max, i])
                        all_lamps[e] = np.append(all_lamps[e], lamp_point)
                extracted_object.append(np.array(sci))
                if self.science_object.lamp_count > 0:
                    for l in range(self.science_object.lamp_count):
                        extracted_object.append(np.array(all_lamps[l]))
                # plt.plot(sci)
                # plt.plot(all_lamps[0])
                # plt.plot(all_lamps[1])
                # plt.show()
                # print(chebyshev, width, x, y)
                sci_pack.append(extracted_object)
            return sci_pack
        else:
            log.error("There are no traces discovered here!!.")
            return []

    def wavelength_calibration(self, extracted_targets):
        if extracted_targets:
            for target in extracted_targets:
                '''First element in target is always the science target'''
                sci_target = target[0]
                lamps = target[1:]
                calibrated_lamps = self.calibrate_lamps_with_template(lamps)
                # print(len(target))
            print(len(extracted_targets))
        else:
            log.error("Empty input data")

    def calibrate_lamps_with_template(self, lamps):
        log.info("Calibrating lamps!")
        for e in range(len(lamps)):
            log.debug(self.science_object.lamp_type[e])
            lamp_header = self.lamps_header[e]
            lamp_header['COMMENT'] = 'Extracted lamp'
            lamp_hdu = fits.PrimaryHDU(np.array(lamps[e]), lamp_header)
            out_file = 'lamp-' + self.science_object.lamp_type[e] + '.fits'
            lamp_hdu.writeto(out_file, clobber=True)

            plt.plot(lamps[e], label=self.science_object.lamp_type[e])
        plt.legend(loc='best')
        plt.show()
        return True


class IdentifiedTarget:
    def __init__(self, amplitude, mean, stddev):
        self.amplitude = amplitude
        self.mean = mean
        self.stddev = stddev


'''
class Templates:
    def __init__(self):'''
