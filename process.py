from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting


class Process:
    """

    """

    def __init__(self, source_path, sci_obj):
        self.science_object = sci_obj
        self.path = source_path
        self.data = fits.getdata(self.path + self.science_object.file_name)
        self.header = fits.getheader(self.path + self.science_object.file_name)
        self.lamps_data = []
        if self.science_object.lamp_count > 0:
            for lamp_index in range(self.science_object.lamp_count):
                lamp_data = fits.getdata(self.path + self.science_object.lamp_file[lamp_index])
                lamp_type = self.science_object.lamp_type[lamp_index]
                self.lamps_data.append([lamp_data, lamp_type])
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

        I doesn't take arguments because it will read the <object> self.science_object that contains all the
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
        """
        Notes:
            Is not worth doing gaussian fits to the subsamples
        Args:
            targets:

        Returns:

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
            # print("width ",width)

            sample_data = self.data[x_min:x_max, :]
            sx, sy = sample_data.shape

            max_positions = []
            max_index = []
            for y in np.linspace(0, sy - ends_pix_spacing, n_samples, dtype=int):
                # print(y)
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
                    # print(sub_argmax, gauss.mean.value)
                    max_positions.append(gauss.mean.value + x_min)
                    max_index.append(y + int(sample_width / 2.))
                else:
                    max_positions.append(sub_argmax + x_min)
                    max_index.append(y + int(sample_width / 2.))
            """chebyshev fitting for defining the trace"""
            # cheb_x_axis = range(sy)
            chebyshev_init = models.Chebyshev1D(2, domain=[0, sy])
            fit_cheb = fitting.LinearLSQFitter()
            cheb = fit_cheb(chebyshev_init, max_index, max_positions)
            # current_trace
            traces.append([cheb, width])

        return traces

    def extract(self, traces):
        spectra = []
        for trace in traces:
            """Defines some variables as well as the width that will be extracted"""
            fitted_cheb, width = trace
            if width % 2 == 1:
                half_width = int((width - 1) / 2)
            else:
                half_width = int(width / 2)
            """Spectrum extraction part"""
            print("data type", type(self.data))
            x, y = self.data.shape
            extracted_spectrum = []
            extracted_lamps = []
            for i in range(y):
                cheb_eval = int(round(fitted_cheb(i)))
                x_min = cheb_eval - half_width
                x_max = cheb_eval + half_width
                """Actual extraction of spectrum"""
                spectrum_section = self.data[x_min:x_max, i]
                extracted_spectrum.append(spectrum_section)
                """Lamp extraction part"""
                if self.lamps_data != [] and False:
                    lamp_spectra = [[]] * self.science_object.lamp_count
                    for e in range(self.science_object.lamp_count):
                        lamp_data, lamp_type = self.lamps_data[e]
                        lamp_section = lamp_data[x_min:x_max, i]
                        lamp_spectra[e].append(lamp_section)
                    extracted_lamps = lamp_spectra

            # plt.imshow(extracted_spectrum,cmap='gray', clim=(5, 150))
            # plt.show()
            # print(fitted_cheb(0),fitted_cheb(1))
            spectra.append([np.array(extracted_spectrum), extracted_lamps])

        return spectra

    def wavelength_calibration(self, extracted):
        object_data, lamps_data = extracted
        for object in object_data:
            print("Data shape ", type(object), len(object))
        one_d_spectrum = np.median(object_data, axis=0)
        print("1D spectrum shape", one_d_spectrum.shape)


class IdentifiedTarget:
    def __init__(self, amplitude, mean, stddev):
        self.amplitude = amplitude
        self.mean = mean
        self.stddev = stddev
