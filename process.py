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
        self.targets = self.identify_spectra()
        self.traces = self.trace(self.targets)

    def identify_spectra(self):
        """Finds the location of the spectrum or spectra in case there are more

        Performs an analysis along the spatial direction to determine the presence of important peaks wich will be
        identified as science targets (or objects)

        I doesn't take arguments because it will read the <object> self.science_object that contains all the
        necessary data.

        Returns:

        """
        x, y = self.data.shape
        sample_data = self.data[:, y / 2. - 100:y / 2. + 100]
        median_data = np.median(sample_data, axis=1)
        x_axis = np.linspace(0, len(median_data), x)
        std = np.std(median_data)
        median = np.median(median_data)
        threshold = median + std
        keep_searching = True
        all_candidates = []
        while keep_searching:
            candidate = []
            candidate_index = []
            found = False
            for i in range(len(median_data)):
                if median_data[i] > threshold:
                    found = True
                    candidate.append(median_data[i])
                    candidate_index.append(i)
                elif found:
                    all_candidates.append([candidate, candidate_index])
                    candidate = []
                    candidate_index = []
                    found = False
                if i == len(median_data) - 1:
                    keep_searching = False
        identified_targets = []
        for single_candidate in all_candidates:
            cmax = np.argmax(single_candidate[0])
            max_val = np.max(single_candidate[0]) / 2.
            max_pos = single_candidate[1][cmax]
            print("max-val", max_val, "max-pos", max_pos)
            gauss_init = models.Gaussian1D(amplitude=max_val, mean=max_pos, stddev=1)
            fit_gaussian = fitting.LevMarLSQFitter()
            gauss = fit_gaussian(gauss_init, x_axis, median_data)
            identified_targets.append(IdentifiedTarget(gauss.amplitude.value, gauss.mean.value, gauss.stddev.value))
            # print(gauss)
            # plt.axhspan(gauss.mean.value - 1.5 * gauss.stddev.value, gauss.mean.value + 1.5 * gauss.stddev.value)
            plt.plot(x_axis, gauss(x_axis))
            # plt.axvspan(gauss.mean.value-2.355/2.*gauss.stddev.value,gauss.mean.value+2.355/2.*gauss.stddev.value)

        """Remove afterwards"""

        plt.plot(median_data)
        # print(np.median(median_data))
        # plt.axhline(threshold)
        # plt.plot(range(len(median_data)),)
        # plt.plot(first_derivative)
        # plt.plot(second_derivative)

        # plt.imshow(data, cmap='gray', clim=(5, 150))
        # plt.show()
        plt.savefig(self.path + 'gauss_' + self.science_object.name + '_' + str(int(gauss.mean.value)) + '.png',
                    dpi=300)
        plt.clf()
        # print(x, y)
        return identified_targets

    def trace(self, targets):
        print(len(targets))
        for target in targets:
            sample_data = self.data[target.mean - 2 * target.stddev:target.mean + 2 * target.stddev, :]
            """The code below does not belong here. It must go to an extraction function"""
            median_spectra = np.median(sample_data, axis=0)
            file_name = self.path + 'spectra_' + self.science_object.name + '_' + str(int(target.mean)) + '.png'
            # plt.imshow(sample_data, cmap='gray', clim=(5, 150))
            # plt.show()
            # plt.clf()
            plt.plot(median_spectra)
            plt.xlabel("Pixel")
            plt.ylabel("Intensity")
            plt.savefig(file_name, dpi=300)
            # plt.show()
            plt.clf()
            print(target.mean)
        return True


class IdentifiedTarget:
    def __init__(self, amplitude, mean, stddev):
        self.amplitude = amplitude
        self.mean = mean
        self.stddev = stddev
