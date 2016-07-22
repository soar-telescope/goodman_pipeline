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
        self.extraction(self.traces)

    def identify_spectra(self):
        """Finds the location of the spectrum or spectra in case there are more

        Performs an analysis along the spatial direction to determine the presence of important peaks wich will be
        identified as science targets (or objects)

        I doesn't take arguments because it will read the <object> self.science_object that contains all the
        necessary data.

        Returns:

        """
        x, y = self.data.shape
        print(x,y)
        sample_data = self.data[:, int(y / 2.) - 100:int(y / 2.) + 100]
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
        """
        Notes:
            Is not worth doing gaussian fits to the subsamples
        Args:
            targets:

        Returns:

        """

        """True enables gaussian fits and false disables it. Gaussian fit was left for experimental purposes"""
        do_gaussian_fit = True
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
            sx,sy = sample_data.shape

            max_positions = []
            max_index = []
            for y in np.linspace(0,sy-ends_pix_spacing,n_samples,dtype=int):
                print(y)
                sub_sample = sample_data[:, y:y+n_samples]

                sub_median = np.median(sub_sample,axis=1)
                sub_x_axis = range(len(sub_median))

                sub_argmax = np.argmax(sub_median)
                sub_max = np.max(sub_median)

                if do_gaussian_fit:
                    """Will leave this here just in case someone wants to experiment"""
                    gauss_init = models.Gaussian1D(amplitude=sub_max, mean=sub_argmax, stddev=1.)
                    fit_gaussian = fitting.LevMarLSQFitter()
                    gauss = fit_gaussian(gauss_init, sub_x_axis, sub_median)
                    #print(sub_argmax, gauss.mean.value)
                    max_positions.append(gauss.mean.value + x_min)
                    max_index.append(y + int(sample_width/2.))
                else:
                    max_positions.append(sub_argmax + x_min)
                    max_index.append(y + int(sample_width/2.))
            """chebyshev fitting for defining the trace"""
            cheb_x_axis = range(sy)
            chebyshev_init = models.Chebyshev1D(2,domain=[0, sy])
            fit_cheb = fitting.LinearLSQFitter()
            cheb = fit_cheb(chebyshev_init,max_index,max_positions)
            # current_trace
            traces.append([cheb,width])
            print(cheb)
            #plt.plot(cheb_x_axis,cheb(cheb_x_axis),color='r')




                #plt.plot(sub_median)
                #print(y)
            #plt.savefig("spectra-profile-" + self.science_object.name + '_' + str(int(target.mean)) + '.png',dpi=300)
            #plt.show()
            #plt.clf()
            #print(sx,sy)


            """The code below does not belong here. It must go to an extraction function"""
            median_spectra = np.median(sample_data, axis=0)
            file_name = self.path + 'spectra_' + self.science_object.name + '_' + str(int(target.mean)) + '.png'
            #plt.imshow(self.data, cmap='gray', clim=(5, 150))
            #plt.axhspan(x_min,x_max)
            #plt.plot(max_index,max_positions)
            # plt.show()
            # plt.clf()
            # plt.plot(median_spectra)
            #plt.xlabel("Pixel")
            #plt.ylabel("Intensity")
            # plt.savefig(file_name, dpi=300)
            #plt.show()
            #plt.clf()
            #print(target.mean)
        return traces

    def extraction(self,traces):
        x, y = self.data.shape
        spectra = []
        for trace in traces:
            fitted_cheb, width = trace
            if width%2 == 1:
                half_width = int((width-1)/2)
            else:
                half_width = int(width/2)
            extracted_spectrum = []
            for i in range(y):
                cheb_eval = int(round(fitted_cheb(i)))
                x_min = cheb_eval - half_width
                x_max = cheb_eval + half_width

                section = self.data[x_min:x_max,i]
                extracted_spectrum.append(section)
            # plt.imshow(extracted_spectrum,cmap='gray', clim=(5, 150))
            # plt.show()
            # print(fitted_cheb(0),fitted_cheb(1))
            spectra.append(extracted_spectrum)
        return spectra


class IdentifiedTarget:
    def __init__(self, amplitude, mean, stddev):
        self.amplitude = amplitude
        self.mean = mean
        self.stddev = stddev
