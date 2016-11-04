"""Module that handles intermediate processess

This module performs intermediate steps in the process to obtain a wavelength calibrated spectrum
the __call__ method returns a list contained uni-dimensional data extracted and the modified instance of ScienceObject
"""
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting
import logging as log


FORMAT = '%(levelname)s:%(filename)s:%(module)s: %(message)s'
log.basicConfig(level=log.DEBUG, format=FORMAT)


class Process(object):
    """Set of tools for extracting Goodman High Throughput Spectrograph data.

    This class needs a path parsed as a string and a ScienceObject object. The ScienceObject class is defined in the
    module redspec.py and in principle the current module was made to service it. This class has been tested to work
    very well independently given that you successfully provide the previously mentioned ScienceObject and respective
    path to the data location.

    """

    def __init__(self, source_path, sci_obj, args):
        self.args = args
        self.science_object = sci_obj
        self.path = source_path
        self.data = fits.getdata(self.path + self.science_object.file_name)
        self.header = self.add_wcs_keys(fits.getheader(self.path + self.science_object.file_name))
        self.lamps_data = []
        self.lamps_header = []
        self.close_targets = False
        self.region = None
        self.targets = None
        self.traces = None
        self.extracted_data = None

    def __call__(self, extract_lamps=True):

        log.info('Processing Science File : %s', self.science_object.file_name)
        if extract_lamps:
            if self.science_object.lamp_count > 0:
                for lamp_index in range(self.science_object.lamp_count):
                    lamp_data = fits.getdata(self.path + self.science_object.lamp_file[lamp_index])
                    # lamp_type = self.science_object.lamp_type[lamp_index]
                    self.lamps_data.append(lamp_data)
                    lamp_header = self.add_wcs_keys(
                        fits.getheader(self.path + self.science_object.lamp_file[lamp_index]))
                    self.lamps_header.append(lamp_header)
            else:
                log.error('There Are no lamps available for the target: %s', self.science_object.name)
                return [None, None]

        self.targets = self.identify_spectra()
        self.traces = self.trace(self.targets)
        self.extracted_data = self.extract(self.traces)

        return [self.extracted_data, self.science_object]

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
        # TODO(simon) Improve target contamination rejections
        data_x, data_y = self.data.shape
        self.region = np.ones(data_x)
        sample_data = np.median(self.data[:, int(data_y / 2.) - 100:int(data_y / 2.) + 100], axis=1)
        # x_axis = np.linspace(0, len(sample_data), x)
        x_axis = range(0, len(sample_data), 1)
        # print(x_axis)
        std = np.nanstd(sample_data)
        log.debug('Standard Deviation: %s', std)
        median = np.nanmedian(sample_data)
        log.debug('Median: %s', median)
        threshold = median + std
        log.debug('Selection Threshold: %s', threshold)
        keep_searching = True
        all_candidates = []
        if self.args.plots_enabled:
            plt.plot(sample_data)
            plt.axhline(threshold)
            plt.show()
            plt.clf()

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
        # print('all ', all_candidates, len(all_candidates))
        if self.args.plots_enabled:
            for cand in all_candidates:
                plt.plot(cand[1], cand[0])
            plt.plot(sample_data, linestyle='--')
            plt.show()
            plt.clf()
        for single_candidate in all_candidates:
            # if len(single_candidate[0]) <= 4:
                # log.debug('Not enough data points to make it to a Science Target')
            try:
                cmax = np.argmax(single_candidate[0])
                max_val = np.max(single_candidate[0])
                max_pos = single_candidate[1][cmax]

                model_to_fit = 'voigt'

                if model_to_fit == 'gauss':
                    gauss_init = models.Gaussian1D(amplitude=max_val, mean=max_pos, stddev=4)
                    fit_gaussian = fitting.LevMarLSQFitter()
                    gauss = fit_gaussian(gauss_init, single_candidate[1], single_candidate[0])
                    identified_targets.append(
                        IdentifiedTarget(gauss.amplitude.value, gauss.mean.value, gauss.stddev.value))
                    fitted_model = gauss
                elif model_to_fit == 'lorentz':
                    lorentz_init = models.Lorentz1D(amplitude=max_val, x_0=max_pos, fwhm=8)
                    log.debug('Amplitude: %s X_0: %s FWHM: %s', max_val, max_pos, 8)
                    fit_lorentz = fitting.LevMarLSQFitter()
                    lorentz = fit_lorentz(lorentz_init, single_candidate[1], single_candidate[0])
                    # print(lorentz)
                    identified_targets.append(
                        IdentifiedTarget(lorentz.amplitude.value, lorentz.x_0.value, lorentz.fwhm.value))
                    fitted_model = lorentz
                elif model_to_fit == 'voigt':
                    voigt_init = models.Voigt1D(x_0=max_pos, amplitude_L=max_val, fwhm_L=8, fwhm_G=8)
                    fit_voigt = fitting.LevMarLSQFitter()
                    # voigt = fit_voigt(voigt_init, x_axis, sample_data)
                    voigt = fit_voigt(voigt_init, single_candidate[1], single_candidate[0])
                    # print(voigt)
                    identified_targets.append(IdentifiedTarget(
                        voigt.amplitude_L.value,
                        voigt.x_0.value,
                        voigt.fwhm_L.value,
                        voigt.fwhm_G.value))
                    fitted_model = voigt

                # x_0 : float
                # Position of the peak
                # amplitude_L : float
                # The Lorentzian amplitude
                # fwhm_L : float
                # The Lorentzian full width at half maximum
                # fwhm_G : float
                # The Gaussian full width at half maximum

                residuals = sample_data - fitted_model(x_axis)
                if self.args.plots_enabled:
                    plt.plot(x_axis, sample_data)
                    plt.plot(x_axis, fitted_model(x_axis))
                    plt.plot(x_axis, residuals)
                    # plt.axvline(gauss.x_0.value)
                    plt.show()
                    plt.clf()
            except TypeError as err:
                log.error(err)
                log.info('Object candidate will be ignored. Most likely does not have enough points.')


        if len(identified_targets) > 0:
            log.info('Identified %s targets.', len(identified_targets))
            if len(identified_targets) > 1:
                # check how close the targets are to each other
                plt.plot(x_axis, sample_data)
                center = (identified_targets[0].mean + identified_targets[1].mean)/2.
                mean_width = (identified_targets[0].stddev + identified_targets[1].stddev)/2.
                plt.axvspan(center - mean_width, center + mean_width, color='g', alpha=.5)
                for target_index in range(0, len(identified_targets)):
                    print identified_targets[target_index].stddev, identified_targets[target_index].fwhm_g
                    half = identified_targets[target_index].stddev
                    plt.axvspan(identified_targets[target_index].mean - half,
                                identified_targets[target_index].mean + half,
                                color='r', alpha=0.5)
                plt.show()
                return identified_targets
            else:
                return identified_targets
        else:
            log.error('No Target identified')

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
        This method also mask the regions that will be extracted as science target (or spectrum) and the regions
        that will be used for background subtractions. At first two regions are defined for each science target
        at each side by an offset defined by the background_offset variable and the same width as the science target's
        aperture. If there are background subtraction regions that fall into a science target they will be ignored.
        In that case only one background region will be used.
        If both subtraction regions are valid they will be averaged and then subtracted.
        If no background region are valid, no background subtraction will happen.

        Note:
            It is not worth doing gaussian fits to the sub samples since it will not improve the result. The only case
            this might be useful is the case when someone does an interpolation along the spatial direction.

        Args:
            targets (list): Each element is a class that stores the parameters of the gaussian fitted to the data in
                        previous steps. This data tells the location in the image of the target or targets.

        Returns:
            traces (list): Every element is a list with two elements. The first is the fitted Chebyshev and the
                       second and last is the width to be extracted.

        """

        # True enables gaussian fits and false disables it. Gaussian fit was left for experimental purposes
        do_gaussian_fit = False
        # Voigt fits better
        do_voigt_fit = True
        # half of number of sigmas to be sub_sampled
        half_n_sigma = 3
        # space allowance in pixels for extreme of images
        ends_pix_spacing = 50
        # Number of samples to be done in the dispersion direction
        n_samples = 50
        # sample width must be smaller or equal than ends_pix_spacing
        sample_width = 50
        # Background offsets from target
        # background_offset = 5

        traces = []
        regions = []

        for target in targets:
            x_min = int(target.mean - half_n_sigma * target.stddev)
            x_max = int(target.mean + half_n_sigma * target.stddev)
            width = x_max - x_min
            background_offset = 1.5 * width
            # target value in mask is -1, for background 0 for non-masked is 1
            self.mask(x_min, x_max, -1)

            regions.append([x_min, x_max, width])

            sample_data = self.data[x_min:x_max, :]
            sample_y = sample_data.shape[1]

            max_positions = []
            max_index = []
            for index_y in np.linspace(0, sample_y - ends_pix_spacing, n_samples, dtype=int):

                sub_sample = sample_data[:, index_y:index_y + n_samples]

                sub_median = np.median(sub_sample, axis=1)
                sub_x_axis = range(len(sub_median))

                sub_argmax = np.argmax(sub_median)
                sub_max = np.max(sub_median)

                if do_gaussian_fit:
                    # Will leave this here just in case someone wants to experiment
                    # A Voigt profile is better for this case
                    gauss_init = models.Gaussian1D(amplitude=sub_max, mean=sub_argmax, stddev=1.)
                    fit_gaussian = fitting.LevMarLSQFitter()
                    gauss = fit_gaussian(gauss_init, sub_x_axis, sub_median)
                    max_positions.append(gauss.mean.value + x_min)
                    max_index.append(index_y + int(sample_width / 2.))
                elif do_voigt_fit:
                    voigt_init = models.Voigt1D(x_0=sub_argmax, amplitude_L=sub_max, fwhm_L=8, fwhm_G=8)
                    fit_voigt = fitting.LevMarLSQFitter()
                    voigt = fit_voigt(voigt_init, sub_x_axis, sub_median)
                    max_positions.append(voigt.x_0.value + x_min)
                    max_index.append(index_y + int(sample_width / 2.))
                else:
                    max_positions.append(sub_argmax + x_min)
                    max_index.append(index_y + int(sample_width / 2.))

            # chebyshev fitting for defining the trace
            chebyshev_init = models.Chebyshev1D(2, domain=[0, sample_y])
            fit_cheb = fitting.LinearLSQFitter()
            cheb = fit_cheb(chebyshev_init, max_index, max_positions)
            traces.append([cheb, width])

            if self.args.plots_enabled:
                plt.imshow(self.data, clim=(5, 150))
                plt.plot(max_index, max_positions, color='r')
                plt.show()

        # Mask Background Extraction zones
        # TODO(simon): Define what to do in case no background extraction zone is suitable for use.
        for region in regions:
            x_min, x_max, width = region
            b_low_min = x_min - background_offset - width
            b_low_max = x_min - background_offset
            self.mask(b_low_min, b_low_max, 0)
            b_high_min = x_max + background_offset
            b_high_max = x_max + background_offset + width
            self.mask(b_high_min, b_high_max, 0)

        #
        # plots the masked regions.
        if self.args.plots_enabled:
            plt.title("Masked Regions for Background extraction")
            plt.xlabel("Spatial Direction")
            plt.ylabel("Intensity")
            for target in targets:
                plt.axvline(target.mean)
            plt.plot(self.data[:, 1000])
            plt.plot(self.region)
            limits = []
            for i in range(len(self.region) - 1):
                if self.region[i] != self.region[i + 1]:
                    if limits == [] or len(limits) % 2 == 0:
                        limits.append(i + 1)
                    else:
                        limits.append(i)
            # print(limits)
            colors = ['red', 'green', 'blue']
            colors_index = 0
            for limit_index in range(0, len(limits), 2):
                plt.axvspan(limits[limit_index], limits[limit_index + 1], color=colors[colors_index], alpha=0.3)
                colors_index += 1
                if colors_index == 3:
                    colors_index = 0
            # Watch out here! it works but might be an error.
            plt.xlim([target.mean - 200, target.mean + 200])
            plt.savefig("background-extraction-zones" + self.science_object.name + "_2.png", dpi=300)
            plt.show()

        #
        return traces

    def mask(self, mask_min, mask_max, value):
        """Masks the region that will be extracted

        The class Process has the attribute region (self.region) which is a section across the spatial direction of the
        image being processed. It is filled with ones (1) meaning that is _free zone_. Every time a new science
        target is found a mask will be created with -1 values. This type of mask has priority over all the rest.
        For every science target two background extraction zones are defined and masked with 0 but if they
        interfere with another science target they will be ignored.


        Args:
            mask_min (int): Starting point to the region to be masked
            mask_max (int): Ending point to the region to be masked
            value (int): Value to be replaced. -1 is for science targets and 0 for background regions.

        Returns:
            True or False: True if the masking suceeds or false if it fails.

        """

        if int(value) not in self.region[mask_min:mask_max] and int(-1) not in self.region[mask_min:mask_max]:
            self.region[mask_min:mask_max] = value
            return True
        else:
            log.warning("Region %s:%s is already masked", mask_min, mask_max)
            return False

    def extract(self, traces):
        """Extracts the spectra and returns an stack of science data and lamps

        Apart of traces this function makes use of the class attribute region which contains a mask for the extraction
        zones identified different for science targets and for background zones. The background subtraction is only
        used for the science targets since it would not make sense to do it for lamps.

        Args:
            traces (list): Contains a list per each science target previously detected. Every list has two elements
                that are: Chebyshev fitted to spectrum peaks along the full extension of the spectrum or along the
                full length of the dispersion axis.

        Returns:
            sci_pack (list): The sci_pack is list that contains two main elements. The first contains all the data and
            the second contains the headers in the same order. In the first element, the first is always the science
            target and the following are lamps that are not calibrated yet. Again, the same order for the headers.

        """
        if len(traces) > 0:
            # Loop through traces
            sci_pack = []
            for trace_index in range(len(traces)):
                # Initial checks
                extracted_object = []
                headers = []
                # TODO(simon) this might be simplified in a pythonic way
                if self.science_object.lamp_count > 0:
                    all_lamps = []
                    for lamp_index in range(self.science_object.lamp_count):
                        lamp = np.array([])
                        all_lamps.append(lamp)
                else:
                    log.warning('There are no lamps available for this Target.')

                # Extraction of data parsed as argument
                chebyshev, width = traces[trace_index]
                # log.debug("Offset for Background: %s", offset)
                if width % 2 == 1:
                    half_width = int((width - 1) / 2)
                else:
                    half_width = int(width / 2)

                # Define background subtraction zones
                # TODO(simon): Make the background subtraction zones to follow the trace
                background = []
                limits = []
                for i in range(len(self.region) - 1):
                    if self.region[i] != self.region[i + 1]:
                        log.debug("Found a mask limit.")
                        limits.append(i + 1)

                for limit_index in range(0, len(limits) - 1, 2):
                    if limits[limit_index + 1] - limits[limit_index] == width and\
                                    -1 not in self.region[limits[limit_index]:limits[limit_index + 1]]:
                        background.append([limits[limit_index], limits[limit_index + 1]])
                        log.debug("Defining background extraction zone [%s:%s] for target.",
                                  limits[limit_index],
                                  limits[limit_index + 1])

                # Getting data shape
                data_y = self.data.shape[1]
                sci = []

                unsubtracted = []
                subtracted_background = []

                # Actual extraction
                # Avoid printing inside this loop since it goes through all the columns
                for i in range(data_y):
                    # Define limits of aperture for spectrum
                    x_min = int(round(chebyshev(i))) - half_width
                    x_max = int(round(chebyshev(i))) + half_width
                    apnum1 = '%s %s %s %s'%(trace_index, 1, x_min, x_max)
                    if i == int(data_y/2):
                        log.debug('Aperture for extraction [%s:%s] at %s', x_min, x_max, i)
                        log.debug('APNUM1 = %s', apnum1)

                    # If there are background extraction zones here are prepared to subtract
                    if len(background) > 1:
                        background_part = []
                        for back in background:
                            part = np.sum(self.data[back[0]:back[1], i])
                            background_part.append(part)
                        background_data = np.mean(background_part)
                    elif len(background) == 1:
                        background_data = np.sum(self.data[background[0][0]:background[0][1], i])
                    else:
                        background_data = 0

                    # Stored for debugging process
                    # print background_data
                    subtracted_background.append(background_data)
                    unsubtracted.append(np.sum(self.data[x_min:x_max, i]))

                    data_point = np.sum(self.data[x_min:x_max, i]) - background_data
                    sci.append(data_point)

                    # Lamp extraction
                    if len(self.lamps_data) > 0:
                        for limit_index in range(self.science_object.lamp_count):
                            lamp_point = np.sum(self.lamps_data[limit_index][x_min:x_max, i])
                            all_lamps[limit_index] = np.append(all_lamps[limit_index], lamp_point)
                # Construction of extracted_object (to be returned)
                extracted_object.append(np.array(sci))
                self.header['APNUM1'] = apnum1
                headers.append(self.header)
                if len(self.lamps_data) > 0:
                    for lamp_index in range(self.science_object.lamp_count):
                        extracted_object.append(np.array(all_lamps[lamp_index]))
                        self.lamps_header[lamp_index]['APNUM1'] = apnum1
                        headers.append(self.lamps_header[lamp_index])
                #
                # Plot background subtraction
                if self.args.plots_enabled:
                    plt.title('Background Subtraction\n' + self.science_object.name)
                    plt.xlabel('Pixels (dispersion direction)')
                    plt.ylabel('Intensity (counts)')
                    plt.plot(sci, label='Background Subtracted')
                    plt.plot(unsubtracted, label='Unsubtracted')
                    plt.plot(subtracted_background, label='Background')
                    # plt.plot(all_lamps[0],label='lamp 1')
                    # plt.plot(all_lamps[1],label='lamp 2')
                    plt.legend(loc='best')
                    plt.tight_layout()
                    plt.savefig('background-subtraction_'
                                + self.science_object.name
                                + '_'
                                + str(int(chebyshev(10)))
                                + '.png', dpi=300)
                    plt.show()

                # print(chebyshev, width, x, y)
                #
                sci_pack.append(extracted_object)
                sci_pack.append(headers)
            return sci_pack
        else:
            log.error("There are no traces discovered here!!.")
            return []

    @staticmethod
    def add_wcs_keys(header):
        """Adds generic keyword to the header

        Linear wavelength solutions require a set of standard fits keywords. Later on they will be updated accordingly

        The main goal of putting them here is to have consistent and nicely ordered headers
        """
        try:
            header['BANDID1'] = 'spectrum - background none, weights none, clean no'
            header['APNUM1'] = '1 1 0 0'
            header['WCSDIM'] = 1
            header['CTYPE1'] = 'LINEAR'
            header['CRVAL1'] = 1
            header['CRPIX1'] = 1
            header['CDELT1'] = 1
            header['CD1_1'] = 1
            header['LTM1_1'] = 1
            header['WAT0_001'] = 'system=equispec'
            header['WAT1_001'] = 'wtype=linear label=Wavelength units=angstroms'
            header['DC-FLAG'] = 0
            header['DCLOG1'] = 'REFSPEC1 = non set'
            return header
        except TypeError as err:
            log.error("Can't add wcs keywords to header")
            log.debug(err)




class IdentifiedTarget(object):
    """Allows for easy storage and manipulation of the targets found.

    """

    def __init__(self, amplitude, mean, stddev, fwhmg=0):
        self.amplitude = amplitude
        self.mean = mean
        self.stddev = stddev
        self.fwhm_g = fwhmg

