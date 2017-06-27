# -*- coding: utf8 -*-
"""Module that handles intermediate processess

This module performs intermediate steps in the process to obtain a wavelength calibrated spectrum
the __call__ method returns a list contained uni-dimensional data extracted and the modified instance of ScienceObject
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


import logging
import os

import numpy as np
import matplotlib
matplotlib.use('GTK3Agg')

from astropy.io import fits
from astropy.modeling import models, fitting
from goodman_ccd.core import add_wcs_keys
from matplotlib import pyplot as plt


# FORMAT = '%(levelname)s:%(filename)s:%(module)s: %(message)s'
# log.basicConfig(level=log.INFO, format=FORMAT)
log = logging.getLogger('redspec.process')


class Process(object):
    """Set of tools for extracting Goodman High Throughput Spectrograph data.

    This class needs a path parsed as an attribute of args (arguments object) and a ScienceObject object.
    The ScienceObject class is defined in the module redspec.py and in principle the current module was made to service
    it. This class has been tested to work very well independently given that you successfully provide the previously
    mentioned ScienceObject and respective path to the data location.

    """

    def __init__(self, sci_obj, args):
        """Initialize Process class

        It needs the sci_obj (ScienceObject instance) that contain information regarding the images being processed,
        such as, filename, science target name, etc. It also uses the arguments of the program that also are a class
        instance. For initialization it only uses the source path, or where the files are located.

        Args:
            sci_obj (object): Instance of ScienceObject.
            args (object): Runtime arguments.
        """
        self.args = args
        self.science_object = sci_obj
        self.path = self.args.source
        self.data = fits.getdata(os.path.join(self.path,
                                              self.science_object.file_name))
        self.header = add_wcs_keys(fits.getheader(os.path.join(self.path,
                                                               self.science_object.file_name)))
        self.lamps_data = []
        self.lamps_header = []
        self.close_targets = False
        self.region = None
        self.targets = None
        self.traces = None
        self.extracted_data = None

    def __call__(self, extract_lamps=True):
        """Call method for Process Class

        This method will call the other class' methods that will identify, trace and extract the data. It will also
        modify attributes of the science_object attribute (instance of ScienceObject) and the return it along with the
        extracted data.

        Args:
            extract_lamps (bool): Whether it is necessary to extract the lamp to find the wavelength solution. If it's
             set to False it is assumed that a previously found wavelength solution will be applied.

        Returns:
            extracted (list): Contains two elements. self.extracted_data (object) and self.science_object (object)

        """

        log.info('Processing Science File : %s', self.science_object.file_name)
        if extract_lamps:
            if self.science_object.lamp_count > 0:
                for lamp_index in range(self.science_object.lamp_count):
                    lamp_data = fits.getdata(os.path.join(self.path, self.science_object.lamp_file[lamp_index]))
                    # lamp_type = self.science_object.lamp_type[lamp_index]
                    self.lamps_data.append(lamp_data)
                    lamp_header = add_wcs_keys(
                        fits.getheader(os.path.join(self.path, self.science_object.lamp_file[lamp_index])))
                    self.lamps_header.append(lamp_header)
            else:
                log.error('There Are no lamps available for the target: %s', self.science_object.name)
                return [None, None]

        self.targets = self.identify_spectra()
        if self.targets is not None:
            self.traces = self.trace(self.targets)
            self.extracted_data = self.extract(self.traces)

        extracted = [self.extracted_data, self.science_object]
        return extracted

    def identify_spectra(self):
        """Identify spectra in 2D (image) data

        This method takes a sample of the image and averages it in the dispersion direction, then analyze the spatial
        direction in search for peaks and count them as candidates to targets. In order to validate the target it will
        do a Voigt profile fit and if the new center deviates too much from the original location it will be discarded.

        Returns:
            identified_targets (list): Each element being and IdentifiedTarget instance

        """
        try:
            y_size, x_size = self.data.shape
        except ValueError:
            log.error('Image is not 2D')
            return None
        self.region = np.ones(y_size)
        half_width = int(0.03 * x_size)
        sample_loc = int(x_size / 2.)
        # print x_size, y_size, int(0.05 * x_size)
        sample = np.median(self.data[:, sample_loc:sample_loc + 2 * half_width], axis=1)
        sample_median = np.median(sample)
        sample_std = np.std(sample)

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
                        # print(voigt)
                        if abs(voigt.x_0.value - spectrum) < 5:
                            identified_targets.append(IdentifiedTarget(
                                amplitude=voigt.amplitude_L.value,
                                mean=voigt.x_0.value,
                                stddev=voigt.fwhm_L.value,
                                fwhmg=voigt.fwhm_G.value,
                                raw_width=sample_width,
                                sample_loc=sample_loc))
                            # fitted_model = voigt
                            if self.args.debug_mode:
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

        if self.args.debug_mode:
            plt.title(self.science_object.name)
            plt.axhline(sample_median, color='m', label='Median')
            plt.axhline(1.1 * sample_median, color='c', label='110% Median')
            plt.plot(sample, label='Data')
            plt.xlabel('Pixel (Spatial Direction)')
            plt.ylabel('Intensity')
            plt.legend(loc='best')
            plt.show()
        # for target in identified_targets:
            # print target
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
        half_n_sigma = 5
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
            # target.print_all()
            # x_min = int(target.mean - half_n_sigma * target.stddev)
            x_min = int(target.mean - target.raw_width / 2.)
            # x_max = int(target.mean + half_n_sigma * target.stddev)
            x_max = int(target.mean + target.raw_width / 2.)
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

                # print sub_median
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
                    if False:
                        plt.plot(sub_x_axis, sub_median, color='r', label='Data')
                        plt.plot(sub_x_axis, voigt(sub_x_axis), color='g', label='Voigt Fit')
                        plt.title(self.science_object.name)
                        plt.xlabel('Pixel (spatial direction)')
                        plt.ylabel('Intensity')
                        plt.legend(loc='best')
                        plt.show()

                else:
                    max_positions.append(sub_argmax + x_min)
                    max_index.append(index_y + int(sample_width / 2.))

            # chebyshev fitting for defining the trace
            if np.std(max_positions) < width:
                chebyshev_init = models.Chebyshev1D(2, domain=[0, sample_y])
                fit_cheb = fitting.LinearLSQFitter()
                cheb = fit_cheb(chebyshev_init, max_index, max_positions)


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

                traces.append([cheb, width, self.region])
                self.science_object.update_no_targets(add_one=True)

                if self.args.debug_mode:
                    fig1 = plt.figure(1)
                    fig1.canvas.set_window_title('Trace')
                    plt.title(self.science_object.name)
                    plt.imshow(self.data, clim=(5, 150), cmap='cubehelix', origin='lower', interpolation='nearest')
                    plt.plot(max_index, max_positions, marker='o', color='y', label='Sampled values')
                    plt.plot(cheb(range(0, sample_y)), linestyle='--', color='r', label='Trace fit (Cheb)')
                    # plt.plot(max_index, max_positions, )
                    plt.xlabel('Dispersion Direction')
                    plt.ylabel('Spatial Direction')
                    plt.legend(loc='best')
                    plt.xlim((0, self.data.shape[1]))
                    plt.ylim((0, self.data.shape[0]))
                    plt.tight_layout()
                    plt.show()

                    half_width = int(0.03 * self.data.shape[1])
                    sample_loc = int(self.data.shape[1] / 2.)
                    # print x_size, y_size, int(0.05 * x_size)
                    sample = np.median(self.data[:, sample_loc:sample_loc + 2 * half_width], axis=1)

                    fig2 = plt.figure(2)
                    fig2.canvas.set_window_title('Masked Regions')

                    plt.title("Masked Regions for Background extraction\n%s" % self.science_object.name)
                    plt.xlabel("Spatial Direction")
                    plt.ylabel("Intensity")
                    for target in targets:
                        plt.axvline(target.mean, color='k', linestyle='--', label='Target Mean')
                    plt.plot(sample, color='k', label='Data Sample')
                    # plt.plot(self.region)
                    limits = []
                    for i in range(len(self.region) - 1):
                        if self.region[i] != self.region[i + 1]:
                            if limits == [] or len(limits) % 2 == 0:
                                limits.append(i + 1)
                            else:
                                limits.append(i)
                    # check that the number of limits is even
                    if len(limits) % 2 == 1:
                        log.error('Uneven number of limits.')
                        log.info('Removing last limit.')
                        limits.pop(-1)
                    # print(limits)
                    # colors = ['red', 'green', 'blue']
                    # colors_index = 0
                    plt.axvspan(None, None, color='r', alpha=0.3, label='Background')
                    plt.axvspan(None, None, color='k', alpha=0.3, label='Object')
                    for limit_index in range(0, len(limits), 2):
                        if self.region[limits[limit_index]] == -1:
                            color = 'k'
                        elif self.region[limits[limit_index]] == 0:
                            color = 'r'
                        else:
                            color = 'w'
                        print(limits, limit_index)
                        plt.axvspan(limits[limit_index], limits[limit_index + 1], color=color, alpha=0.3)
                        # colors_index += 1
                        # if colors_index == 3:
                            # colors_index = 0
                    # Watch out here! it works but might be an error.
                    # plt.xlim([target.mean - 200, target.mean + 200])
                    # plt.savefig("background-extraction-zones" + self.science_object.name + "_2.png", dpi=300)
                    plt.legend(loc='best')
                    plt.show()

            else:
                log.error('Target at %s discarded', int(target.mean))

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
            log.warning("Region %d:%d is already masked", mask_min, mask_max)
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
            sci_pack = SciencePack()

            for trace_index in range(len(traces)):
                # Initial checks
                history_headers = []
                new_header = self.header.copy()
                if self.science_object.lamp_count > 0:
                    all_lamps = [np.array([]) for lamp_index in range(self.science_object.lamp_count)]
                else:
                    log.warning('There are no lamps available for this Target.')

                # Extraction of data parsed as argument
                chebyshev, width, _ = traces[trace_index]
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
                        hist = "Defining background extraction zone [%s:%s] for target." % (limits[limit_index],
                                                                                            limits[limit_index + 1])
                        history_headers.append(hist)
                        log.debug(hist)

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
                    apnum1 = '%s %s %s %s'%(trace_index + 1, 1, x_min, x_max)
                    if i == int(data_y/2):
                        hist = 'Aperture for extraction [%s:%s] at %s' % (x_min, x_max, i)
                        history_headers.append(hist)
                        log.debug(hist)
                        log.debug('APNUM1 = %s', apnum1)

                    # If there are background extraction zones here are prepared to subtract
                    if len(background) > 1:
                        background_part = []
                        for back in background:
                            part = np.sum(self.data[back[0]:back[1], i])
                            # background_median = np.median(self.data[back[0]:back[1], i])
                            # part2 = abs(x_max - x_min) * background_median
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
                    # print('DATA POINT ', np.sum(self.data[x_min:x_max, i]), background_data)
                    sci.append(data_point)

                    # Lamp extraction
                    if len(self.lamps_data) > 0:
                        for limit_index in range(self.science_object.lamp_count):
                            lamp_point = np.sum(self.lamps_data[limit_index][x_min:x_max, i])
                            all_lamps[limit_index] = np.append(all_lamps[limit_index], lamp_point)
                # Construction of extracted_object (to be returned)
                # extracted_object.append(np.array(sci))
                sci_pack.add_data(np.array(sci))
                # if int(trace_index + 1) > 1:
                #     new_header.rename_keyword('APNUM1', 'APNUM%s' % str(int(trace_index + 1)))
                new_header['APNUM1'] = apnum1
                if history_headers != []:
                    for hist in history_headers:
                        new_header['HISTORY'] = hist
                # headers.append(new_header)
                sci_pack.add_header(new_header)
                if len(self.lamps_data) > 0:
                    for lamp_index in range(self.science_object.lamp_count):
                        # extracted_object.append(np.array(all_lamps[lamp_index]))
                        sci_pack.add_lamp(np.array(all_lamps[lamp_index]))
                        self.lamps_header[lamp_index]['APNUM1'] = apnum1
                        sci_pack.add_lamp_header(self.lamps_header[lamp_index])
                        # headers.append(self.lamps_header[lamp_index])
                #
                # Plot background subtraction
                if self.args.debug_mode:
                    # sci_sample = sci[int(len(sci) / 2.) - 30:int(len(sci) / 2.) + 30]
                    fig = plt.figure(1)
                    fig.canvas.set_window_title('Subtraction')
                    plt.title('Background Subtraction\n' + self.science_object.name)
                    plt.xlabel('Pixels (dispersion direction)')
                    plt.ylabel('Intensity (counts)')
                    plt.xlim((0, len(sci)))
                    # plt.yscale('log')
                    plt.plot(sci, color='k', alpha=1, label='Background Subtracted')
                    plt.plot(unsubtracted, color='k', alpha=0.5, label='Unsubtracted')
                    plt.plot(subtracted_background, color='r', label='Background')
                    # plt.plot(all_lamps[0],label='lamp 1')
                    # plt.plot(all_lamps[1],label='lamp 2')
                    plt.legend(loc='best')
                    plt.tight_layout()
                    # plt.savefig('background-subtraction_'
                    #            + self.science_object.name
                    #            + '_'
                    #            + str(int(chebyshev(10)))
                    #            + '.png', dpi=300)
                    plt.show()

            return sci_pack
        else:
            log.error("There are no traces discovered here!!.")
            return None


class IdentifiedTarget(object):
    """Allows for easy storage and manipulation of the targets found.

    """

    def __init__(self, amplitude, mean, stddev, fwhmg=0, raw_width=None, sample_loc=None):
        """Initialization of class

        This class stores all the relevant information of an identified target in an image. The idea is to make use
        of the tools that classes provide to make it easy to handle.

        Args:
            amplitude (float): Peak value of the target's sample
            mean (float): Voigt fit center
            stddev (float): Full width at half maximum of Lorentzian component
            fwhmg (float): Full width at half maximum of Gaussian component
            raw_width (int): Width of the sample in dispersion axis.
            sample_loc (int): Location of the sample in the dispersion direction.
        """
        self.amplitude = amplitude
        self.mean = mean
        self.stddev = stddev
        self.fwhm_g = fwhmg
        self.raw_width = raw_width
        self.sample_loc = sample_loc

    def print_all(self):
        """Prints all the attributes of the identified target

        This is mostly for debugging purposes.

        """
        log.debug('Amplitude: %s', self.amplitude)
        log.debug('Mean: %s', self.mean)
        log.debug('Std Deviation: %s', self.stddev)
        log.debug('FWHM (gauss): %s', self.fwhm_g)
        log.debug('Raw Width: %s', self.raw_width)
        log.debug('Sample Location: %s', self.sample_loc)


class SciencePack(object):
    """Science data packaging

    This class is designed to work as a packer of extracted data. The attributes are separated list for data, headers
    lamps data and lamps header. Their only matching mechanism is their index. therefore if you have added three science
    spectra and two comparison lamps, in order to get, say, de second one, you would have to get SciencePack.data[1]
    SciencePack.headers[1]. The same thing for the lamps.
    """

    def __init__(self):
        """SciencePack initialization

        A SciencePack instance will be initialized having four empty list, lists that will store data an headers.

        Since data and headers are added in separated steps, there is no checks as for whether its been added
        consistently. Developer should be careful with this.

        Attributes:
            self.data (list): Stores science target data
            self.headers (list): Stores science target headers
            self.lamps_data (list): Stores comparison lamps data
            self.lamps_headers (list): Stores comparison lamps headers.

        """
        self.data = []
        self.headers = []
        self.lamps_data = []
        self.lamps_headers = []

    def add_data(self, new_data):
        """Appends science data"""
        self.data.append(new_data)

    def add_header(self, new_header):
        """Appends science header"""
        self.headers.append(new_header)

    def add_lamp(self, new_lamp):
        """Appends comparison lamp data"""
        self.lamps_data.append(new_lamp)

    def add_lamp_header(self, new_lamp_header):
        """Appends comparison lamp header"""
        self.lamps_headers.append(new_lamp_header)

    def check_consistency(self):
        """Check that all stored data is consistent

        There should be the same number of _data_ and their respective
        _headers_, the same for _lamps\_data_ and _lamps\_headers_

        """
        if len(self.data) == len(self.headers):
            log.debug('Science data and headers are consistent')
        else:
            log.error('Science data and headers are not consistent')
            log.error('Data: {:d}, Headers: {:d}'.format(len(self.data),
                                                         len(self.headers)))
        if len(self.lamps_data) == len(self.lamps_headers):
            log.debug('Lamps data and headers are consistent')
        else:
            log.error('Lamps data and headers are not consistent')
            log.error('Lamps: {:d}, Headers: {:d}'.format(len(self.data),
                                                         len(self.headers)))