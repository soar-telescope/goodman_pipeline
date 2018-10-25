# -*- coding: utf8 -*-
"""Contains the tools to produce a wavelength solution

This module gets the extracted data to produce a wavelength solution, linearize
the spectrum and write the solution to the image's header following the FITS
standard.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import ccdproc
import glob
import logging
import os
import re
import sys

import astropy.units as u

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip
from ccdproc import CCDData
from matplotlib.backends.backend_pdf import PdfPages
from scipy import signal

from ..wcs.wcs import WCS

from ..core import (write_fits)
from ..core import (ReferenceData)


SHOW_PLOTS = False


class WavelengthCalibration(object):
    """Wavelength Calibration Class

    The WavelengthCalibration class is instantiated for each of the science
    images, which are treated as a "science object". In this first release it
    can find a wavelength solution for a given comparison lamp using an
    interactive GUI based on Matplotlib. Although it works very good, for the
    next release there is a plan for creating an independent GUI based on QT in
    order to work better in different screen sizes and other topic such as
    showing warnings, messages and help.

    This class takes 1D spectrum with no wavelength calibration and returns fits
    files with wavelength solutions using the FITS standard for linear
    solutions. Goodman spectra are slightly non-linear therefore they are
    linearized and smoothed before they are returned for the user.

    """

    def __init__(self, args):
        """Wavelength Calibration Class Initialization

        A WavelengthCalibration class is instantiated for each science target
        being processed, i.e. every science image.

        Notes:
            This class violates some conventions as for length and number of
            attributes is concerned. Solving this is part of a prioritary plans
            for next release.

        Args:
            args (Namespace): Runtime arguments.

        """

        self.log = logging.getLogger(__name__)
        self.args = args
        self.poly_order = 3
        self.wcs = WCS()
        self.wsolution = None
        self.wcal_lamp_file = None
        self.sci_target_file = None
        self.n_points = None
        self.n_rejections = None
        self.rms_error = None
        self.cross_corr_tolerance = 5
        # print(self.args.reference_dir)
        self.reference_data = ReferenceData(self.args.reference_dir)
        # self.science_object = science_object
        self.slit_offset = None
        self.interpolation_size = 200
        self.line_search_method = 'derivative'
        # Instrument configuration and spectral characteristics
        self.pixel_size = 15 * u.micrometer
        self.pixel_scale = 0.15 * u.arcsec
        self.goodman_focal_length = 377.3 * u.mm
        self.grating_frequency = None
        self.grating_angle = None
        self.camera_angle = None
        self.serial_binning = None
        self.parallel_binning = None

        self.pixel_count = None
        self.alpha = None
        self.beta = None
        self.center_wavelength = None
        self.blue_limit = None
        self.red_limit = None
        # Interactive wavelength finding
        self.reference_marks_x = []
        self.reference_marks_y = []
        self.raw_data_marks_x = []
        self.raw_data_marks_y = []
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
        self.legends = None
        self.points_ref = None
        self.points_raw = None
        self.line_raw = None
        self.ref_filling_value = 1000
        self.raw_filling_value = 1000
        self.events = True
        self.first = True
        self.evaluation_comment = None
        # self.binning = self.lamp.header[]
        self.pixel_center = []
        # this data must come parsed
        self.path = self.args.source
        # self.science_pack = sci_pack
        # self.sci_filename = self.science_object.file_name
        # self.history_of_lamps_solutions = {}
        self.reference_solution = None

    def __call__(self,
                 ccd,
                 comp_list,
                 object_number=None,
                 wsolution_obj=None,
                 corr_tolerance=15):
        """Call method for the WavelengthSolution Class

        It takes extracted data and produces wavelength calibrated 1D FITS file.
        The call method takes care of the order and logic needed to call the
        different methods. A wavelength solution can be recycled for the next
        science object. In that case, the wavelength solution is parsed as an
        argument and then there is no need to calculate it again. The recycling
        part has to be implemented in the caller function.

        Args:
            ccd (CCDData) a :class:`~astropy.nddata.CCDData` instance
            comp_list (list): Comparison lamps for the science target that will
                be processed here. Every element of this list is an instance of
                :class:`~astropy.nddata.CCDData`.
            object_number (int): In case of multiple detections in a single
                image this number will be added as a suffix before `.fits` in
                order to allow for multiple 1D files. Default value is None.
            wsolution_obj (object): Mathematical model of the wavelength
                solution if exist. If it doesnt is a None
            corr_tolerance (int): `cross_corr_tolerance` stands for cross
                correlation tolerance, in other words, how far the cross
                correlation can be from the global cross correlation. It usually
                increases with the frequency of the grating.

        Returns:
            wavelength_solution (object): The mathematical model of the
                wavelength solution. If it fails to create it will return a
                None element.

        """
        assert isinstance(ccd, CCDData)
        assert isinstance(comp_list, list)

        self.cross_corr_tolerance = corr_tolerance
        self.sci_target_file = ccd.header['GSP_FNAM']

        self.i_fig = None

        self.log.info('Processing Science Target: '
                      '{:s}'.format(ccd.header['OBJECT']))
        if comp_list is not None:
            wavelength_solutions = []
            reference_lamp_names = []
            for self.lamp in comp_list:
                try:
                    self.calibration_lamp = self.lamp.header['GSP_FNAM']
                except KeyError:
                    self.calibration_lamp = ''

                self.raw_pixel_axis = range(self.lamp.shape[0])

                self.lamp_name = self.lamp.header['OBJECT']

                self.log.info('Processing Comparison Lamp: '
                              '{:s}'.format(self.lamp_name))

                self.lines_center = self._get_lines_in_lamp()
                self.spectral = self._get_spectral_characteristics()

                self._automatic_wavelength_solution(
                        corr_tolerance=self.cross_corr_tolerance)

                if self.wsolution is not None:

                    linear_x_axis, self.lamp.data = self._linearize_spectrum(
                        self.lamp.data)

                    self.lamp = self.wcs.write_gsp_wcs(ccd=self.lamp,
                                                       model=self.wsolution)

                    self.lamp = self.add_wavelength_solution(
                        ccd=self.lamp,
                        x_axis=linear_x_axis)

                    self.wcal_lamp_file = self._save_wavelength_calibrated(
                        ccd=self.lamp,
                        original_filename=self.calibration_lamp,
                        index=object_number,
                        lamp=True)

                    wavelength_solutions.append(self.wsolution)
                    reference_lamp_names.append(self.wcal_lamp_file)
                else:
                    self.log.error('It was not possible to get a wavelength '
                                   'solution from lamp '
                                   '{:s} {:s}.'.format(
                                       self.lamp.header['GSP_FNAM'],
                                       self.lamp.header['OBJECT']))
                    continue

            if len(wavelength_solutions) > 1:
                self.log.warning("The current version of the pipeline does not "
                                 "combine multiple solution instead it saves a "
                                 "single version of the science file for each "
                                 "wavelength solution calculated.")
                for i in range(len(wavelength_solutions)):
                    # TODO (simon): Combine Multiple solutions
                    self.wsolution = wavelength_solutions[i]
                    self.wcal_lamp_file = reference_lamp_names[i]
                    self._save_science_data(ccd=ccd, index=i + 1)

            elif len(wavelength_solutions) == 1:
                self.wsolution = wavelength_solutions[0]
                self.wcal_lamp_file = reference_lamp_names[0]
                self._save_science_data(ccd=ccd)
            else:
                self.log.error("No wavelength solution.")

        else:
            print('Data should be saved anyways')

    def add_wavelength_solution(self,
                                ccd,
                                x_axis,
                                evaluation_comment=None):
        """Add wavelength solution to the new FITS header

        Defines FITS header keyword values that will represent the wavelength
        solution in the header so that the image can be read in any other
        astronomical tool. (e.g. IRAF)

        Notes:
            This method also saves the data to a new FITS file, This should be
            in separated methods to have more control on either process.

        Args:
            ccd (CCDData) Instance of :class:`~astropy.nddata.CCDData`
            x_axis:
            evaluation_comment (str): A comment with information regarding the
              quality of the wavelength solution

        Returns:
            ccd (CCDData) A :class:`~astropy.nddata.CCDData` instance with
              linear wavelength solution on it.

        """
        # TODO (simon): Move this to WCS class
        # print(self.n_points, self.n_rejections, self.rms_error)
        # if not all([self.n_points, self.n_rejections, self.rms_error]):
        #     self.rms_error, self.n_points, self.n_rejections = \
        #         self.evaluate_solution()

        ccd.header.set('GSP_WRMS', value=self.rms_error)
        ccd.header.set('GSP_WPOI', value=self.n_points)
        ccd.header.set('GSP_WREJ', value=self.n_rejections)

        if evaluation_comment is None:
            self.evaluation_comment = 'Lamp Solution RMSE = {:.3f} ' \
                                      'Npoints = {:d}, ' \
                                      'NRej = {:d}'.format(self.rms_error,
                                                           self.n_points,
                                                           self.n_rejections)

        new_crpix = 1
        new_crval = x_axis[new_crpix - 1]
        new_cdelt = x_axis[new_crpix] - x_axis[new_crpix - 1]

        ccd.header['BANDID1'] = 'spectrum - background none, weights none, ' \
                                'clean no'
        # ccd.header['APNUM1'] = '1 1 1452.06 1454.87'
        ccd.header['WCSDIM'] = 1
        ccd.header['CTYPE1'] = 'LINEAR  '
        ccd.header['CRVAL1'] = new_crval
        ccd.header['CRPIX1'] = new_crpix
        ccd.header['CDELT1'] = new_cdelt
        ccd.header['CD1_1'] = new_cdelt
        ccd.header['LTM1_1'] = 1.
        ccd.header['WAT0_001'] = 'system=equispec'
        ccd.header['WAT1_001'] = 'wtype=linear label=Wavelength units=angstroms'
        ccd.header['DC-FLAG'] = 0
        # print(self.calibration_lamp)
        ccd.header['DCLOG1'] = 'REFSPEC1 = {:s}'.format(self.calibration_lamp)

        return ccd

    def _automatic_wavelength_solution(self, corr_tolerance=15):
        """Finds a Wavelength Solution Automatically

        This method uses a library of previously wavelength-calibrated
        comparison lamps. It will only process them if they are the exact match.
        A workflow summary is presented below:
          - Identify the exactly matching reference comparison lamp. If it
            doesn't exist it will return None. If it does exist the reference
            lamp will be loaded and it's wavelength solution read.
          - Identify lines in the new lamp, the lamp data has been already
            loaded at the initialization of the class
          - According to the lines detected it will split both spectrum in the
            same number of pieces and same respective sizes and then will do
            cross correlation for each of them.
          - The line's pixel value is stored
          - Using the reference lamp's wavelength solution mathematical model,
            the corresponding value in angstrom is calculated using the offset
            obtained from the cross correlation something like this:
            angstrom = model(pixel + offset)
          - As a first order filter one-iteration of a two-sigma clipping is
            applied to the cross-correlation offsets, this is necessary to
            eliminate mismatched lines.
          - A new wavelength solution is calculated using the points collected
            above.
          - Using the Angstrom values previously found and the detected lines
            plus the newly calculated solution, the differences in angstrom are
            calculated to which values a new sigma-clipping is applied, again
            one iteration two-sigmas, since the distributions are not
            necessarily normal distributions.
          - Once these values are cleaned of rejected values the final solution
            is calculated.

        Returns:
            None in case it is not possible to find a suitable template lamp or
            if is not possible to calculate the solution.

        """
        # TODO (simon): Implement the use of the binning value
        try:
            reference_lamp_ccd = self.reference_data.get_reference_lamp(
                header=self.lamp.header)

            self.log.debug('Found reference lamp: '
                           '{:s}'.format(reference_lamp_ccd.header['GSP_FNAM']))
        except NotImplementedError:

            self.log.warning('This configuration is not supported in '
                             'automatic mode.')

            # TODO (simon): Evaluate if send this to interactive mode
            return None

        reference_lamp_wav_axis, reference_lamp_ccd.data = \
            self.wcs.read_gsp_wcs(ccd=reference_lamp_ccd)

        if self.serial_binning != 1:
            reference_lamp_wav_axis, reference_lamp_ccd.data = \
                self._bin_reference_data(wavelength=reference_lamp_wav_axis,
                                         intensity=reference_lamp_ccd.data)

            self.wcs.binning = self.serial_binning

        '''detect lines in comparison lamp (not reference)'''
        lamp_lines_pixel = self._get_lines_in_lamp()
        lamp_lines_angst = self.wcs.model(lamp_lines_pixel)

        pixel_values = []
        angstrom_values = []
        correlation_values = []
        angstrom_differences = []

        self.log.debug('Length {:d}'.format(len(self.lamp.data)))
        self.log.debug('NLines {:d}'.format(len(lamp_lines_pixel)))

        self.log.debug('Length / NLines {:.3f}'.format(
            len(self.lamp.data) / float(len(lamp_lines_pixel))))

        global_cross_corr = self._cross_correlation(reference_lamp_ccd.data,
                                                    self.lamp.data)

        half_width = np.max(
            [int((len(self.lamp.data) / float(len(lamp_lines_pixel)))),
             4 * global_cross_corr])

        for i in range(len(lamp_lines_pixel)):
            line_value_pixel = lamp_lines_pixel[i]
            line_value_angst = lamp_lines_angst[i]

            xmin = int(max(0, round(line_value_pixel - half_width)))

            xmax = int(min(round(line_value_pixel + half_width),
                           len(self.lamp.data)))

            if xmin >= xmax:
                continue
            # print(xmin, xmax, self.lamp.data.size)
            # TODO (simon): Convolve to match wider lines such as those from
            # TODO (cont): the slit of 5 arseconds
            ref_sample = reference_lamp_ccd.data[xmin:xmax]
            # ref_wavele = reference_lamp_wav_axis[xmin:xmax]
            lamp_sample = self.lamp.data[xmin:xmax]

            correlation_value = self._cross_correlation(ref_sample, lamp_sample)

            self.log.debug('Cross correlation value '
                           '{:s} vs {:s}'.format(str(global_cross_corr),
                                                 str(correlation_value)))

            if - corr_tolerance < (global_cross_corr - correlation_value) < \
                    corr_tolerance:
                """record value for reference wavelength"""
                # print(global_cross_corr - correlation_value)
                angstrom_value_model = self.wcs.model(
                    line_value_pixel + correlation_value)

                # print(correlation_value, angstrom_value_model)
                correlation_values.append(correlation_value)
                angstrom_differences.append(angstrom_value_model -
                                            line_value_angst)
                angstrom_values.append(angstrom_value_model)
                # print(angstrom_values)
                pixel_values.append(line_value_pixel)
            else:
                self.log.debug("Local cross correlation value {:.3f} is too far"
                               " from global cross correlation value "
                               "{:.3f}".format(correlation_value,
                                               global_cross_corr))

            if False:
                # print(global_cross_corr, correlation_value)
                plt.ion()
                plt.title('Samples after cross correlation\n Shift {:.3f}'
                          ''.format(correlation_value))
                plt.xlabel('Pixel Axis')
                plt.ylabel('Intensity')

                plt.plot(ref_sample,
                         color='k',
                         label='Reference Sample')

                plt.plot([x + correlation_value for x in
                          range(len(lamp_sample))],
                         lamp_sample,
                         label='New Lamp Sample')

                plt.legend(loc='best')
                plt.draw()
                plt.pause(1)
                plt.clf()
                plt.ioff()

        # This is good and necessary as a first approach for some very wrong
        # correlation results
        clipped_values = sigma_clip(correlation_values,
                                    sigma=3,
                                    iters=1,
                                    cenfunc=np.ma.median)
        # print(clipped_values)

        if np.ma.is_masked(clipped_values):
            _pixel_values = list(pixel_values)
            _angstrom_values = list(angstrom_values)
            # print(_angstrom_values)
            pixel_values = []
            angstrom_values = []
            for i in range(len(clipped_values)):
                if clipped_values[i] is not np.ma.masked:
                    pixel_values.append(_pixel_values[i])
                    # print(_angstrom_values[i][0])
                    angstrom_values.append(_angstrom_values[i])

        # Create a wavelength solution
        self.log.info('Creating Wavelength Solution')

        self.wsolution = self.wcs.fit(physical=pixel_values,
                                      wavelength=angstrom_values,
                                      model_name='chebyshev',
                                      degree=self.poly_order)

        if self.wsolution is None:
            self.log.error('Failed to find wavelength solution using reference '
                           'file: {:s}'.format(self.calibration_lamp))
            return None

        # finding differences in order to improve the wavelength solution
        wavelength_differences = [angstrom_values[i] -
                                  self.wsolution(pixel_values[i]) for i in
                                  range(len(pixel_values))]

        clipped_differences = sigma_clip(wavelength_differences,
                                         sigma=2,
                                         iters=3,
                                         cenfunc=np.ma.median)

        if np.ma.is_masked(clipped_differences):

            self.log.debug('Cleaning pixel to angstrom match to improve '
                           'wavelength solution')

            _pixel_values = list(pixel_values)
            _angstrom_values = list(angstrom_values)
            pixel_values = []
            angstrom_values = []
            for i in range(len(clipped_differences)):
                if clipped_differences[i] is not np.ma.masked:
                    pixel_values.append(_pixel_values[i])
                    angstrom_values.append(_angstrom_values[i])
            self.log.info('Re-fitting wavelength solution')

            self.wsolution = self.wcs.fit(physical=pixel_values,
                                          wavelength=angstrom_values,
                                          model_name='chebyshev',
                                          degree=self.poly_order)

        self._evaluate_solution(clipped_differences=clipped_differences)

        if self.args.plot_results or self.args.debug_with_plots or \
                self.args.save_plots:  # pragma: no cover
            plt.close('all')
            plt.switch_backend('Qt5Agg')
            # print(self.i_fig)
            self.i_fig = None
            if self.i_fig is None:
                self.i_fig = plt.figure()
                self.i_fig.canvas.set_window_title(
                    'Automatic Wavelength Solution')
                self.ax1 = self.i_fig.add_subplot(111)
                self.ax1.set_rasterization_zorder(1)

                mng = plt.get_current_fig_manager()
                mng.window.showMaximized()
            # else:
            #     print("clear figure")
            #     self.i_fig.clf()
            #     self.i_fig.canvas.set_window_title(
            #         'Blank')
            if not self.args.debug_with_plots:
                plt.ion()
                # plt.show()
            else:
                plt.ioff()

            self.ax1.plot([], color='m', label='Pixels')
            self.ax1.plot([], color='c', label='Angstrom')
            for val in pixel_values:
                self.ax1.axvline(self.wsolution(val), color='m', zorder=0)
            for val2 in angstrom_values:
                self.ax1.axvline(val2, color='c', linestyle='--', zorder=0)

            self.ax1.plot(reference_lamp_wav_axis,
                          reference_lamp_ccd.data,
                          label='Reference',
                          color='k',
                          alpha=1, zorder=0)

            self.ax1.plot(self.wsolution(self.raw_pixel_axis),
                          self.lamp.data,
                          label='Last Solution',
                          color='r',
                          alpha=0.7, zorder=0)

            try:
                wavmode = self.lamp.header['wavmode']
            except KeyError as error:
                self.log.debug(error)
                wavmode = ''

            self.ax1.set_xlabel('Wavelength (Angstrom)')
            self.ax1.set_ylabel('Intensity (ADU)')

            self.ax1.set_title('Automatic Wavelength Solution\n'
                               + self.lamp.header['OBJECT']
                               + ' ' + wavmode + '\n'
                               + 'RMS Error: {:.3f}'.format(self.rms_error))

            self.ax1.legend(loc='best')
            self.i_fig.tight_layout()

            if self.args.save_plots:

                plots_path = os.path.join(self.args.destination, 'plots')
                if not os.path.isdir(plots_path):
                    os.path.os.makedirs(plots_path)
                # saves pdf files of the wavelength solution plot
                out_file_name = 'automatic-solution_' + self.lamp.header[
                    'GSP_FNAM']
                out_file_name = re.sub('.fits', '', out_file_name)

                file_count = len(glob.glob(
                    os.path.join(self.args.destination,
                                 out_file_name + '*'))) + 1

                out_file_name += '_RMS_{:.3f}_{:03d}.pdf'.format(self.rms_error,
                                                                 file_count)
                pdf_pages = PdfPages(
                    os.path.join(plots_path, out_file_name))
                plt.savefig(pdf_pages, format='pdf')
                pdf_pages.close()

                plot_name = os.path.join(plots_path,
                                         re.sub('pdf', 'png', out_file_name))

                plt.savefig(plot_name, rasterized=True, format='png', dpi=300)

                plt.ioff()
                plt.clf()
            if self.args.debug_with_plots or self.args.plot_results:  # pragma: no cover

                manager = plt.get_current_fig_manager()

                if plt.get_backend() == u'GTK3Agg':
                    manager.window.maximize()
                elif plt.get_backend() == u'Qt5Agg':
                    manager.window.showMaximized()

                if self.args.debug_with_plots:
                    plt.show()
                elif self.args.plot_results:
                    plt.draw()
                    plt.pause(1)
                    plt.ioff()
                    plt.close()
                    # plt.close(self.i_fig)
            # else:
            #     plt.close('all')

    def _bin_reference_data(self, wavelength, intensity):
        """Bins a 1D array

        This method reduces the size of an unbinned array by binning.
        The function to combine data is `numpy.mean`.

        Args:
            wavelength (array): Wavelength axis
            intensity (array): Intensity

        Returns:
            Binned wavelength and intensity arrays.

        """
        if self.serial_binning != 1:
            b_wavelength = ccdproc.block_reduce(wavelength,
                                                self.serial_binning,
                                                np.mean)
            b_intensity = ccdproc.block_reduce(intensity,
                                               self.serial_binning,
                                               np.mean)
            return b_wavelength, b_intensity
        else:
            return wavelength, intensity

    def _cross_correlation(self, reference, new_array, mode='full', plot=False):
        """Do cross correlation of two arrays

        Args:
            reference (array): Reference array.
            new_array (array): Array to be matched.
            mode (str): Correlation mode for `scipy.signal.correlate`.
            plot (bool): Switch debugging plots on or off.

        Returns:
            correlation_value (int): Shift value in pixels.

        """
        cyaxis2 = new_array
        if float(re.sub('[A-Za-z" ]', '', self.lamp.header['SLIT'])) > 3:

            box_width = float(
                re.sub('[A-Za-z" ]',
                       '',
                       self.lamp.header['SLIT'])) / (0.15 * self.serial_binning)

            self.log.debug('BOX WIDTH: {:f}'.format(box_width))
            box_kernel = Box1DKernel(width=box_width)
            max_before = np.max(reference)
            cyaxis1 = convolve(reference, box_kernel)
            max_after = np.max(cyaxis1)
            cyaxis1 *= max_before / max_after

        else:
            kernel_stddev = float(
                re.sub('[A-Za-z" ]',
                       '',
                       self.lamp.header['SLIT'])) / (0.15 * self.serial_binning)

            gaussian_kernel = Gaussian1DKernel(stddev=kernel_stddev)
            cyaxis1 = convolve(reference, gaussian_kernel)
            cyaxis2 = convolve(new_array, gaussian_kernel)

        ccorr = signal.correlate(cyaxis1, cyaxis2, mode=mode)

        max_index = np.argmax(ccorr)

        x_ccorr = np.linspace(-int(len(ccorr) / 2.),
                              int(len(ccorr) / 2.),
                              len(ccorr))

        correlation_value = x_ccorr[max_index]
        if plot:
            plt.ion()
            plt.title('Cross Correlation')
            plt.xlabel('Lag Value')
            plt.ylabel('Correlation Value')
            plt.plot(x_ccorr, ccorr)
            plt.draw()
            plt.pause(2)
            plt.clf()
            plt.ioff()
        return correlation_value

    def _evaluate_solution(self, clipped_differences):
        """Calculates Root Mean Square Error for the wavelength solution.

        Args:
            clipped_differences (ndarray): Numpy masked array of differences
              between reference line values in angstrom and the value calculated
              using the model of the wavelength solution.

        Returns:
            Root Mean Square Error, number of points and number of points
              rejected in the calculation of the wavelength solution.

        """
        self.n_points = len(clipped_differences)
        self.n_rejections = np.ma.count_masked(clipped_differences)
        square_differences = []
        for i in range(len(clipped_differences)):
            if clipped_differences[i] is not np.ma.masked:
                square_differences.append(clipped_differences[i] ** 2)
        self.rms_error = np.sqrt(
            np.sum(square_differences) / len(square_differences))

        self.log.info('Wavelength solution RMS Error : {:.3f}'.format(
            self.rms_error))

        return self.rms_error, self.n_points, self.n_rejections

    def _get_lines_in_lamp(self, ccddata_lamp=None):
        """Identify peaks in a lamp spectrum

        Uses `scipy.signal.argrelmax` to find peaks in a spectrum i.e emission
        lines, then it calls the _recenter_lines method that will recenter them
        using a "center of mass", because, not always the maximum value (peak)
        is the center of the line.

        Returns:
            lines_candidates (list): A common list containing pixel values at
                approximate location of lines.

        """

        if ccddata_lamp is None:
            lamp_data = self.lamp.data
            lamp_header = self.lamp.header
            raw_pixel_axis = self.raw_pixel_axis
        elif isinstance(ccddata_lamp, CCDData):
            # print(ccddata_lamp.data.shape)
            lamp_data = ccddata_lamp.data
            lamp_header = ccddata_lamp.header
            raw_pixel_axis = range(len(lamp_data))
        else:
            self.log.error('Error receiving lamp')
            return None

        no_nan_lamp_data = np.asarray(np.nan_to_num(lamp_data))

        filtered_data = np.where(
            np.abs(no_nan_lamp_data > no_nan_lamp_data.min() +
                   0.03 * no_nan_lamp_data.max()),
            no_nan_lamp_data,
            None)

        # replace None to zero and convert it to an array
        none_to_zero = [0 if it is None else it for it in filtered_data]
        filtered_data = np.array(none_to_zero)

        _upper_limit = no_nan_lamp_data.min() + 0.03 * no_nan_lamp_data.max()
        slit_size = np.float(re.sub('[a-zA-Z" ]', '', lamp_header['slit']))

        serial_binning, parallel_binning = [
            int(x) for x in lamp_header['CCDSUM'].split()]

        new_order = int(round(float(slit_size) / (0.15 * serial_binning)))
        self.log.debug('New Order:  {:d}'.format(new_order))

        # print(round(new_order))
        peaks = signal.argrelmax(filtered_data, axis=0, order=new_order)[0]

        if slit_size >= 5.:

            lines_center = self._recenter_broad_lines(
                lamp_data=no_nan_lamp_data,
                lines=peaks,
                order=new_order)
        else:
            # lines_center = peaks
            lines_center = self._recenter_lines(no_nan_lamp_data, peaks)

        if self.args.debug_with_plots:  # pragma: no cover
            # print(new_order, slit_size, )
            plt.close('all')
            fig, ax = plt.subplots()
            # ax = fig.add_subplot(111)
            fig.canvas.set_window_title('Lines Detected')

            mng = plt.get_current_fig_manager()
            mng.window.showMaximized()

            ax.set_title('Lines detected in Lamp\n'
                         '{:s}'.format(lamp_header['OBJECT']))
            ax.set_xlabel('Pixel Axis')
            ax.set_ylabel('Intensity (counts)')

            # Build legends without data to avoid repetitions
            ax.plot([], color='k', label='Comparison Lamp Data')

            ax.plot([], color='k', linestyle=':',
                    label='Spectral Line Detected')

            ax.axhline(_upper_limit, color='r')

            for line in peaks:
                ax.axvline(line, color='k', linestyle=':')

            ax.plot(raw_pixel_axis, no_nan_lamp_data, color='k')
            ax.legend(loc='best')
            plt.tight_layout()
            plt.show()

        return lines_center

    def _get_spectral_characteristics(self):
        """Calculates some Goodman's specific spectroscopic values.

        From the header value for Grating, Grating Angle and Camera Angle it is
        possible to estimate what are the wavelength values at the edges as well
        as in the center. It was necessary to add offsets though, since the
        formulas provided are slightly off. The values are only an estimate.

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
        # TODO (simon): find a definite solution for this, this only work
        # TODO (simon): (a little) for one configuration
        blue_correction_factor = -50 * u.angstrom
        red_correction_factor = -37 * u.angstrom

        self.grating_frequency = float(re.sub('[A-Za-z_-]',
                                              '',
                                              self.lamp.header['GRATING'])
                                       ) / u.mm

        self.grating_angle = float(self.lamp.header['GRT_ANG']) * u.deg
        self.camera_angle = float(self.lamp.header['CAM_ANG']) * u.deg

        # serial binning - dispersion binning
        # parallel binning - spatial binning
        self.serial_binning, self.parallel_binning = [
            int(x) for x in self.lamp.header['CCDSUM'].split()]

        self.pixel_count = len(self.lamp.data)
        # Calculations
        # TODO (simon): Check whether is necessary to remove the
        # TODO (simon): self.slit_offset variable
        self.alpha = self.grating_angle.to(u.rad)
        self.beta = self.camera_angle.to(u.rad) - self.grating_angle.to(u.rad)

        self.center_wavelength = (np.sin(self.alpha) +
                                  np.sin(self.beta)) / self.grating_frequency

        limit_angle = np.arctan(
            self.pixel_count *
            ((self.pixel_size * self.serial_binning) /
             self.goodman_focal_length) / 2)

        self.blue_limit = (
            (np.sin(self.alpha) + np.sin(self.beta - limit_angle.to(u.rad))) /
            self.grating_frequency).to(u.angstrom) + blue_correction_factor

        self.red_limit = (
            (np.sin(self.alpha) + np.sin(self.beta + limit_angle.to(u.rad))) /
            self.grating_frequency).to(u.angstrom) + red_correction_factor

        pixel_one = 0
        pixel_two = 0
        self.log.debug(
            'Center Wavelength : {:.3f} Blue Limit : '
            '{:.3f} Red Limit : {:.3f} '.format(
                self.center_wavelength.to(u.angstrom),
                self.blue_limit,
                self.red_limit))

        spectral_characteristics = {'center': self.center_wavelength,
                                    'blue': self.blue_limit,
                                    'red': self.red_limit,
                                    'alpha': self.alpha,
                                    'beta': self.beta,
                                    'pix1': pixel_one,
                                    'pix2': pixel_two}
        return spectral_characteristics

    def get_wsolution(self):
        """Returns the mathematical model of the wavelength solution

        The wavelength solution is a callable mathematical function from
        astropy.modeling.models. By obtaining this mathematical model the user
        can use its own method to apply it to a given data.

        Returns:
            A callable mathematical function. None if the wavelength solution
            doesn't exist.

        """
        if self.wsolution is not None:
            return self.wsolution
        else:
            self.log.error("Wavelength Solution doesn't exist!")
            return None

    def _linearize_spectrum(self, data, plots=False):
        """Produces a linearized version of the spectrum

        Storing wavelength solutions in a FITS header is not simple at all for
        non-linear solutions therefore is easier for the final user and for the
        development code to have the spectrum linearized. It first finds a
        spline representation of the data, then creates a linear wavelength axis
        (angstrom) and finally it resamples the data from the spline
        representation to the linear wavelength axis.

        It also applies a median filter of kernel size three to smooth the
        linearized spectrum. Sometimes the splines produce funny things when
        the original data is too steep.

        Args:
            data (Array): The non-linear spectrum
            plots (bool): Whether to show the plots or not

        Returns:
            linear_data (list): Contains two elements: Linear wavelength axis
            and the smoothed linearized data itself.

        """
        # for data_point in data:
        #     print(data_point)
        # print('data ', data)
        pixel_axis = range(len(data))
        if np.nan in data:
            print("there are nans")
            sys.exit(0)

        # print(pixel_axis)
        if self.wsolution is not None:
            x_axis = self.wsolution(pixel_axis)
            try:
                plt.imshow(data)
                plt.show()
            except TypeError:
                pass
            new_x_axis = np.linspace(x_axis[0], x_axis[-1], len(data))
            tck = scipy.interpolate.splrep(x_axis, data, s=0)
            linearized_data = scipy.interpolate.splev(new_x_axis,
                                                      tck,
                                                      der=0)

            smoothed_linearized_data = signal.medfilt(linearized_data)
            # print('sl ', smoothed_linearized_data)
            if plots:  # pragma: no cover
                fig6 = plt.figure(6)
                plt.xlabel('Wavelength (Angstrom)')
                plt.ylabel('Intensity (Counts)')
                fig6.canvas.set_window_title('Linearized Data')

                plt.plot(x_axis,
                         data,
                         color='k',
                         label='Data')

                plt.plot(new_x_axis,
                         linearized_data,
                         color='r',
                         linestyle=':',
                         label='Linearized Data')

                plt.plot(new_x_axis,
                         smoothed_linearized_data,
                         color='m',
                         alpha=0.5,
                         label='Smoothed Linearized Data')

                fig6.tight_layout()
                plt.legend(loc=3)
                plt.show()

                fig7 = plt.figure(7)
                plt.xlabel('Pixels')
                plt.ylabel('Angstroms')
                fig7.canvas.set_window_title('Wavelength Solution')
                plt.plot(x_axis, color='b', label='Non linear wavelength-axis')
                plt.plot(new_x_axis, color='r', label='Linear wavelength-axis')
                fig7.tight_layout()
                plt.legend(loc=3)
                plt.show()

            linear_data = [new_x_axis, smoothed_linearized_data]
            return linear_data

    def _recenter_lines(self, data, lines, plots=False):
        """Finds the centroid of an emission line

        For every line center (pixel value) it will scan left first until the
        data stops decreasing, it assumes it is an emission line and then will
        scan right until it stops decreasing too. Defined those limits it will
        use the line data in between and calculate the centroid.

        Notes:
            This method is used to recenter relatively narrow lines only, there
            is a special method for dealing with broad lines.

        Args:
            data (ndarray): numpy.ndarray instance. or the data attribute of a
                :class:`~astropy.nddata.CCDData` instance.
            lines (list): A line list in pixel values.
            plots (bool): If True will plot spectral line as well as the input
                center and the recentered value.

        Returns:
            A list containing the recentered line positions.

        """
        new_center = []
        x_size = data.shape[0]
        median = np.median(data)
        for line in lines:
            # TODO (simon): Check if this definition is valid, so far is not
            # TODO (cont..): critical
            left_limit = 0
            right_limit = 1
            condition = True
            left_index = int(line)

            while condition and left_index - 2 > 0:

                if (data[left_index - 1] > data[left_index]) and \
                        (data[left_index - 2] > data[left_index - 1]):

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

                if (data[right_index + 1] > data[right_index]) and \
                        (data[right_index + 2] > data[right_index + 1]):

                    condition = False
                    right_limit = right_index

                elif data[right_index] < median:
                    condition = False
                    right_limit = right_index

                else:
                    right_limit = right_index
                right_index += 1
            index_diff = [abs(line - left_index), abs(line - right_index)]

            sub_x_axis = range(line - min(index_diff),
                               (line + min(index_diff)) + 1)

            sub_data = data[line - min(index_diff):(line + min(index_diff)) + 1]
            centroid = np.sum(sub_x_axis * sub_data) / np.sum(sub_data)

            # checks for asymmetries
            differences = [abs(data[line] - data[left_limit]),
                           abs(data[line] - data[right_limit])]

            if max(differences) / min(differences) >= 2.:
                if plots:  # pragma: no cover
                    plt.axvspan(line - 1, line + 1, color='g', alpha=0.3)
                new_center.append(line)
            else:
                new_center.append(centroid)
        if plots:  # pragma: no cover
            fig, ax = plt.subplots(1, 1)
            fig.canvas.set_window_title('Lines Detected in Lamp')
            ax.axhline(median, color='b')

            ax.plot(self.raw_pixel_axis,
                    data,
                    color='k',
                    label='Lamp Data')

            for line in lines:

                ax.axvline(line + 1,
                           color='k',
                           linestyle=':',
                           label='First Detected Center')

            for center in new_center:

                ax.axvline(center,
                           color='k',
                           linestyle='.-',
                           label='New Center')

            plt.show()
        return new_center

    @staticmethod
    def _recenter_broad_lines(lamp_data, lines, order):
        """Recenter broad lines

        Notes:
            This method is used to recenter broad lines only, there is a special
            method for dealing with narrower lines.

        Args:
            lamp_data (ndarray): numpy.ndarray instance. It contains the lamp
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

    def _save_science_data(self, ccd, index=None):
        """Save science data"""
        ccd = ccd.copy()
        linear_x_axis, ccd.data = self._linearize_spectrum(ccd.data)

        ccd = self.wcs.write_gsp_wcs(ccd=ccd,
                                     model=self.wsolution)

        ccd = self.add_wavelength_solution(
            ccd=ccd,
            x_axis=linear_x_axis)

        self._save_wavelength_calibrated(
            ccd=ccd,
            original_filename=ccd.header['GSP_FNAM'],
            index=index)

        # wavelength_solution = WavelengthSolution(
        #     solution_type='non_linear',
        #     model_name='chebyshev',
        #     model_order=self.poly_order,
        #     model=self.wsolution,
        #     ref_lamp=self.calibration_lamp,
        #     eval_comment=self.evaluation_comment,
        #     header=ccd.header)

        if self.args.plot_results or self.args.debug_with_plots or \
                self.args.save_plots:  # pragma: no cover

            plt.close(1)
            if self.args.plot_results:
                plt.ion()
                # plt.show()
            elif self.args.debug_with_plots:
                plt.ioff()

            wavelength_axis = self.wsolution(range(ccd.data.size))

            object_name = ccd.header['OBJECT']
            grating = ccd.header['GRATING']

            fig_title = 'Wavelength Calibrated Data : ' \
                        '{:s}\n{:s}'.format(object_name, grating)

            fig, ax1 = plt.subplots(1)
            fig.canvas.set_window_title(ccd.header['GSP_FNAM'])
            # ax1 = fig.add_subplot(111)

            mng = plt.get_current_fig_manager()
            mng.window.showMaximized()

            ax1.set_title(fig_title)
            ax1.set_xlabel('Wavelength (Angstrom)')
            ax1.set_ylabel('Intensity (ADU)')
            ax1.set_xlim((wavelength_axis[0], wavelength_axis[-1]))
            # plt.close(1)

            ax1.plot(wavelength_axis,
                     ccd.data,
                     color='k',
                     label='Data')

            ax1.legend(loc='best')
            fig.tight_layout()
            if self.args.save_plots:
                self.log.info('Saving plots')
                plots_dir = os.path.join(self.args.destination,
                                         'plots')
                if not os.path.isdir(plots_dir):
                    os.mkdir(plots_dir)
                plot_name = re.sub('.fits',
                                   '.png',
                                   ccd.header['GSP_FNAM'])
                plot_path = os.path.join(plots_dir, plot_name)
                # print(plot_path)
                plt.savefig(plot_path, dpi=300)
                self.log.info('Saved plot as {:s} file '
                              'DPI=300'.format(plot_name))

            if self.args.debug_with_plots or self.args.plot_results:  # pragma: no cover
                manager = plt.get_current_fig_manager()
                if plt.get_backend() == u'GTK3Agg':
                    manager.window.maximize()
                elif plt.get_backend() == u'Qt5Agg':
                    manager.window.showMaximized()

                if self.args.debug_with_plots:
                    plt.show()
                elif self.args.plot_results:
                    plt.draw()
                    plt.pause(2)
                    plt.ioff()

                # return wavelength_solution

    def _save_wavelength_calibrated(self,
                                    ccd,
                                    original_filename,
                                    index=None,
                                    lamp=False):
        if index is None:
            f_end = '.fits'
        else:
            f_end = '_ws_{:d}.fits'.format(index)

        new_filename = os.path.join(self.args.destination,
                                    self.args.output_prefix +
                                    original_filename.replace('.fits', f_end))

        if lamp:
            self.log.info('Wavelength-calibrated {:s} file saved to: '
                          '{:s} for science file {:s}'
                          ''.format(ccd.header['OBSTYPE'],
                                    os.path.basename(new_filename),
                                    self.sci_target_file))

            ccd.header.set('GSP_SCTR',
                           value=self.sci_target_file,
                           after='GSP_FLAT')
        else:
            self.log.info('Wavelength-calibrated {:s} file saved to: '
                          '{:s} using reference lamp {:s}'
                          ''.format(ccd.header['OBSTYPE'],
                                    os.path.basename(new_filename),
                                    self.wcal_lamp_file))
            ccd.header.set(
                'GSP_LAMP',
                value=self.wcal_lamp_file,
                comment='Reference lamp used to obtain wavelength solution',
                after='GSP_FLAT')

        write_fits(ccd=ccd,
                   full_path=new_filename,
                   parent_file=original_filename)

        return os.path.basename(new_filename)


class WavelengthSolution(object):
    """Contains all relevant information of a given wavelength solution

    Stores the mathematical model that allows to convert from pixel to angstrom
    as well as more detailed information regarding the type of solution and the
    quality.

    """
    def __init__(self,
                 solution_type=None,
                 model_name=None,
                 model_order=0,
                 model=None,
                 ref_lamp=None,
                 eval_comment='',
                 header=None):
        """Init method for the WavelengthSolution class

        Args:
            solution_type (str): Type of wavelength solution.
            model_name (str): Mathematical model name.
            model_order (int): Order of the mathematical model in case it is a
                polynomial which in most cases it is.
            model (object): Instance of astropy.modeling.Model, represents the
                transformation from pixel to angstrom.
            ref_lamp (str): File name of reference lamp used to find the
                wavelength solution
            eval_comment (str): Text describing the qualitative evaluation of
                the wavelength solution.
            header (object): Instance of astropy.io.fits.header.Header
        """
        self.log = logging.getLogger(__name__)
        self.dtype_dict = {None: -1,
                           'linear': 0,
                           'log_linear': 1,
                           'non_linear': 2}

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
        self.solution_name = self.set_solution_name(header)

    def set_spectral_features(self, header):
        """Creates dictionary that defines the instrument configuration

        Both Blue and Red Camera produce slightly different FITS headers being
        the red camera the one that provides more precise and better
        information. This method will recognize the camera and create the
        dictionary accordingly.

        Notes:
            As of August 2017 both headers are FITS compliant and contain the
            same keywords.

        Args:
            header (Header): Instance of astropy.io.fits.header.Header.

        Returns:
            A dictionary that contains key information regarding the kind of
            spectroscopic data, basically related to the instrument
            configuration.
            The keywords it reads are: INSTCONF, GRATING, ROI, FILTER, FILTER2,
            SLIT, WAVMODE, CAM_ANG, GRT_ANG.

        """
        # TODO (simon): Use CAM_TARG and GRT_TARG instead of CAM_ANG and GRT_ANG
        if header is None:
            self.log.error('Header has not been parsed')
        else:
            try:
                self.log.debug('{:s} Camera'.format(header['INSTCONF']))

                spectral_dict = {'camera': header['INSTCONF'],
                                 'grating': header['GRATING'],
                                 'roi': header['ROI'],
                                 'filter1': header['FILTER'],
                                 'filter2': header['FILTER2'],
                                 'slit': header['SLIT'],
                                 'instconf': header['INSTCONF'],
                                 'wavmode': header['WAVMODE'],
                                 'cam_ang': header['CAM_ANG'],
                                 'grt_ang': header['GRT_ANG']}

                return spectral_dict

            except KeyError:
                self.log.debug('(Old) Blue Camera')

                spectral_dict = {'camera': 'blue',
                                 'grating': header['GRATING'],
                                 'ccdsum': header['CCDSUM'],
                                 'filter1': header['FILTER'],
                                 'filter2': header['FILTER2'],
                                 'slit': header['SLIT'],
                                 'serial_bin': header['PARAM18'],
                                 'parallel_bin': header['PARAM22'],
                                 'cam_ang': header['CAM_ANG'],
                                 'grt_ang': header['GRT_ANG']}

                return spectral_dict

    def check_compatibility(self, header=None):
        """Checks compatibility of new data

        A wavelength solution is stored as an object (this class). As an
        attribute of this class there is a dictionary that contains critical
        parameters of the spectrum that were used to obtain this solution.
        In order to apply the same solution to another spectrum its header has
        to be parsed and then the parameters are compared in order of
        importance.

        Args:
            header (object): FITS header instance from
                astropy.io.fits.header.Header.

        Returns:
            True if the new data is compatible with the solution or False if its
            not.

        """
        if header is not None:
            new_dict = self.set_spectral_features(header)
            for key in new_dict.keys():
                if self.spectral_dict['camera'] == 'red':
                    if key in ['grating', 'roi', 'instconf', 'wavmode'] and \
                                    new_dict[key] != self.spectral_dict[key]:

                        self.log.info('Keyword: {:s} does not Match'.format(
                            key.upper()))

                        self.log.info('{:s} - Solution: {:s} - New '
                                      'Data: {:s}'.format(
                                           key.upper(),
                                           self.spectral_dict[key],
                                           new_dict[key]))

                        return False

                    elif key in ['cam_ang',  'grt_ang'] and \
                            abs(new_dict[key] - self.spectral_dict[key]) > 1:

                        self.log.debug('Keyword: {:s} Lamp: {:s} Data: '
                                       '{:s}'.format(key,
                                                     self.spectral_dict[key],
                                                     new_dict[key]))

                        self.log.info('Solution belong to a different '
                                      'Instrument Configuration.')

                        return False
                    # else:
                    #     return True
                elif self.spectral_dict['camera'] == 'blue':
                    if key in ['grating',
                               'ccdsum',
                               'serial_bin',
                               'parallel_bin'] and \
                                    new_dict[key] != self.spectral_dict[key]:

                        self.log.debug('Keyword: {:s} does not '
                                       'Match'.format(key.upper()))

                        self.log.info('{:s} - Solution: {:s} - New Data: {:s}',
                                      key.upper(),
                                      self.spectral_dict[key],
                                      new_dict[key])

                        return False

                    elif key in ['cam_ang',
                                 'grt_ang'] and abs(
                                float(new_dict[key]) -
                                float(self.spectral_dict[key])) > 1:

                        self.log.debug('Keyword: {:s} Lamp: {:s} Data: {:s}',
                                       key,
                                       self.spectral_dict[key],
                                       new_dict[key])

                        self.log.info('Solution belong to a different '
                                      'Instrument Configuration.')

                        return False
                    # else:
                    #     return True
            return True
        else:
            self.log.error('Header has not been parsed')
            return False

    def set_solution_name(self, header):
        """Defines a name for the solution

        Using the header's information define a string that could be used as a
        keyword to identify a particular solution in the event that multiple
        solutions or instances of this class are stored somewhere/somehow.

        Args:
            header (Header): FITS header instance from
                astropy.io.fits.header.Header.

        Returns:
            Wavelength solution name.

        """
        name_text = ''
        # Grating part of the flat name
        if header['grating'] == '<NO GRATING>':
            try:
                if header['wavmode'] != 'Imaging':
                    name_text += '_nogrt'
            except KeyError:
                self.log.error('KeyError: Blue Camera')
        else:
            grating = header['grating'].split('_')[1]
            name_text += '_' + grating
            # Mode for the grating part of the flat name
            try:
                mode = header['wavmode'].split(' ')[1]
                name_text += '_' + mode.upper()
            except KeyError:
                self.log.error('KeyError: Blue Camera')
            except IndexError:
                # it means it is Custom mode
                mode = header['wavmode']
                # TODO (simon): Include other modes
                if mode == 'Custom':
                    grating_frequency = int(re.sub('[a-zA-Z-]', '', grating))
                    alpha = float(header['grt_ang'])
                    beta = float(header['cam_ang']) - float(header['grt_ang'])

                    center_wavelength = (1e6 / grating_frequency) * (
                        np.sin(alpha * np.pi / 180.) +
                        np.sin(beta * np.pi / 180.))

                    self.log.debug(center_wavelength)

                    name_text += '_' + \
                                 mode.upper() + \
                                 '_{:d}nm'.format(int(round(center_wavelength)))
                else:
                    # print(mode)
                    self.log.error('WAVMODE: {:s} not supported'.format(mode))

        # First filter wheel part of the flat name
        if header['filter'] != '<NO FILTER>':
            name_text += '_' + header['filter']
        # Second filter wheel part of the flat name
        if header['filter2'] != '<NO FILTER>':
            name_text += '_' + header['filter2']
        name_text += '_' + re.sub('[a-zA-Z" ]', '', header['slit'])

        wsolution_name = 'ws' + name_text

        self.log.debug(wsolution_name)
        return wsolution_name


if __name__ == '__main__':
    sys.exit('This can not be run on its own.')
