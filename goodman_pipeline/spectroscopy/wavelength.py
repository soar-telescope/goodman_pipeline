# -*- coding: utf8 -*-
"""Contains the tools to produce a wavelength solution

This module gets the extracted data to produce a wavelength solution, linearize
the spectrum and write the solution to the image's header following the FITS
standard.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import glob
import json
import logging
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np

from astropy.stats import sigma_clip
from ccdproc import CCDData
from matplotlib.backends.backend_pdf import PdfPages

from ..wcs.wcs import WCS

from ..core import (add_linear_wavelength_solution,
                    bin_reference_data,
                    cross_correlation,
                    evaluate_wavelength_solution,
                    get_lines_in_lamp,
                    linearize_spectrum,
                    write_fits)

from ..core import (ReferenceData, NoMatchFound)

log = logging.getLogger(__name__)


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

    def __init__(self):
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
        self.poly_order = 3
        self.wcs = WCS()
        self.wsolution = None
        self.wcal_lamp_file = None
        self.sci_target_file = None
        self.n_points = None
        self.n_rejections = None
        self.rms_error = None
        self.cross_corr_tolerance = 5
        self.reference_data_dir = None
        self.reference_data = None
        self.calibration_lamp = ''
        self.wcal_lamp_file = ''

        # Instrument configuration and spectral characteristics
        self.serial_binning = None
        self.parallel_binning = None

    def __call__(self,
                 ccd,
                 comp_list,
                 save_data_to,
                 reference_data,
                 object_number=None,
                 corr_tolerance=15,
                 output_prefix='w',
                 plot_results=False,
                 save_plots=False,
                 plots=False,
                 json_output=False):
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
            save_data_to (str): Path to save processed data.
            object_number (int): In case of multiple detections in a single
            image this number will be added as a suffix before `.fits` in
            order to allow for multiple 1D files. Default value is None.
            corr_tolerance (int): `cross_corr_tolerance` stands for cross
            correlation tolerance, in other words, how far the cross
            correlation can be from the global cross correlation. It usually
            increases with the frequency of the grating.
            output_prefix (str):  Prefix to add to files.
            plot_results (bool): Present a plot showing the wavelength
            calibrated data.
            save_plots (bool): Save any plot shown. They are saved under
            `<path>/<save_data_to>/plots/` where `<path>/<save_data_to>` is
            the full path to the folder that `save_data_to` is pointing.
            plots (bool): Show plots during operation.


        Returns:
            wavelength_solution (object): The mathematical model of the
                wavelength solution. If it fails to create it will return a
                None element.

        """
        assert isinstance(ccd, CCDData)
        assert isinstance(comp_list, list)

        json_payload = {'wavelength_solution': [],
                        'warning': '',
                        'error': ''}

        if os.path.isdir(reference_data):
            if self.reference_data_dir != reference_data:
                self.reference_data_dir = reference_data
                self.reference_data = ReferenceData(
                    reference_dir=self.reference_data_dir)

        self.cross_corr_tolerance = corr_tolerance
        self.sci_target_file = ccd.header['GSP_FNAM']

        self.i_fig = None

        log.info('Processing Science Target: '
                 '{:s}'.format(ccd.header['OBJECT']))

        if len(comp_list) == 0:
            log.warning("No comparison lamps were provided for file {}"
                        "".format(self.sci_target_file))
            log.error("Ending processing of {}".format(self.sci_target_file))
            if json_output:
                json_payload['error'] ='Unable to process without reference lamps'
                return json_payload
            else:
                return
        else:
            wavelength_solutions = []
            reference_lamp_names = []
            for self.lamp in comp_list:
                self.calibration_lamp = self.lamp.header['GSP_FNAM']
                log.info('Using reference lamp {}'.format(self.calibration_lamp))

                self.raw_pixel_axis = range(self.lamp.shape[0])

                self.lamp_name = self.lamp.header['OBJECT']

                log.info('Processing Comparison Lamp: '
                         '{:s}'.format(self.lamp_name))

                self.lines_center = get_lines_in_lamp(
                    ccd=self.lamp, plots=plots)
                try:
                    self._automatic_wavelength_solution(
                        save_data_to=save_data_to,
                        corr_tolerance=self.cross_corr_tolerance)
                except NoMatchFound as message:
                    raise NoMatchFound(message)

                if self.wsolution is not None:
                    ccd.header.set('GSP_WRMS', value=self.rms_error)
                    ccd.header.set('GSP_WPOI', value=self.n_points)
                    ccd.header.set('GSP_WREJ', value=self.n_rejections)

                    linear_x_axis, self.lamp.data = linearize_spectrum(
                        self.lamp.data,
                        wavelength_solution=self.wsolution)

                    self.lamp = self.wcs.write_gsp_wcs(ccd=self.lamp,
                                                       model=self.wsolution)

                    self.lamp = add_linear_wavelength_solution(
                        ccd=self.lamp,
                        x_axis=linear_x_axis,
                        reference_lamp=self.calibration_lamp)

                    self.wcal_lamp_file = self._save_wavelength_calibrated(
                        ccd=self.lamp,
                        original_filename=self.calibration_lamp,
                        save_data_to=save_data_to,
                        output_prefix=output_prefix,
                        index=object_number,
                        lamp=True)

                    wavelength_solutions.append(self.wsolution)
                    reference_lamp_names.append(self.wcal_lamp_file)
                else:
                    log.error('It was not possible to get a wavelength '
                              'solution from lamp '
                              '{:s} {:s}.'.format(
                               self.lamp.header['GSP_FNAM'],
                               self.lamp.header['OBJECT']))
                    continue

            if len(wavelength_solutions) > 1:
                warning_message = str("The current version of the pipeline "
                                      "does not combine multiple solution "
                                      "instead it saves a single version of "
                                      "the science file for each wavelength "
                                      "solution calculated.")
                log.warning(warning_message)
                all_solution_info = []
                for i in range(len(wavelength_solutions)):
                    # TODO (simon): Combine Multiple solutions
                    self.wsolution = wavelength_solutions[i]
                    self.wcal_lamp_file = reference_lamp_names[i]
                    ccd = self.wcs.write_gsp_wcs(ccd=ccd, model=self.wsolution)
                    saved_file_name = self._save_science_data(
                        ccd=ccd,
                        wavelength_solution=self.wsolution,
                        save_to=save_data_to,
                        index=i + 1,
                        plot_results=plot_results,
                        save_plots=save_plots,
                        plots=plots)
                    all_solution_info.append({
                        'solution_info': {'rms_error': "{:.4f}".format(self.rms_error),
                                          'npoints': "{:d}".format(self.n_points),
                                          'nrjections': "{:d}".format(self.n_rejections)},
                        'file_name': saved_file_name,
                        'reference_lamp': self.wcal_lamp_file})

                if json_output:
                    json_payload['warning'] = warning_message
                    json_payload['wavelength_solution'] = all_solution_info
                    return json_payload

            elif len(wavelength_solutions) == 1:
                self.wsolution = wavelength_solutions[0]
                self.wcal_lamp_file = reference_lamp_names[0]
                ccd = self.wcs.write_gsp_wcs(ccd=ccd, model=self.wsolution)

                saved_file_name = self._save_science_data(
                    ccd=ccd,
                    wavelength_solution=self.wsolution,
                    save_to=save_data_to,
                    plot_results=plot_results,
                    save_plots=save_plots,
                    index=object_number,
                    plots=plots)
                if json_output:
                    json_payload['wavelength_solution'] = [
                        {'solution_info': {'rms_error': "{:.4f}".format(self.rms_error),
                                           'npoints': "{:d}".format(self.n_points),
                                           'nrjections': "{:d}".format(self.n_rejections)},
                         'file_name': saved_file_name,
                         'reference_lamp': self.wcal_lamp_file}]

                    return json_payload
            else:
                log.error("No wavelength solution.")
                if json_output:
                    json_payload['error'] = "Unable to obtain wavelength solution"
                    return json_payload

    def _automatic_wavelength_solution(self,
                                       save_data_to,
                                       corr_tolerance=15,
                                       plot_results=False,
                                       save_plots=False,
                                       plots=False):
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
        try:
            reference_lamp_ccd = self.reference_data.get_reference_lamp(
                header=self.lamp.header)

            log.debug('Found reference lamp: '
                      '{:s}'.format(reference_lamp_ccd.header['GSP_FNAM']))
        except NoMatchFound as error:
            raise NoMatchFound(error)
        except NotImplementedError as error:
            raise NotImplemented(error)

        # TODO (simon): Evaluate possibility to read iraf wcs. [#304]

        reference_lamp_wav_axis, reference_lamp_ccd.data = \
            self.wcs.read_gsp_wcs(ccd=reference_lamp_ccd)

        self.serial_binning, self.parallel_binning = [
            int(x) for x in self.lamp.header['CCDSUM'].split()]

        if self.serial_binning != 1:
            reference_lamp_wav_axis, reference_lamp_ccd.data = \
                bin_reference_data(wavelength=reference_lamp_wav_axis,
                                   intensity=reference_lamp_ccd.data,
                                   serial_binning=self.serial_binning)

            self.wcs.binning = self.serial_binning

        '''detect lines in comparison lamp (not reference)'''
        lamp_lines_pixel = get_lines_in_lamp(ccd=self.lamp,
                                             plots=plots)
        lamp_lines_angst = self.wcs.model(lamp_lines_pixel)

        pixel_values = []
        angstrom_values = []
        correlation_values = []
        angstrom_differences = []

        log.debug('Length {:d}'.format(len(self.lamp.data)))
        log.debug('NLines {:d}'.format(len(lamp_lines_pixel)))

        log.debug('Length / NLines {:.3f}'.format(
            len(self.lamp.data) / float(len(lamp_lines_pixel))))

        slit_size = float(re.sub('["A-Za-z_ ]', '', self.lamp.header['SLIT']))

        global_cross_corr = cross_correlation(
            reference=reference_lamp_ccd.data,
            compared=self.lamp.data,
            slit_size=slit_size,
            serial_binning=self.serial_binning)

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

            slit_size = float(re.sub('["A-Za-z_ ]', '', self.lamp.header['SLIT']))

            correlation_value = cross_correlation(
                reference=ref_sample,
                compared=lamp_sample,
                slit_size=slit_size,
                serial_binning=self.serial_binning)

            log.debug('Cross correlation value '
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
                log.debug("Local cross correlation value {:.3f} is too far "
                          "from {:.3f}".format(correlation_value,
                                               global_cross_corr))

            if plots:  # pragma: no cover
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
                                    maxiters=1,
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
        log.info('Creating Wavelength Solution')

        self.wsolution = self.wcs.fit(physical=pixel_values,
                                      wavelength=angstrom_values,
                                      model_name='chebyshev',
                                      degree=self.poly_order)

        if self.wsolution is None:
            log.error('Failed to find wavelength solution using reference '
                      'file: {:s}'.format(self.calibration_lamp))
            return None

        # finding differences in order to improve the wavelength solution
        wavelength_differences = [angstrom_values[i] -
                                  self.wsolution(pixel_values[i]) for i in
                                  range(len(pixel_values))]

        clipped_differences = sigma_clip(wavelength_differences,
                                         sigma=2,
                                         maxiters=3,
                                         cenfunc=np.ma.median)

        if np.ma.is_masked(clipped_differences):

            log.debug('Cleaning pixel to angstrom match to improve '
                      'wavelength solution')

            _pixel_values = list(pixel_values)
            _angstrom_values = list(angstrom_values)
            pixel_values = []
            angstrom_values = []
            for i in range(len(clipped_differences)):
                if clipped_differences[i] is not np.ma.masked:
                    pixel_values.append(_pixel_values[i])
                    angstrom_values.append(_angstrom_values[i])
            log.info('Re-fitting wavelength solution')

            self.wsolution = self.wcs.fit(physical=pixel_values,
                                          wavelength=angstrom_values,
                                          model_name='chebyshev',
                                          degree=self.poly_order)

        self.rms_error, self.n_points, self.n_rejections = \
            evaluate_wavelength_solution(
                clipped_differences=clipped_differences)

        if plot_results or plots or \
                save_plots:  # pragma: no cover
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

            if not plots:
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
                log.debug(error)
                wavmode = ''

            self.ax1.set_xlabel('Wavelength (Angstrom)')
            self.ax1.set_ylabel('Intensity (ADU)')

            self.ax1.set_title('Automatic Wavelength Solution\n'
                               + self.lamp.header['OBJECT']
                               + ' ' + wavmode + '\n'
                               + 'RMS Error: {:.3f}'.format(self.rms_error))

            self.ax1.legend(loc='best')
            self.i_fig.tight_layout()

            if save_plots:

                plots_path = os.path.join(save_data_to, 'plots')
                if not os.path.isdir(plots_path):
                    os.path.os.makedirs(plots_path)
                # saves pdf files of the wavelength solution plot
                out_file_name = 'automatic-solution_' + self.lamp.header[
                    'GSP_FNAM']
                out_file_name = re.sub('.fits', '', out_file_name)

                file_count = len(glob.glob(
                    os.path.join(save_data_to,
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
            if plots or plot_results:  # pragma: no cover

                manager = plt.get_current_fig_manager()

                if plt.get_backend() == u'GTK3Agg':
                    manager.window.maximize()
                elif plt.get_backend() == u'Qt5Agg':
                    manager.window.showMaximized()

                if plots:
                    plt.show()
                elif plot_results:
                    plt.draw()
                    plt.pause(1)
                    plt.ioff()
                    plt.close()

    def _save_science_data(self,
                           ccd,
                           wavelength_solution,
                           save_to,
                           index=None,
                           plot_results=False,
                           save_plots=False,
                           plots=False):
        """Save wavelength calibrated data

        The spectrum is linearized, then the linear solution is recorded in the
        ccd's header and finally it calls the method
        :func:`~wavelength.WavelengthCalibration._save_wavelength_calibrated`
        which performs the actual saving to a file.

        Args:
            ccd (CCDData): Instance of :class:`~astropy.nddata.CCDData` with a
            1D spectrum.
            wavelength_solution (object): A :class:`~astropy.modeling.Model`
            save_to (str): Path to save location
            index (int): If there are more than one target, they are identified
            by this index.
            plot_results (bool): Whether to show plots or not.
            save_plots (bool): Whether to save plots to files.
            plots

        Returns:
            File name of saved file.

        """
        ccd = ccd.copy()
        linear_x_axis, ccd.data = linearize_spectrum(
            data=ccd.data,
            wavelength_solution=wavelength_solution)

        ccd = add_linear_wavelength_solution(
            ccd=ccd,
            x_axis=linear_x_axis,
            reference_lamp=self.calibration_lamp)

        save_file_name = self._save_wavelength_calibrated(
            ccd=ccd,
            original_filename=ccd.header['GSP_FNAM'],
            save_data_to=save_to,
            index=index)

        if plot_results or plots or save_plots:  # pragma: no cover

            plt.close(1)
            if plot_results:
                plt.ion()
                # plt.show()
            elif plots:
                plt.ioff()

            wavelength_axis = wavelength_solution(range(ccd.data.size))

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
            if save_plots:
                log.info('Saving plots')
                plots_dir = os.path.join(save_to,
                                         'plots')
                if not os.path.isdir(plots_dir):
                    os.mkdir(plots_dir)
                plot_name = re.sub('.fits',
                                   '.png',
                                   ccd.header['GSP_FNAM'])
                plot_path = os.path.join(plots_dir, plot_name)
                # print(plot_path)
                plt.savefig(plot_path, dpi=300)
                log.info('Saved plot as {:s} file '
                              'DPI=300'.format(plot_name))

            if plots or plot_results:  # pragma: no cover
                manager = plt.get_current_fig_manager()
                if plt.get_backend() == u'GTK3Agg':
                    manager.window.maximize()
                elif plt.get_backend() == u'Qt5Agg':
                    manager.window.showMaximized()

                if plots:
                    plt.show()
                elif plot_results:
                    plt.draw()
                    plt.pause(2)
                    plt.ioff()

        return save_file_name

    def _save_wavelength_calibrated(self,
                                    ccd,
                                    original_filename,
                                    save_data_to,
                                    output_prefix='w',
                                    index=None,
                                    lamp=False):
        if index is None:
            f_end = '.fits'
        else:
            f_end = '_ws_{:d}.fits'.format(index)

        file_full_path = os.path.join(save_data_to,
                                    output_prefix +
                                    original_filename.replace('.fits', f_end))

        if lamp:
            log.info('Wavelength-calibrated {:s} file saved to: '
                     '{:s} for science file {:s}'
                     ''.format(ccd.header['OBSTYPE'],
                               os.path.basename(file_full_path),
                               self.sci_target_file))

            ccd.header.set('GSP_SCTR',
                           value=self.sci_target_file,
                           after='GSP_FLAT')
        else:
            log.info('Wavelength-calibrated {:s} file saved to: '
                     '{:s} using reference lamp {:s}'
                     ''.format(ccd.header['OBSTYPE'],
                               os.path.basename(file_full_path),
                               self.wcal_lamp_file))
            ccd.header.set(
                'GSP_LAMP',
                value=self.wcal_lamp_file,
                comment='Reference lamp used to obtain wavelength solution',
                after='GSP_FLAT')

        write_fits(ccd=ccd,
                   full_path=file_full_path,
                   parent_file=original_filename)

        return file_full_path


if __name__ == '__main__':  # pragma: no cover
    sys.exit('This can not be run on its own.')
