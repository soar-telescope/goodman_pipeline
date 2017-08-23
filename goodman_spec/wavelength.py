# -*- coding: utf8 -*-
"""Contains the tools to produce a wavelength solution

This module gets the extracted data to produce a wavelength solution, linearize
the spectrum and write the solution to the image's header following the FITS
standard.
"""

# TODO Reformat file - It is confusing at the moment
# TODO Reformat _ Re-order imports (first "import ...", then "from ... import ..." alphabetically)
# TODO (simon): Discuss this because there are other rules that will probably conflict with this request.
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import astropy.units as u
import glob
import logging
import numpy as np
# import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import os
import re
import sys
import scipy.interpolate

from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
from ccdproc import CCDData
from matplotlib.backends.backend_pdf import PdfPages
from scipy import signal

from .wsbuilder import (ReadWavelengthSolution, WavelengthFitter)
from .linelist import ReferenceData
from goodman_ccd.core import (spectroscopic_extraction,
                              search_comp_group,
                              add_wcs_keys,
                              NoTargetException,
                              NoMatchFound,
                              NotEnoughLinesDetected,
                              CriticalError)



# FORMAT = '%(levelname)s:%(filename)s:%(module)s: 	%(message)s'
# log.basicConfig(level=log.INFO, format=FORMAT)
log = logging.getLogger('redspec.wavelength')

SHOW_PLOTS = False


def process_spectroscopy_data(data_container, args, extraction_type='simple'):
    """Does spectroscopic processing

    This is a high level method that manages all subprocess regarding the
    spectroscopic reduction.

    It returns a True value but it should be modified to return something more
    usable in case it is integrated in other application.

    Args:
        data_container (object): Instance of goodman_ccd.core.NightDataContainer
            class that is used to store classified data.
        args (object): Instance of arparse.Namespace that contains all the
            arguments of the pipeline.
        extraction_type (str): string that defines the type of extraction to be
            performed. 'simple' or 'optimal'. This is required by the extraction
            function.

    Returns:
        a True value.

    """

    assert data_container.is_empty is False
    assert any(extraction_type == option for option in ['simple',
                                                        'optimal'])
    # to be returned
    extracted = None
    comps = None

    full_path = data_container.full_path

    for sub_container in [groups for groups in [data_container.spec_groups,
                          data_container.object_groups] if groups is not None]:
        for group in sub_container:
            # instantiate WavelengthCalibration here for each group.
            get_wsolution = WavelengthCalibration(args=args)
            # this will contain only obstype == OBJECT
            object_group = group[group.obstype == 'OBJECT']
            # this has to be initialized here
            comp_group = None
            comp_ccd_list = []
            if 'COMP' in group.obstype.unique():
                log.debug('Group has comparison lamps')
                comp_group = group[group.obstype == 'COMP']
                log.debug('Adding comparison lamps group to data container')
                data_container.add_comp_group(comp_group=comp_group)
            else:
                log.debug('Group does not have comparison lamps')
                if data_container.comp_groups is not None:
                    log.debug('There are comparison lamp group candidates')

                    try:

                        comp_group = search_comp_group(
                            object_group=object_group,
                            comp_groups=data_container.comp_groups)

                        log.warning('This comparison lamp might not be optimal '
                                    'if you are doing radial velocity studies')

                    except NoMatchFound:

                        log.error('It was not possible to find a comparison '
                                  'group')
                else:
                    log.warning('Data will be extracted but not calibrated')

            for spec_file in object_group.file.tolist():
                log.info('Processing Science File: {:s}'.format(spec_file))
                file_path = os.path.join(full_path, spec_file)
                ccd = CCDData.read(file_path, unit=u.adu)
                ccd.header = add_wcs_keys(header=ccd.header)
                ccd.header['OFNAME'] = (spec_file, 'Original File Name')
                if comp_group is not None and comp_ccd_list == []:
                    for comp_file in comp_group.file.tolist():
                        comp_path = os.path.join(full_path, comp_file)
                        comp_ccd = CCDData.read(comp_path, unit=u.adu)
                        comp_ccd.header = add_wcs_keys(header=comp_ccd.header)
                        comp_ccd.header['OFNAME'] = (comp_file,
                                                     'Original File Name')
                        comp_ccd_list.append(comp_ccd)
                        # plt.imshow(comp_ccd.data)
                        # plt.show()
                else:
                    log.debug('Comp Group is None or comp list already exist')

                try:
                    extracted, comps = spectroscopic_extraction(
                        ccd=ccd,
                        extraction=extraction_type,
                        comp_list=comp_ccd_list,
                        nfind=args.max_n_targets,
                        plots=SHOW_PLOTS)

                    if args.debug_mode:
                        fig = plt.figure(0)
                        fig.clf()
                        fig.canvas.set_window_title('Extracted Data')

                        manager = plt.get_current_fig_manager()

                        if plt.get_backend() == u'GTK3Agg':
                            manager.window.maximize()
                        elif plt.get_backend() == u'Qt4Agg':
                            manager.window.showMaximized()

                        for edata in extracted:
                            plt.plot(edata.data, label=edata.header['OBJECT'])
                            if comps != []:
                                for comp in comps:
                                    plt.plot(comp.data,
                                             label=comp.header['OBJECT'])
                        plt.legend(loc='best')
                        if plt.isinteractive():
                            plt.draw()
                            plt.pause(1)
                        else:
                            plt.show()

                    object_number = None
                    for i in range(len(extracted)):
                        extracted_ccd = extracted[i]

                        # this is for numbering the last file.
                        if len(extracted) == 1:
                            object_number = None
                        else:
                            object_number = i + 1

                        wsolution_obj = get_wsolution(ccd=extracted_ccd,
                                                      comp_list=comps,
                                                      object_number=object_number)
                        # return wsolution_obj
                except NoTargetException:
                    log.error('No target was identified')
                    break
    return True


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
            args (object): Runtime arguments.

        """

        # TODO - Documentation missing
        self.args = args
        self.poly_order = 2
        self.wsolution = None
        self.rms_error = None
        self.reference_data = ReferenceData(self.args.reference_dir)
        # self.science_object = science_object
        self.slit_offset = None
        self.interpolation_size = 200
        self.line_search_method = 'derivative'
        """Instrument configuration and spectral characteristics"""
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
        """Interactive wavelength finding"""
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
        self.filling_value = 1000
        self.events = True
        self.first = True
        self.evaluation_comment = None
        # self.binning = self.lamp_header[]
        self.pixelcenter = []
        """this data must come parsed"""
        self.path = self.args.source
        # self.science_pack = sci_pack
        # self.sci_filename = self.science_object.file_name
        # self.history_of_lamps_solutions = {}
        self.reference_solution = None

    def __call__(self, ccd, comp_list, object_number=None, wsolution_obj=None):
        """Call method for the WavelengthSolution Class

        It takes extracted data and produces wavelength calibrated 1D FITS file.
        The call method takes care of the order and logic needed to call the
        different methods. A wavelength solution can be recycled for the next
        science object. In that case, the wavelength solution is parsed as an
        argument and then there is no need to calculate it again. The recycling
        part has to be implemented in the caller function.

        Args:
            ccd (object): a ccdproc.CCDData instance
            comp_list (list): Comparison lamps for the science target that will
                be processed here. Every element of this list is an instance of
                ccdproc.CCDData.
            object_number (int): In case of multiple detections in a single
                image this number will be added as a suffix before `.fits` in
                order to allow for multiple 1D files. Default value is None.
            wsolution_obj (object): Mathematical model of the wavelength
                solution if exist. If it doesnt is a None

        Returns:
            wavelength_solution (object): The mathematical model of the
                wavelength solution. If it fails to create it will return a
                None element.

        """
        assert isinstance(ccd, CCDData)
        assert isinstance(comp_list, list)

        self.i_fig = None

        log.info('Processing Science Target: {:s}'.format(ccd.header['OBJECT']))
        if comp_list is not None:
            for lamp_ccd in comp_list:

                try:
                    self.calibration_lamp = lamp_ccd.header['OFNAME']
                except KeyError:
                    self.calibration_lamp = ''

                self.lamp_data = lamp_ccd.data
                self.raw_pixel_axis = range(len(self.lamp_data))
                # self.raw_pixel_axis = range(len(self.lamp_data))
                self.lamp_header = lamp_ccd.header.copy()
                self.lamp_name = self.lamp_header['OBJECT']

                log.info('Processing Comparison Lamp: '
                         '{:s}'.format(self.lamp_name))

                self.data1 = self.interpolate(self.lamp_data)
                # self.lines_limits = self.get_line_limits()
                # self.lines_center = self.get_line_centers(self.lines_limits)
                self.lines_center = self.get_lines_in_lamp()
                self.spectral = self.get_spectral_characteristics()
                object_name = ccd.header['OBJECT']
                if self.args.interactive_ws:
                    self.interactive_wavelength_solution(object_name=object_name)
                else:
                    log.warning('Automatic Wavelength Solution might fail to '
                                'provide accurate solutions')
                    self.automatic_wavelength_solution()
                    # self.wsolution = self.wavelength_solution()
                if self.wsolution is not None:
                    # TODO (simon): plug in a record system
                    # record = '{:s} {:.3f} {:.3f}'.format(
                    #     self.lamp_header['GRATING'],
                    #     self.lamp_header['GRT_TARG'],
                    #     self.lamp_header['CAM_TARG'])
                    #
                    # for par in self.wsolution.parameters:
                    #     record += ' {:.5f}'.format(par)
                    #
                    # os.system("echo \'{:s}\' >> parametros.txt".format(record))

                    self.linear_lamp = self.linearize_spectrum(self.lamp_data)

                    self.lamp_header = self.add_wavelength_solution(
                        new_header=self.lamp_header,
                        spectrum=self.linear_lamp,
                        original_filename=self.calibration_lamp,
                        index=object_number)

                    # print(ccd.header)

                    self.linearized_sci = self.linearize_spectrum(ccd.data)

                    self.header = self.add_wavelength_solution(
                        new_header=ccd.header,
                        spectrum=self.linearized_sci,
                        original_filename=ccd.header['OFNAME'],
                        index=object_number)


                    wavelength_solution = WavelengthSolution(
                        solution_type='non_linear',
                        model_name='chebyshev',
                        model_order=self.poly_order,
                        model=self.wsolution,
                        ref_lamp=self.calibration_lamp,
                        eval_comment=self.evaluation_comment,
                        header=self.header)

                    if self.args.plot_results or self.args.debug_mode or \
                            self.args.save_plots:

                        if not self.args.debug_mode:
                            plt.ion()
                            # plt.show()
                        else:
                            plt.ioff()

                        wavelength_axis = self.wsolution(range(ccd.data.size))

                        object_name = ccd.header['OBJECT']
                        grating = ccd.header['GRATING']

                        fig_title = 'Wavelength Calibrated Data : ' \
                                    '{:s}\n{:s}'.format(object_name, grating)

                        fig = plt.figure()
                        fig.canvas.set_window_title(ccd.header['OFNAME'])
                        ax1 = fig.add_subplot(111)
                        manager = plt.get_current_fig_manager()
                        if plt.get_backend() == u'GTK3Agg':
                            manager.window.maximize()
                        elif plt.get_backend() == u'Qt4Agg':
                            manager.window.showMaximized()

                        ax1.set_title(fig_title)
                        ax1.set_xlabel('Wavelength (Angstrom)')
                        ax1.set_ylabel('Intensity (ADU)')
                        ax1.set_xlim((wavelength_axis[0], wavelength_axis[-1]))

                        ax1.plot(wavelength_axis,
                                 ccd.data,
                                 color='k',
                                 label='Data')

                        ax1.legend(loc='best')
                        plt.tight_layout()
                        if self.args.save_plots:
                            log.info('Saving plots')
                            plots_dir = os.path.join(self.args.destiny, 'plots')
                            if not os.path.isdir(plots_dir):
                                os.mkdir(plots_dir)
                            plot_name = re.sub('.fits',
                                               '.png',
                                               ccd.header['OFNAME'])
                            plot_path = os.path.join(plots_dir, plot_name)
                            # print(plot_path)
                            plt.savefig(plot_path, dpi=300)
                            log.info('Saved plot as {:s} file '
                                     'DPI=300'.format(plot_name))

                        if self.args.debug_mode:
                            plt.show()
                        else:
                            plt.draw()
                            plt.pause(2)
                            plt.ioff()

                    return wavelength_solution
                else:
                    log.error('It was not possible to get a wavelength '
                              'solution from this lamp.')
                    return None
        else:
            print('Data should be saved anyways')

        # else:
        #     self.wsolution = wsolution_obj.wsolution
        #     self.calibration_lamp = wsolution_obj.reference_lamp
        #     self.evaluation_comment = wsolution_obj.evaluation_comment
        #     # print('wavelengthSolution ', self.wsolution)
        #     # print('Evaluation Comment', self.evaluation_comment)
        #     # repeat for all sci
        #     for target_index in range(self.science_object.no_targets):
        #         log.debug('Processing target {:d}'.format(target_index + 1))
        #         new_data = self.science_pack.data[target_index]
        #         new_header = self.science_pack.headers[target_index]
        #         if self.science_object.no_targets > 1:
        #             new_index = target_index + 1
        #         else:
        #             new_index = None
        #         self.linearized_sci = self.linearize_spectrum(new_data)
        #         self.header = self.add_wavelength_solution(new_header,
        #                                                    self.linearized_sci,
        #                                                    self.sci_filename,
        #                                                    self.evaluation_comment,
        #                                                    index=new_index)
        #     return wsolution_obj

    def get_wsolution(self):
        """Get the mathematical model of the wavelength solution

        The wavelength solution is a callable mathematical function from
        astropy.modeling.models. By obtaining this mathematical model the user
        can use its own method to apply it to a given data.

        Returns:
            wsolution (callable): A callable mathematical function. None if the
                wavelength solution doesn't exist.

        """
        if self.wsolution is not None:
            return self.wsolution
        else:
            log.error("Wavelength Solution doesn't exist!")
            return None

    def get_calibration_lamp(self):
        """Get the name of the calibration lamp used to obtain the solution

        Returns:
            calibration_lamp (str): Filename of calibration lamp used to obtain
                wavelength solution

        """
        if self.wsolution is not None and self.calibration_lamp is not None:
            return self.calibration_lamp
        else:
            log.error('Wavelength solution has not been calculated yet.')

    def get_lines_in_lamp(self, ccddata_lamp=None):
        """Identify peaks in a lamp spectrum

        Uses scipy.signal.argrelmax to find peaks in a spectrum i.e emission
        lines, then it calls the recenter_lines method that will recenter them
        using a "center of mass", because, not always the maximum value (peak)
        is the center of the line.

        Returns:
            lines_candidates (list): A common list containing pixel values at
                approximate location of lines.

        """

        if ccddata_lamp is None:
            lamp_data = self.lamp_data
            lamp_header = self.lamp_header
            raw_pixel_axis = self.raw_pixel_axis
        elif isinstance(ccddata_lamp, CCDData):
            # print(ccddata_lamp.data.shape)
            lamp_data = ccddata_lamp.data
            lamp_header = ccddata_lamp.header
            raw_pixel_axis = range(len(lamp_data))
        else:
            log.error('Error receiving lamp')

        no_nan_lamp_data = np.nan_to_num(lamp_data)

        filtered_data = np.where(
            np.abs(no_nan_lamp_data > no_nan_lamp_data.min() +
                   0.03 * no_nan_lamp_data.max()),
            no_nan_lamp_data,
            None)

        _upper_limit = no_nan_lamp_data.min() + 0.03 * no_nan_lamp_data.max()
        slit_size = re.sub('[a-zA-Z" ]', '', lamp_header['slit'])

        serial_binning, parallel_binning = [
            int(x) for x in lamp_header['CCDSUM'].split()]

        new_order = int(round(float(slit_size) / (0.15 * serial_binning)))
        log.debug('New Order:  {:d}'.format(new_order))

        # print(round(new_order))
        peaks = signal.argrelmax(filtered_data, axis=0, order=new_order)[0]

        if slit_size >= 5:

            lines_center = self.recenter_broad_lines(lamp_data=no_nan_lamp_data,
                                                     lines=peaks,
                                                     slit_size=slit_size,
                                                     order=new_order)
        else:
            # lines_center = peaks
            lines_center = self.recenter_lines(no_nan_lamp_data, peaks)

        if self.args.debug_mode:
            fig = plt.figure(1)
            fig.canvas.set_window_title('Lines Detected')
            plt.title('Lines detected in Lamp\n{:s}'.format(lamp_header['OBJECT']))
            plt.xlabel('Pixel Axis')
            plt.ylabel('Intensity (counts)')

            # Build legends without data to avoid repetitions
            plt.plot([], color='k', label='Comparison Lamp Data')
            plt.plot([], color='k', linestyle=':', label='Spectral Line Detected')
            plt.axhline(_upper_limit, color='r')

            for line in peaks:
                plt.axvline(line + 1, color='k', linestyle=':')

            # plt.axhline(median + stddev, color='g')
            # for rc_line in lines_center:
            #     plt.axvline(rc_line, color='r')

            plt.plot(raw_pixel_axis, no_nan_lamp_data, color='k')
            plt.legend(loc='best')
            plt.show()

        return lines_center

    def recenter_lines(self, data, lines, plots=False):
        """Finds the centroid of an emission line

        For every line center (pixel value) it will scan left first until the
        data stops decreasing, it assumes it is an emission line and then will
        scan right until it stops decreasing too. Defined those limits it will
        use the line data in between and calculate the centroid.

        Notes:
            This method is used to recenter relatively narrow lines only, there
            is a special method for dealing with broad lines.

        Args:
            data (array): numpy.ndarray instance. or the data attribute of a
                ccdproc.CCDData instance.
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
                if plots:
                    plt.axvspan(line - 1, line + 1, color='g', alpha=0.3)
                new_center.append(line + 1)
            else:
                new_center.append(centroid + 1)
        if plots:
            fig = plt.figure(1)
            fig.canvas.set_window_title('Lines Detected in Lamp')
            plt.axhline(median, color='b')

            plt.plot(self.raw_pixel_axis,
                     data,
                     color='k',
                     label='Lamp Data')

            for line in lines:

                plt.axvline(line + 1,
                            color='k',
                            linestyle=':',
                            label='First Detected Center')

            for center in new_center:

                plt.axvline(center,
                            color='k',
                            linestyle='.-',
                            label='New Center')

            plt.show()
        return new_center

    def recenter_broad_lines(self, lamp_data, lines, order):
        """Recenter broad lines

        Notes:
            This method is used to recenter broad lines only, there is a special
            method for dealing with narrower lines.

        Args:
            lamp_data (array): numpy.ndarray instance. It contains the lamp
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
            # if self.args.debug_mode:
            #     plt.plot(x_axis, lamp_sample)
            #     plt.plot(x_axis, gaussian_model(x_axis))
            #     plt.plot(x_axis, fitted_gaussian(x_axis), color='k')
            #     plt.axvline(line)
            #     plt.show()
        return new_line_centers

    def get_spectral_characteristics(self):
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

        self.grating_frequency = float(
            re.sub('[A-Za-z_-]',
                   '',
                   self.lamp_header['GRATING'])) / u.mm

        # print('Grating Frequency ' + '{:d}'.format(int(self.grating_frequency)))
        self.grating_angle = float(self.lamp_header['GRT_ANG']) * u.deg
        self.camera_angle = float(self.lamp_header['CAM_ANG']) * u.deg

        self.serial_binning, self.parallel_binning = [
            int(x) for x in self.lamp_header['CCDSUM'].split()]

        self.pixel_count = len(self.lamp_data)
        # Calculations
        # TODO (simon): Check whether is necessary to remove the self.slit_offset variable
        self.alpha = self.grating_angle.to(u.rad)
        self.beta = self.camera_angle.to(u.rad) - self.grating_angle.to(u.rad)

        self.center_wavelength = (np.sin(self.alpha) +
                                  np.sin(self.beta)) / self.grating_frequency

        limit_angle = np.arctan(
            self.pixel_count *
            (self.pixel_size / self.goodman_focal_length) / 2)

        self.blue_limit = ((np.sin(self.alpha) +
                            np.sin(self.beta - limit_angle.to(u.rad))) /
                           self.grating_frequency).to(u.angstrom) + \
                          blue_correction_factor

        self.red_limit = ((np.sin(self.alpha) +
                           np.sin(self.beta +
                                  limit_angle.to(u.rad))) /
                          self.grating_frequency).to(u.angstrom) +\
                         red_correction_factor

        pixel_one = 0
        pixel_two = 0
        log.debug('Center Wavelength : {:.3f} Blue Limit : '
                  '{:.3f} Red Limit : {:.3f} '.format(self.center_wavelength,
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

    def interpolate(self, spectrum):
        """Creates an interpolated version of the input spectrum

        This method creates an interpolated version of the input array, it is
        used mainly for a spectrum but it can also be used with any
        unidimensional array, assuming you are happy with the interpolation_size
        attribute defined for this class. The reason for doing interpolation is
        that it allows to find the lines and its respective center more
        precisely. The default interpolation size is 200 (two hundred) points.

        Args:
            spectrum (array): an uncalibrated spectrum or any unidimensional
                array.

        Returns:
            Two dimensional array containing x-axis and interpolated array.
                The x-axis preserves original pixel values.

        """
        x_axis = range(spectrum.size)
        first_x = x_axis[0]
        last_x = x_axis[-1]

        new_x_axis = np.linspace(first_x,
                                 last_x,
                                 spectrum.size * self.interpolation_size)

        tck = scipy.interpolate.splrep(x_axis, spectrum, s=0)
        new_spectrum = scipy.interpolate.splev(new_x_axis, tck, der=0)
        return [new_x_axis, new_spectrum]

    def recenter_line_by_data(self, data_name, x_data):
        """Finds a better center for a click-selected line

        This method is called by another method that handles click events. An
        argument is parsed that will tell which plot was clicked and what is
        the x-value in data coordinates. Then the closest pixel center will be
        found and from there will extract a 20 pixel wide sample of the data
        (this could be a future improvement: the width of the extraction should
        depend on the FWHM of the lines). The sample of the data is used to
        calculate a centroid (center of mass) which is a good approximation but
        could be influenced by data shape or if the click was too far
        (unquantified yet). That is why for the reference data, a database of
        laboratory line center will be
        used and for raw data the line centers are calculated earlier in the
        process, independently of any human input.

        It wil also plot the sample, the centroid and the reference line center
        at the fourth subplot, bottom right corner.

        Args:
            data_name (str): 'reference' or 'raw-data' is where the click was
                done
            x_data (float): click x-axis value in data coordinates

        Returns:
            reference_line_value (float): The value of the line center that will
            be used later to do the wavelength fit

        """
        if data_name == 'reference':
            pseudo_center = np.argmin(abs(self.reference_solution[0] - x_data))

            reference_line_index = np.argmin(
                abs(self.reference_data.get_line_list_by_name(
                    self.lamp_name) - x_data))

            reference_line_value = self.reference_data.get_line_list_by_name(
                self.lamp_name)[reference_line_index]

            sub_x = self.reference_solution[0][
                    pseudo_center - 10: pseudo_center + 10]

            sub_y = self.reference_solution[1][
                    pseudo_center - 10: pseudo_center + 10]

            center_of_mass = np.sum(sub_x * sub_y) / np.sum(sub_y)
            # print 'centroid ', center_of_mass
            # plt.figure(3)
            # if self.ax4_plots is not None or self.ax4_com is not None or
            # self.ax4_rlv is not None:
            try:
                self.ax4.cla()
                self.ax4.relim()
            except NameError as err:
                log.error(err)

            self.ax4.set_title('Reference Data Clicked Line')
            self.ax4.set_xlabel('Wavelength (Angstrom)')
            self.ax4.set_ylabel('Intensity (Counts)')

            self.ax4_plots = self.ax4.plot(sub_x,
                                           sub_y,
                                           color='k',
                                           label='Data')

            self.ax4_rlv = self.ax4.axvline(reference_line_value,
                                            linestyle='-',
                                            color='r',
                                            label='Reference Line Value')

            self.ax4_com = self.ax4.axvline(center_of_mass,
                                            linestyle='--',
                                            color='b',
                                            label='Centroid')

            self.ax4.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
            self.ax4.legend(loc=3, framealpha=0.5)
            self.i_fig.canvas.draw()
            # return center_of_mass
            return reference_line_value
        elif data_name == 'raw-data':
            pseudo_center = np.argmin(abs(self.raw_pixel_axis - x_data))
            raw_line_index = np.argmin(abs(self.lines_center - x_data))
            raw_line_value = self.lines_center[raw_line_index]
            # print(raw_line_value)
            sub_x = self.raw_pixel_axis[pseudo_center - 10: pseudo_center + 10]
            sub_y = self.lamp_data[pseudo_center - 10: pseudo_center + 10]
            center_of_mass = np.sum(sub_x * sub_y) / np.sum(sub_y)
            # print 'centroid ', center_of_mass
            # plt.figure(3)
            # if self.ax4_plots is not None or self.ax4_com is not None or
            # self.ax4_rlv is not None:
            try:
                self.ax4.cla()
                self.ax4.relim()
            except NameError as err:
                log.error(err)
            self.ax4.set_title('Raw Data Clicked Line')
            self.ax4.set_xlabel('Pixel Axis')
            self.ax4.set_ylabel('Intensity (Counts)')

            self.ax4_plots = self.ax4.plot(sub_x,
                                           sub_y,
                                           color='k',
                                           label='Data')

            self.ax4_rlv = self.ax4.axvline(raw_line_value,
                                            linestyle='-',
                                            color='r',
                                            label='Line Center')

            self.ax4_com = self.ax4.axvline(center_of_mass,
                                            linestyle='--',
                                            color='b',
                                            label='Centroid')

            self.ax4.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
            self.ax4.legend(loc=3, framealpha=0.5)
            self.i_fig.canvas.draw()
            # return center_of_mass
            return raw_line_value
        else:
            log.error('Unrecognized data name')

    def predicted_wavelength(self, pixel):
        """Find the predicted wavelength value for a given pixel

        It is possible to estimate the wavelength position of any pixel given
        the instrument configuration.

        Notes:
            The equations are not precise enough so the value returned here has
            to be used as an estimate only.

        Args:
            pixel (int): Pixel number.

        Returns:
            Wavelength value in angstrom.

        """
        # TODO (simon): Update with bruno's new calculations
        alpha = self.alpha
        beta = self.beta
        # pixel_count = self.pixel_count
        binning = self.serial_binning

        grating_frequency = self.grating_frequency

        wavelength = 10 * (1e6 / grating_frequency) * \
                     (np.sin(alpha * np.pi / 180.) +
                      np.sin((beta * np.pi / 180.) +
                             np.arctan((pixel * binning - 2048) *
                                       0.015 / 377.2)))
        return wavelength

    def automatic_wavelength_solution(self):
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
            # reference_lamp_file = self.reference_data.get_best_reference_lamp(
            #     header=self.lamp_header)
            reference_lamp_file = self.reference_data.get_exact_lamp(
                header=self.lamp_header)

            log.debug('Found reference lamp: {:s}'.format(reference_lamp_file))

            reference_lamp_data = CCDData.read(reference_lamp_file, unit=u.adu)
        except NotImplementedError:

            log.warning('This configuration is not supported in '
                        'automatic mode.')

            # TODO (simon): Evaluate if send this to interactive mode
            return None

        read_wsolution = ReadWavelengthSolution(
            header=reference_lamp_data.header, data=reference_lamp_data.data)

        reference_lamp_wav_axis, reference_lamp_data.data = read_wsolution()
        reference_lamp_copy = reference_lamp_data.copy()

        if False:
            box_kernel = Box1DKernel(width=33)
            test = convolve(reference_lamp_copy.data, box_kernel)

        else:
            gaussian_kernel = Gaussian1DKernel(stddev=3.)

            reference_lamp_copy.data = convolve(reference_lamp_copy.data,
                                                gaussian_kernel)

        # Initialize wavelength builder class
        wavelength_solution = WavelengthFitter(model='chebyshev',
                                               degree=self.poly_order)
        # self.wsolution = wavelength_solution.ws_fit(pixel, auto_angs)

        '''detect lines in comparison lamp (not reference)'''
        lamp_lines_pixel = self.get_lines_in_lamp()
        lamp_lines_angst = read_wsolution.math_model(lamp_lines_pixel)

        pixel_values = []
        angstrom_values = []
        correlation_values = []
        angstrom_differences = []

        log.debug('Length {:d}'.format(len(self.lamp_data)))
        log.debug('NLines {:d}'.format(len(lamp_lines_pixel)))

        log.debug('Length / NLines {:.3f}'.format(
            len(self.lamp_data) / float(len(lamp_lines_pixel))))

        half_width = int((len(self.lamp_data) /
                          float(len(lamp_lines_pixel))) / 2.)

        for i in range(len(lamp_lines_pixel)):
            line_value_pixel = lamp_lines_pixel[i]
            line_value_angst = lamp_lines_angst[i]
            # half_width = 100

            xmin = int(max(0, round(line_value_pixel - half_width)))

            xmax = int(min(round(line_value_pixel + half_width),
                           len(self.lamp_data)))

            if xmin >= xmax:
                continue
            # print(xmin, xmax, self.lamp_data.size)
            # TODO (simon): Convolve to match wider lines such as those from
            # TODO (cont): the slit of 5 arseconds
            ref_sample = reference_lamp_data.data[xmin:xmax]
            ref_wavele = reference_lamp_wav_axis[xmin:xmax]
            lamp_sample = self.lamp_data[xmin:xmax]

            correlation_value = self.cross_correlation(ref_sample, lamp_sample)
            log.debug('Cross correlation value '
                      '{:s}'.format(str(correlation_value)))

            """record value for reference wavelength"""
            angstrom_value_model = read_wsolution.math_model(
                line_value_pixel + correlation_value)

            correlation_values.append(correlation_value)
            angstrom_differences.append(angstrom_value_model - line_value_angst)
            angstrom_values.append(angstrom_value_model)
            pixel_values.append(line_value_pixel)

            if False:
                plt.ion()
                plt.title('Samples after cross correlation')
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
                                    sigma=2,
                                    iters=1,
                                    cenfunc=np.ma.median)

        if np.ma.is_masked(clipped_values):
            _pixel_values = list(pixel_values)
            _angstrom_values = list(angstrom_values)
            pixel_values = []
            angstrom_values = []
            for i in range(len(clipped_values)):
                if clipped_values[i] is not np.ma.masked:
                    pixel_values.append(_pixel_values[i])
                    angstrom_values.append(_angstrom_values[i])

        # Create a wavelength solution
        log.info('Creating Wavelength Solution')

        self.wsolution = wavelength_solution.ws_fit(pixel_values,
                                                    angstrom_values)

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
                                         iters=1,
                                         cenfunc=np.ma.median)

        if np.ma.is_masked(clipped_differences):

            log.debug('Cleaning pixel to angstrom match to improve wavelength '
                      'solution')

            _pixel_values = list(pixel_values)
            _angstrom_values = list(angstrom_values)
            pixel_values = []
            angstrom_values = []
            for i in range(len(clipped_differences)):
                if clipped_differences[i] is not np.ma.masked:
                    pixel_values.append(_pixel_values[i])
                    angstrom_values.append(_angstrom_values[i])
            log.info('Re-fitting wavelength solution')

            self.wsolution = wavelength_solution.ws_fit(pixel_values,
                                                        angstrom_values)

        self.evaluate_solution()

        if self.args.plot_results or self.args.debug_mode:
            plt.switch_backend('Qt4Agg')
            if self.i_fig is None:
                self.i_fig = plt.figure(figsize=(15, 10))
                self.i_fig.canvas.set_window_title(
                    'Automatic Wavelength Solution')
            else:
                self.i_fig.clf()

            if not self.args.debug_mode:
                plt.ion()
                # plt.show()
            else:
                plt.ioff()

            manager = plt.get_current_fig_manager()

            if plt.get_backend() == u'GTK3Agg':
                manager.window.maximize()
            elif plt.get_backend() == u'Qt4Agg':
                manager.window.showMaximized()

            self.ax1 = self.i_fig.add_subplot(111)

            self.ax1.plot([], color='m', label='Pixels')
            self.ax1.plot([], color='c', label='Angstrom')
            for val in pixel_values:
                self.ax1.axvline(self.wsolution(val), color='m')
            for val2 in angstrom_values:
                self.ax1.axvline(val2, color='c', linestyle='--')

            self.ax1.plot(reference_lamp_wav_axis,
                          reference_lamp_data.data,
                          label='Reference',
                          color='k',
                          alpha=1)

            self.ax1.plot(self.wsolution(self.raw_pixel_axis),
                          self.lamp_data,
                          label='Last Solution',
                          color='r',
                          alpha=0.7)

            try:
                wavmode = self.lamp_header['wavmode']
            except KeyError as error:
                log.debug(error)
                wavmode = ''

            self.ax1.set_xlabel('Wavelength (Angstrom)')
            self.ax1.set_ylabel('Intensity (ADU)')

            self.ax1.set_title('Automatic Wavelength Solution\n'
                               + self.lamp_header['OBJECT']
                               + ' ' + wavmode + '\n'
                               + 'RMS Error: {:.3f}'.format(self.rms_error))

            self.ax1.legend(loc='best')
            self.i_fig.tight_layout()

            if self.args.save_plots:
                # saves pdf files of the wavelength solution plot
                out_file_name = 'automatic-solution_' + self.lamp_header[
                    'OBJECT']

                file_count = len(glob.glob(os.path.join(self.args.destiny,
                                                        out_file_name + '*')))

                out_file_name += '_{:04d}.pdf'.format(file_count)
                pdf_pages = PdfPages(
                    os.path.join(self.args.destiny, out_file_name))
                plt.savefig(pdf_pages, format='pdf')
                pdf_pages.close()

                # saves png images
                plots_path = os.path.join(self.args.destiny, 'plots')
                if not os.path.isdir(plots_path):
                    os.path.os.makedirs(plots_path)
                plot_name = os.path.join(plots_path, out_file_name + '.png')
                plt.savefig(plot_name, dpi=300)
            if self.args.debug_mode:
                plt.show()
            else:
                plt.draw()
                plt.pause(1)
                plt.ioff()
                plt.close()
                # plt.close(self.i_fig)

    def interactive_wavelength_solution(self, object_name=''):
        """Find the wavelength solution interactively

        This method uses the graphical capabilities of matplotlib in particular
        the user interface (UI) capabilities such as click and key-pressed
        events. We could say that it implements a matplotlib Graphical User
        Interface. There are some limitation though, for instance, the access to
        a better formatted help text.

        It will display four plots of the same height but their width have a
        4:1 ratio being the leftmost the wider.

        Top Left Plot: Displays the spectrum of the lamp we want to calibrate,
        or as it is called in the plot, the Raw Data, therefore the axis are
        intensity (ADU) versus dispersion (pixels). There will be red vertical
        dashed lines that are lines detected by pipeline. If you see the a line
        without a red dashed line it means it was not detected.

        Bottom Left Plot: Displays the reference lamp data, therefore the axis
        are intensity (ADU) versus wavelength (Angstrom). In this case there are
        vertical dotted lines plotted with a relatively low alpha value. They
        represent reference laboratory lines obtained from the NIST database at
        https://physics.nist.gov/PhysRefData/ASD/lines_form.html.

        Top Right Plot: This is a quick reference help text. The formatting of
        the text is quite difficult here due to matplotlib's own limitations.
        This text is also displayed on the terminal.

        Bottom Right Plot: This is a dynamic plot since it shows different kind
        of information. First it shows short messages, such as warnings, errors,
        or general information. Once you select a line it will show a zoomed
        version of the line plus the reference value chosen as well as were the
        click was done. And finally, it shows the dispersion of the wavelength
        solution. This could be one of the most useful plots in the entire
        screen.

        Usage: A quick reference of the instructions is shown in the Top Right
        Plot, but in general you have to select the matching lines in both data
        plots this will create a table of pixels versus angstrom which at its
        roots a wavelength solution. Then this will have to be translated to a
        mathematical model. Having this mathematical model the wavelength
        solution can be applied, evaluated etc. Among the features implemented
        you can:
          - Delete a line or a matching pair.
          - Find more lines automatically once you have a decent wavelength
            solution.
          - Remove all the automatically added lines at once, this is useful in
            some cases when adding more lines actually makes the solution worst.
          - Find the RMS Error of the solution
          - Remove ALL the points at a given point to start over.
          - Fit a wavelength solution



        Notes:
            This method uses the Qt4Agg backend, in theory it could also work
            with GTK3Agg but it is being forced to use Qt4Agg.

        """
        plt.switch_backend('Qt4Agg')

        # disable full screen to allow the use of f for fitting the solution

        plt.rcParams['keymap.fullscreen'] = [u'ctrl+f']

        try:
            reference_file = self.reference_data.get_best_reference_lamp(
                header=self.lamp_header)
        except NotImplementedError:
            reference_file = self.reference_data.get_reference_lamps_by_name(
                lamp_name=self.lamp_header['OBJECT'])
            log.warning('Could not find a perfect match for reference data')
            # reference_file = None
            # log.critical('Could not find a comparison lamp in the reference.')

        # reference_file = self.reference_data.get_reference_lamps_by_name(
        #     self.lamp_name)

        if reference_file is not None:
            log.info('Using reference file: {:s}'.format(reference_file))
            reference_plots_enabled = True
            ref_data = fits.getdata(reference_file)
            ref_header = fits.getheader(reference_file)
            fits_ws_reader = ReadWavelengthSolution(ref_header, ref_data)
            self.reference_solution = fits_ws_reader()
        else:
            reference_plots_enabled = False
            log.error('Please Check the OBJECT Keyword of your reference data')

        # ------- Plots -------
        self.i_fig, ((self.ax1, self.ax2), (self.ax3, self.ax4)) = \
            plt.subplots(2,
                         2,
                         gridspec_kw={'width_ratios': [4, 1]})

        self.i_fig.canvas.set_window_title('Science Target: {:s}'.format(
            object_name))

        manager = plt.get_current_fig_manager()
        if plt.get_backend() == u'GTK3Agg':
            manager.window.maximize()
        elif plt.get_backend() == u'Qt4Agg':
            manager.window.showMaximized()

        self.ax1.set_title('Raw Data - {:s}\n{:s} - {:s}'.format(
            self.lamp_name,
            self.lamp_header['GRATING'],
            self.lamp_header['SLIT']))

        self.ax1.set_xlabel('Pixels')
        self.ax1.set_ylabel('Intensity (counts)')
        self.ax1.plot([], linestyle='--', color='r', label='Detected Lines')
        for idline in self.lines_center:
            self.ax1.axvline(idline, linestyle='--', color='r')

        self.ax1.plot(self.raw_pixel_axis,
                      self.lamp_data,
                      color='k',
                      label='Raw Data')

        self.ax1.set_xlim((0, len(self.lamp_data)))
        self.ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
        self.ax1.legend(loc=2)

        # Update y limits to have an extra 5% to top and bottom
        ax1_ylim = self.ax1.get_ylim()
        ax1_y_range = ax1_ylim[1] - ax1_ylim[0]

        self.ax1.set_ylim((ax1_ylim[0] - 0.05 * ax1_y_range,
                           ax1_ylim[1] + 0.05 * ax1_y_range))

        self.ax3.set_title('Reference Data - {:s}'.format(self.lamp_name))
        self.ax3.set_xlabel('Wavelength (Angstrom)')
        self.ax3.set_ylabel('Intensity (counts)')
        self.ax3.set_xlim((self.blue_limit.value, self.red_limit.value))
        self.ax3.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

        self.ax3.plot([],
                      linestyle=':',
                      color='r',
                      label='Reference Line Values')

        for rline in self.reference_data.get_line_list_by_name(self.lamp_name):
            self.ax3.axvline(rline, linestyle=':', color='r')

        if reference_plots_enabled:

            self.ax3.plot(self.reference_solution[0],
                          self.reference_solution[1],
                          color='k',
                          label='Reference Lamp Data')

        self.ax3.legend(loc=2)

        # Update y limits to have an extra 5% to top and bottom
        ax3_ylim = self.ax3.get_ylim()
        ax3_y_range = ax3_ylim[1] - ax3_ylim[0]

        self.ax3.set_ylim((ax3_ylim[0] - 0.05 * ax3_y_range,
                           ax3_ylim[1] + 0.05 * ax3_y_range))
        # print(ax3_ylim)

        self.display_help_text()

        # zoomed plot
        self.display_onscreen_message('Use middle button click to select a '
                                      'line')

        plt.subplots_adjust(left=0.05,
                            right=0.99,
                            top=0.96,
                            bottom=0.04,
                            hspace=0.17,
                            wspace=0.11)

        self.raw_data_bb = self.ax1.get_position()
        self.reference_bb = self.ax3.get_position()
        self.contextual_bb = self.ax4.get_position()

        # if self.click_input_enabled:
        self.i_fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.i_fig.canvas.mpl_connect('key_press_event', self.key_pressed)
        # print self.wsolution
        plt.show()
        return True

    def cross_correlation(self, reference, new_array, mode='full'):
        """Do cross correlation to two arrays

        Args:
            reference (array): Reference array.
            new_array (array): Array to be matched.
            mode (str): Correlation mode for `scipy.signal.correlate`.

        Returns:
            correlation_value (int): Shift value in pixels.

        """
        # print(reference, new_array)
        cyaxis2 = new_array
        if float(re.sub('[A-Za-z" ]', '', self.lamp_header['SLIT'])) > 3:

            box_width = float(
                re.sub('[A-Za-z" ]', '', self.lamp_header['SLIT'])) / 0.15

            log.debug('BOX WIDTH: {:f}'.format(box_width))
            box_kernel = Box1DKernel(width=box_width)
            max_before = np.max(reference)
            cyaxis1 = convolve(reference, box_kernel)
            max_after = np.max(cyaxis1)
            cyaxis1 *= max_before / max_after


        else:
            gaussian_kernel = Gaussian1DKernel(stddev=2.)
            cyaxis1 = convolve(reference, gaussian_kernel)
            cyaxis2 = convolve(new_array, gaussian_kernel)
        # plt.plot(cyaxis1, color='k', label='Reference')
        # plt.plot(cyaxis2, color='r', label='New Array')
        # plt.plot(reference, color='g')
        # plt.show()
        try:
            ccorr = signal.correlate(cyaxis1, cyaxis2, mode=mode)
        except ValueError:
            print(cyaxis1, cyaxis2)
        # print('Corr ', ccorr)
        max_index = np.argmax(ccorr)

        x_ccorr = np.linspace(-int(len(ccorr) / 2.),
                              int(len(ccorr) / 2.),
                              len(ccorr))

        correlation_value = x_ccorr[max_index]
        if False:
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

    def on_click(self, event):
        """Handles Click events for Interactive Mode

        Calls the method register_mark

        Args:
            event (object): Click event
        """
        if event.button == 2:
            self.register_mark(event)
        # TODO (simon): Make sure the text below is useless
        # else:
        #     print(event.button)
        # elif event.button == 3:
        #     if len(self.reference_marks_x) == len(self.raw_data_marks_x):
        #         self.click_input_enabled = False
        #         log.info('Leaving interactive mode')
        #     else:
        #         if len(self.reference_marks_x) < len(self.raw_data_marks_x):
        #             log.info('There is {:d} click missing in the '
        #                      'Reference plot'.format(
        #                 len(self.raw_data_marks_x) -
        #                 len(self.reference_marks_x)))
        #         else:
        #             log.info('There is {:d} click missing in the '
        #                      'New Data plot'.format(
        #                 len(self.reference_marks_x) -
        #                 len(self.raw_data_marks_x)))

    def key_pressed(self, event):
        """Key event handler

        There are several key events that need to be taken care of.
        See a brief description of each one of them below:

        F1 or ?: Prints a help message
        F2 or f: Fit wavelength solution model.
        F3 or a: Find new lines.
        F4: Evaluate solution
        F6 or l: Linearize data although this is done automatically after the
        wavelength function is fit
        d: deletes closest point
        ctrl+d: deletes all recorded marks
        ctrl+z: Reverts the action of F3 or a.
        Middle Button Click or m: records data location.
        Enter: Close figure and apply solution if exists.
        Shift+Enter: Close the program with sys.exit(0)

        Notes:
            This method must be simplified

        Args:
            event (object): Key pressed event

        """
        self.events = True
        if event.key == 'f1' or event.key == '?':
            log.info('Print help regarding interactive mode')
            print("F1 or ?: Prints Help.")
            print("F2 or f: Fit wavelength solution model.")
            print("F3 or a: Find new lines.")
            print("F4: Evaluate solution")
            print("F6 or l: Linearize data (for testing not definitive)")
            print("d: deletes closest point")
            # print("l : resample spectrum to a linear dispersion axis")
            print("ctrl+d: deletes all recorded marks")
            print("ctrl+z: Go back to previous solution "
                  "(deletes automatic added points")
            print('Middle Button Click: records data location.')
            print("Enter: Close figure and apply solution if exists.")
        elif event.key == 'f2' or event.key == 'f':
            log.debug('Calling function to fit wavelength Solution')
            self.fit_pixel_to_wavelength()
            self.plot_raw_over_reference()
        elif event.key == 'f3' or event.key == 'a':
            if self.wsolution is not None:
                self.find_more_lines()
                self.update_marks_plot('reference')
                self.update_marks_plot('raw_data')
        elif event.key == 'f4':
            if self.wsolution is not None and len(self.raw_data_marks_x) > 0:
                self.evaluate_solution(plots=True)
        elif event.key == 'f5' or event.key == 'd':
            # TODO (simon): simplify this code.

            figure_x, figure_y = \
                self.i_fig.transFigure.inverted().transform((event.x, event.y))

            if self.raw_data_bb.contains(figure_x, figure_y):
                log.debug('Deleting raw point')
                # print abs(self.raw_data_marks_x - event.xdata) a[:] =

                closer_index = int(np.argmin(
                    [abs(list_val - event.xdata) for list_val in
                     self.raw_data_marks_x]))

                # print 'Index ', closer_index
                if len(self.raw_data_marks_x) == len(self.reference_marks_x):
                    self.raw_data_marks_x.pop(closer_index)
                    self.raw_data_marks_y.pop(closer_index)
                    self.reference_marks_x.pop(closer_index)
                    self.reference_marks_y.pop(closer_index)
                    self.update_marks_plot('reference')
                    self.update_marks_plot('raw_data')
                else:
                    if closer_index == len(self.raw_data_marks_x) - 1:
                        self.raw_data_marks_x.pop(closer_index)
                        self.raw_data_marks_y.pop(closer_index)
                        self.update_marks_plot('raw_data')
                    else:
                        self.raw_data_marks_x.pop(closer_index)
                        self.raw_data_marks_y.pop(closer_index)
                        self.reference_marks_x.pop(closer_index)
                        self.reference_marks_y.pop(closer_index)
                        self.update_marks_plot('reference')
                        self.update_marks_plot('raw_data')

            elif self.reference_bb.contains(figure_x, figure_y):
                log.debug('Deleting reference point')
                # print 'reference ', self.reference_marks_x, self.re
                # print self.reference_marks_x
                # print abs(self.reference_marks_x - event.xdata)

                closer_index = int(np.argmin(
                    [abs(list_val - event.xdata) for list_val in
                     self.reference_marks_x]))

                if len(self.raw_data_marks_x) == len(self.reference_marks_x):
                    self.raw_data_marks_x.pop(closer_index)
                    self.raw_data_marks_y.pop(closer_index)
                    self.reference_marks_x.pop(closer_index)
                    self.reference_marks_y.pop(closer_index)
                    self.update_marks_plot('reference')
                    self.update_marks_plot('raw_data')
                else:
                    if closer_index == len(self.reference_marks_x) - 1:
                        self.reference_marks_x.pop(closer_index)
                        self.reference_marks_y.pop(closer_index)
                        self.update_marks_plot('reference')
                    else:
                        self.raw_data_marks_x.pop(closer_index)
                        self.raw_data_marks_y.pop(closer_index)
                        self.reference_marks_x.pop(closer_index)
                        self.reference_marks_y.pop(closer_index)
                        self.update_marks_plot('reference')
                        self.update_marks_plot('raw_data')
                    
            elif self.contextual_bb.contains(figure_x, figure_y):
                log.warning("Can't delete points from here because points "
                            "represent each detected line in the raw data.")
                # closer_index_ref = int(np.argmin(
                #     [abs(list_val - event.ydata) for list_val in
                #      self.reference_marks_x]))
                #
                # closer_index_raw = int(np.argmin(
                #     [abs(list_val - event.xdata) for list_val in
                #      self.raw_data_marks_x]))
                #
                # log.debug('Raw Index {:d}, Ref Index {:d}'.format(
                #     closer_index_raw,
                #     closer_index_ref))
                # self.raw_data_marks_x.pop(closer_index_raw)
                # self.raw_data_marks_y.pop(closer_index_raw)
                # self.reference_marks_x.pop(closer_index_raw)
                # self.reference_marks_y.pop(closer_index_raw)
                # self.update_marks_plot('reference')
                # self.update_marks_plot('raw_data')
                # self.update_marks_plot('eval_plots')
        # elif event.key == 'i':
        #     figure_x, figure_y = self.i_fig.transFigure.inverted().transform(
        #         (event.x, event.y))
        #     if self.contextual_bb.contains(figure_x, figure_y):
        #         log.debug("Trying to identify point.")
        #
        #         closer_x_index = int(np.argmin(
        #             [abs(list_val - event.xdata) for list_val in
        #              self.raw_data_marks_x]))
        #         print(self.raw_data_marks_x[closer_x_index])
        #         print(self.reference_marks_x[closer_x_index])


        elif event.key == 'f6' or event.key == 'l':
            log.info('Linearize and smoothing spectrum')
            if self.wsolution is not None:
                self.linearize_spectrum(self.lamp_data, plots=True)

        elif event.key == 'ctrl+z':
            log.info('Deleting automatic added points. If exist.')

            if self.raw_data_marks_x is not [] and \
                            self.reference_marks_x is not []:

                to_remove = []
                for i in range(len(self.raw_data_marks_x)):
                    # print self.raw_data_marks[i], self.filling_value
                    if self.raw_data_marks_y[i] == self.filling_value:
                        to_remove.append(i)
                        # print to_remove
                to_remove = np.array(sorted(to_remove, reverse=True))
                if len(to_remove) > 0:
                    for index in to_remove:
                        self.raw_data_marks_x.pop(index)
                        self.raw_data_marks_y.pop(index)
                        self.reference_marks_x.pop(index)
                        self.reference_marks_y.pop(index)
                    self.update_marks_plot('reference')
                    self.update_marks_plot('raw_data')
                    # else:
                    # print self.raw_click_plot, self.ref_click_plot, 'mmm'

        elif event.key == 'ctrl+d':
            try:
                log.info('Deleting all recording Clicks')
                self.display_onscreen_message(
                    message='All points deleted')
                self.reference_marks_x = []
                self.reference_marks_y = []
                self.raw_data_marks_x = []
                self.raw_data_marks_y = []
                self.update_marks_plot('delete')
                self.plot_raw_over_reference(remove=True)
                log.info('All points deleted!')
            except:
                log.error('No points deleted')

        elif event.key == 'enter':
            if self.wsolution is not None:
                log.info('Closing figure')
                plt.close('all')
            else:
                message = 'There is still no wavelength solution!'
                log.info(message)
                self.display_onscreen_message(message)

        elif event.key == 'm':
            self.register_mark(event)

        elif event.key == 'ctrl+q':
            log.info('Pressed Ctrl+q. Closing the program')
            sys.exit(0)

        else:
            log.debug("No action for key pressed: {:s}".format(event.key))
            pass

    def register_mark(self, event):
        """Marks a line

        Detects where the click was done or m-key was pressed and calls the
        corresponding method. It handles the middle button click and m-key being
        pressed. There are two regions of interest as for where a click was
        done. The raw and reference data respectively. For any of such regions
        it will call the method that recenter the line and once the desired
        value is returned it will be appended to the list that contains all the
        correspondent line positions, raw (pixels) and reference (angstrom)

        Args:
            event (object): Click or m-key pressed event
        """
        if event.xdata is not None and event.ydata is not None:
            figure_x, figure_y = \
                self.i_fig.transFigure.inverted().transform((event.x, event.y))

            if self.reference_bb.contains(figure_x, figure_y):
                # self.reference_marks.append([event.xdata, event.ydata])
                self.reference_marks_x.append(
                    self.recenter_line_by_data('reference', event.xdata))

                self.reference_marks_y.append(event.ydata)
                self.update_marks_plot('reference')
            elif self.raw_data_bb.contains(figure_x, figure_y):
                # self.raw_data_marks.append([event.xdata, event.ydata])
                self.raw_data_marks_x.append(
                    self.recenter_line_by_data('raw-data', event.xdata))

                self.raw_data_marks_y.append(event.ydata)
                self.update_marks_plot('raw_data')
            else:
                log.debug('{:f} {:f} Are not contained'.format(figure_x,
                                                               figure_y))
        else:
            log.error('Clicked Region is out of boundaries')

    def find_more_lines(self):
        """Method to add more lines given that a wavelength solution already
        exists

        This method is part of the interactive wavelength solution mechanism.
        If a wavelength solution exist it uses the line centers in pixels to
        estimate their respective wavelength and then search for the closest
        value in the list of reference lines for the elements in the comparison
        lamp. Then it filters the worst of them by doing sigma clipping.
        Finally it adds them to the class' attributes that contains the list of
        reference points.

        Better results are obtained if the solution is already decent. Visual
        inspection also improves final result.
        """
        new_physical = []
        new_wavelength = []
        square_differences = []
        if self.wsolution is not None:
            wlines = self.wsolution(self.lines_center)
            for i in range(len(wlines)):
                # [abs(list_val - wlines[i]) for list_val in \
                # self.reference_data.get_line_list_by_name(self.lamp_name)]

                closer_index = np.argmin(
                    abs(self.reference_data.get_line_list_by_name(
                        self.lamp_name) - wlines[i]))

                rline = self.reference_data.get_line_list_by_name(
                    self.lamp_name)[closer_index]

                rw_difference = wlines[i] - rline
                # print('Difference w - r ', rw_difference, rline)
                square_differences.append(rw_difference ** 2)
                new_physical.append(self.lines_center[i])
                new_wavelength.append(rline)
            clipped_differences = sigma_clip(square_differences,
                                             sigma=2,
                                             iters=3)

            if len(new_wavelength) == len(new_physical) == \
                    len(clipped_differences):

                for i in range(len(new_wavelength)):
                    if clipped_differences[i] is not \
                            np.ma.masked and new_wavelength[i] not in \
                            self.reference_marks_x:

                        self.reference_marks_x.append(new_wavelength[i])
                        self.reference_marks_y.append(self.filling_value)
                        self.raw_data_marks_x.append(new_physical[i])
                        self.raw_data_marks_y.append(self.filling_value)
        return True

    def update_marks_plot(self, action=None):
        """Update the points that represent marks on lamp plots

        When you mark a line a red dot marks the position of the line at the
        exact y location of the click, for the x location it will use the value
        obtained by means of the recentering method. There are three possible
        actions: Update the reference plot's click, the raw data marks or
        delete them all.

        Args:
            action (str): A string that could be 'reference', 'raw_data' or
            'delete' depending on the action desired
        """
        # print(type(action), type(pixel_axis), type(differences))
        if action == 'reference':
            if self.points_ref is not None:
                try:
                    self.points_ref.remove()
                    self.ax3.relim()
                except:
                    pass
            self.points_ref, = self.ax3.plot(self.reference_marks_x,
                                             self.reference_marks_y,
                                             linestyle='None',
                                             marker='o',
                                             color='r')
            self.i_fig.canvas.draw()
        elif action == 'raw_data':
            # print self.points_raw
            # print dir(self.points_raw)
            if self.points_raw is not None:
                try:
                    self.points_raw.remove()
                    self.ax1.relim()
                except ValueError as err:
                    log.error(err)
            self.points_raw, = self.ax1.plot(self.raw_data_marks_x,
                                             self.raw_data_marks_y,
                                             linestyle='None',
                                             marker='o',
                                             color='r')
            self.i_fig.canvas.draw()
        elif action == 'delete':
            if self.points_raw is not None and self.points_ref is not None:
                self.points_raw.remove()
                self.ax1.relim()
                self.points_ref.remove()
                self.ax3.relim()
                self.i_fig.canvas.draw()
        else:
            log.error('Unknown Action {:s}'.format(action))

    def plot_raw_over_reference(self, remove=False):
        """Overplot raw data over reference lamp using current wavelength
        solution model

        Once the wavelength solution is obtained this method is called to apply
        the already mentioned solution to the raw data and then overplot it on
        the reference lamp plot. This is very useful to have a visual idea of
        how far or close the solution is.

        Args:
            remove (bool): True or False depending whether you want to remove
            the overplotted lamp or not
        """
        if self.wsolution is not None:
            if self.line_raw is not None:
                try:
                    self.line_raw.remove()
                    self.ax3.relim()
                except:
                    pass
            if not remove:
                # TODO(simon): catch TypeError Exception and correct what is
                # TODO (cont): causing it
                self.line_raw, = self.ax3.plot(
                    self.wsolution(self.raw_pixel_axis),
                    self.lamp_data,
                    linestyle='-',
                    color='r',
                    alpha=0.4,
                    label='New Wavelength Solution')

            self.ax3.legend(loc=2)
            self.i_fig.canvas.draw()

    def evaluate_solution(self, plots=False):
        """Calculate the Root Mean Square Error of the solution

        Once the wavelength solution is obtained it has to be evaluated. The
        line centers found for the raw comparison lamp will be converted to,
        according to the new solution, angstrom. Then for each line the closest
        reference line value is obtained. The difference is stored. Then this
        differences are cleaned by means of a sigma clipping method that will
        rule out any outlier or any line that is not well matched. Then, using
        the sigma clipped differences the Root Mean Square error is calculated.

        It also creates a plot in the bottom right subplot of the interactive
        window, showing an scatter plot plus some information regarding the
        quality of the fit.

        Args:
            plots (bool): Whether to create the plot or not

        Returns:
            results (list): Contains three elements: rms_error (float),
            npoints (int), n_rejections (int)

        """
        if self.wsolution is not None:
            differences = np.array([])
            wavelength_line_centers = self.wsolution(self.lines_center)

            for wline in wavelength_line_centers:
                closer_index = np.argmin(
                    abs(self.reference_data.get_line_list_by_name(
                        self.lamp_name) - wline))

                rline = self.reference_data.get_line_list_by_name(
                    self.lamp_name)[closer_index]

                rw_difference = wline - rline
                # print 'Difference w - r ', rw_difference, rline
                differences = np.append(differences, rw_difference)

            clipping_sigma = 2.
            # print(differences)
            clipped_differences = sigma_clip(differences,
                                             sigma=clipping_sigma,
                                             iters=5,
                                             cenfunc=np.ma.median)

            once_clipped_differences = sigma_clip(differences,
                                                  sigma=clipping_sigma,
                                                  iters=1,
                                                  cenfunc=np.ma.median)

            npoints = len(clipped_differences)
            n_rejections = np.ma.count_masked(clipped_differences)
            square_differences = []
            for i in range(len(clipped_differences)):
                if clipped_differences[i] is not np.ma.masked:
                    square_differences.append(clipped_differences[i] ** 2)
            old_rms_error = None
            if self.rms_error is not None:
                old_rms_error = float(self.rms_error)
            self.rms_error = np.sqrt(
                np.sum(square_differences) / len(square_differences))

            log.info('RMS Error : {:.3f}'.format(self.rms_error))

            if plots:
                if self.ax4_plots is not None or \
                                self.ax4_com is not None or \
                                self.ax4_rlv is not None:

                    try:
                        self.ax4.cla()
                        self.ax4.relim()
                    except NameError as err:
                        log.error(err)

                self.ax4.set_title('RMS Error {:.3f} \n'
                                   '{:d} points ({:d} '
                                   'rejected)'.format(self.rms_error,
                                                      npoints,
                                                      n_rejections))
                
                self.ax4.set_ylim(once_clipped_differences.min(),
                                  once_clipped_differences.max())

                # self.ax4.set_ylim(- rms_error, 2 * rms_error)
                self.ax4.set_xlim(np.min(self.lines_center),
                                  np.max(self.lines_center))

                self.ax4_rlv = self.ax4.scatter(self.lines_center,
                                                differences,
                                                marker='x',
                                                color='k',
                                                label='Rejected Points')

                self.ax4_com = self.ax4.axhspan(clipped_differences.min(),
                                                clipped_differences.max(),
                                                color='k',
                                                alpha=0.4,
                                                edgecolor=None,
                                                label='{:.1f} Sigma'.format(
                                                    clipping_sigma))

                self.ax4_plots = self.ax4.scatter(self.lines_center,
                                                  clipped_differences,
                                                  label='Differences')

                if self.rms_error is not None and old_rms_error is not None:
                    # increment_color = 'white'
                    rms_error_difference = self.rms_error - old_rms_error

                    if rms_error_difference > 0.001:
                        increment_color = 'red'
                    elif rms_error_difference < -0.001:
                        increment_color = 'green'
                    else:
                        increment_color = 'white'

                    message = r'$\Delta$ RMSE {:+.3f}'.format(
                        rms_error_difference)

                    # self.display_onscreen_message(message=message,
                    #                               color=increment_color)

                    self.ax4.text(0.05, 0.95,
                                  message,
                                  verticalalignment='top',
                                  horizontalalignment='left',
                                  transform=self.ax4.transAxes,
                                  color=increment_color,
                                  fontsize=15)

                self.ax4.set_xlabel('Pixel Axis (Pixels)')
                self.ax4.set_ylabel('Difference (Angstroms)')

                self.ax4.legend(loc=3, framealpha=0.5)
                self.i_fig.canvas.draw()

            results = [self.rms_error, npoints, n_rejections]
            return results
        else:
            log.error('Solution is still non-existent!')

    def fit_pixel_to_wavelength(self):
        """Does the fit to find the wavelength solution

        Once you have four data points on each side (raw and reference or pixel
        and angstrom) it calculates the fit using a Chebyshev model of third
        degree. This was chosen because it worked better compared to the rest.
        There is a slight deviation from linearity in all Goodman data,
        therefore a linear model could not be used, also is said that a Spline
        of third degree is "too flexible" which I also experienced and since the
        deviation from linearity is not extreme it seemed that it was not
        necessary to implement.

        This method checks that the data that will be used as input to calculate
        the fit have the same dimensions and warns the user in case is not.

        Returns:
            None (None): An empty return is created to finish the execution of
            the method when a fit will not be possible

        """
        if len(self.reference_marks_x) and len(self.raw_data_marks_x) > 0:

            if len(self.reference_marks_x) < 4 or \
                            len(self.raw_data_marks_x) < 4:

                message = 'Not enough marks! Minimum 4 each side.'
                self.display_onscreen_message(message)
                return

            if len(self.reference_marks_x) != len(self.raw_data_marks_x):
                if len(self.reference_marks_x) < len(self.raw_data_marks_x):
                    n = len(self.raw_data_marks_x) - len(self.reference_marks_x)
                    if n == 1:
                        message_text = '{:d} Reference Click is ' \
                                       'missing!.'.format(n)
                    else:
                        message_text = '{:d} Reference Clicks are ' \
                                       'missing!.'.format(n)
                else:
                    n = len(self.reference_marks_x) - len(self.raw_data_marks_x)
                    if n == 1:
                        message_text = '{:d} Raw Click is missing!.'.format(n)
                    else:
                        message_text = '{:d} Raw Clicks are missing!.'.format(n)
                self.display_onscreen_message(message_text)
            else:
                pixel = []
                angstrom = []
                for i in range(len(self.reference_marks_x)):
                    pixel.append(self.raw_data_marks_x[i])
                    angstrom.append(self.reference_marks_x[i])

                wavelength_solution = WavelengthFitter(model='chebyshev',
                                                       degree=self.poly_order)

                self.wsolution = wavelength_solution.ws_fit(pixel, angstrom)
                self.evaluate_solution(plots=True)

        else:
            log.error('Clicks record is empty')
            self.display_onscreen_message(message='Clicks record is empty')
            if self.wsolution is not None:
                self.wsolution = None

    def linearize_spectrum(self, data, plots=False):
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
        # print(pixel_axis)
        if self.wsolution is not None:
            x_axis = self.wsolution(pixel_axis)
            new_x_axis = np.linspace(x_axis[0], x_axis[-1], len(data))
            tck = scipy.interpolate.splrep(x_axis, data, s=0)
            linearized_data = scipy.interpolate.splev(new_x_axis,
                                                      tck,
                                                      der=0)

            smoothed_linearized_data = signal.medfilt(linearized_data)
            # print('sl ', smoothed_linearized_data)
            if plots:
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

                plt.tight_layout()
                plt.legend(loc=3)
                plt.show()

                fig7 = plt.figure(7)
                plt.xlabel('Pixels')
                plt.ylabel('Angstroms')
                fig7.canvas.set_window_title('Wavelength Solution')
                plt.plot(x_axis, color='b', label='Non linear wavelength-axis')
                plt.plot(new_x_axis, color='r', label='Linear wavelength-axis')
                plt.tight_layout()
                plt.legend(loc=3)
                plt.show()

            linear_data = [new_x_axis, smoothed_linearized_data]
            return linear_data

    def add_wavelength_solution(self,
                                new_header,
                                spectrum,
                                original_filename,
                                evaluation_comment=None,
                                index=None):
        """Add wavelength solution to the new FITS header

        Defines FITS header keyword values that will represent the wavelength
        solution in the header so that the image can be read in any other
        astronomical tool. (e.g. IRAF)

        Notes:
            This method also saves the data to a new FITS file, This should be
            in separated methods to have more control on either process.

        Args:
            new_header (object): An Astropy header object
            spectrum (Array): A numpy array that corresponds to the processed
            data
            original_filename (str): Original Image file name
            evaluation_comment (str): A comment with information regarding the
            quality of the wavelength solution
            index (int): If in one 2D image there are more than one target the
            index represents the target number.

        Returns:
            new_header (object): An Astropy header object. Although not
            necessary since there is no further processing

        """
        if evaluation_comment is None:
            rms_error, n_points, n_rejections = self.evaluate_solution()
            self.evaluation_comment = 'Lamp Solution RMSE = {:.3f} ' \
                                      'Npoints = {:d}, ' \
                                      'NRej = {:d}'.format(rms_error,
                                                           n_points,
                                                           n_rejections)

            new_header['HISTORY'] = self.evaluation_comment
        else:
            new_header['HISTORY'] = evaluation_comment

        new_crpix = 1
        new_crval = spectrum[0][new_crpix - 1]
        new_cdelt = spectrum[0][new_crpix] - spectrum[0][new_crpix - 1]

        new_header['BANDID1'] = 'spectrum - background none, weights none, ' \
                                'clean no'
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
        new_header['DCLOG1'] = 'REFSPEC1 = {:s}'.format(self.calibration_lamp)

        # print(new_header['APNUM*'])
        if index is None:
            f_end = '.fits'
        else:
            f_end = '_{:d}.fits'.format(index)
        # idea
        #  remove .fits from original_filename
        # define a base original name
        # modify in to _1, _2 etc in case there are multitargets
        # add .fits

        new_filename = self.args.destiny + \
                       self.args.output_prefix + \
                       original_filename.replace('.fits', '') + \
                       f_end

        #  print('spectrum[0]')
        # print(spectrum[0])
        # print('spectrum[1]')
        # print(spectrum[1])
        # print(len(spectrum))

        fits.writeto(new_filename, spectrum[1], new_header, clobber=True)
        log.info('Created new file: {:s}'.format(new_filename))
        # print new_header
        return new_header

    def display_onscreen_message(self, message='', color='red'):
        """Uses the fourth subplot to display a message

        Displays a warning message on the bottom right subplot of the
        interactive window. It is capable to break down the message in more
        than one line if necessary.

        Args:
            message (str): The message to be displayed
            color (str): Color name for the font's color

        """
        full_message = [message]
        if len(message) > 30:
            full_message = []
            split_message = message.split(' ')
            line_length = 0
            # new_line = ''
            e = 0
            for i in range(len(split_message)):
                # print(i, len(split_message))
                line_length += len(split_message[i]) + 1
                if line_length >= 30:
                    new_line = ' '.join(split_message[e:i])
                    # print(new_line, len(new_line))
                    full_message.append(new_line)
                    # new_line = ''
                    line_length = 0
                    e = i
                if i == len(split_message) - 1:
                    new_line = ' '.join(split_message[e:])
                    # print(new_line, len(new_line))
                    full_message.append(new_line)
        try:
            self.ax4.cla()
            self.ax4.relim()
            self.ax4.set_xticks([])
            self.ax4.set_yticks([])
            self.ax4.set_title('Message')
            for i in range(len(full_message)):
                self.ax4.text(0.05, 0.95 - i * 0.05,
                              full_message[i],
                              verticalalignment='top',
                              horizontalalignment='left',
                              transform=self.ax4.transAxes,
                              color=color,
                              fontsize=15)
            self.i_fig.canvas.draw()
        except:
            pass
        return

    def display_help_text(self):
        """Shows static text on the top right subplot

        This will print static help text on the top right subplot of the
        interactive window.

        Notes:
            This is really hard to format and having a proper GUI should help
            to have probably richer formatted text on the screen.

        """
        self.ax2.set_title('Help')
        self.ax2.set_xticks([])
        self.ax2.set_yticks([])

        self.ax2.text(1, 11, 'F1 or ?:',
                      fontsize=13)

        self.ax2.text(1.46, 11, 'Prints Help.',
                      fontsize=13)

        self.ax2.text(1, 10.5, 'F2 or f:',
                      fontsize=13)

        self.ax2.text(1.46, 10.5, 'Fit Wavelength Solution to points',
                      fontsize=13)

        self.ax2.text(1.46, 10, 'collected',
                      fontsize=13)

        self.ax2.text(1, 9.5, 'F3 or a:',
                      fontsize=13)

        self.ax2.text(1.46, 9.5, 'Find new lines, use when the solution',
                      fontsize=13)

        self.ax2.text(1.46, 9, 'is already decent.',
                      fontsize=13)

        self.ax2.text(1, 8.5, 'F4:',
                      fontsize=13)

        self.ax2.text(1.46, 8.5, 'Evaluate Solution',
                      fontsize=13)

        self.ax2.text(1, 8, 'F6 or l:',
                      fontsize=13)

        self.ax2.text(1.46, 8, 'Linearize Data',
                      fontsize=13)

        self.ax2.text(1, 7.5, 'd :',
                      fontsize=13)

        self.ax2.text(1.46, 7.5, 'Delete Closest Point',
                      fontsize=13)

        self.ax2.text(1, 7, 'Ctrl+d:',
                      fontsize=13)

        self.ax2.text(1.5, 7, 'Delete all recorded marks.',
                      fontsize=13)

        self.ax2.text(1, 6, 'Ctrl+z:',
                      fontsize=13)

        self.ax2.text(1.5, 6, 'Remove all automatic added points.',
                      fontsize=13)

        self.ax2.text(1.5, 5.5, 'Undo what F3 does.',
                      fontsize=13)

        self.ax2.text(1, 5, 'Middle Button Click:',
                      fontsize=13)

        self.ax2.text(1.46, 4.5, 'Finds and records line position',
                      fontsize=13)

        self.ax2.text(1, 4, 'Enter :',
                      fontsize=13)

        self.ax2.text(1.46, 4, 'Close Figure and apply wavelength',
                      fontsize=13)

        self.ax2.text(1.46, 3.5, 'solution',
                      fontsize=13)

        self.ax2.set_ylim((0, 12))
        self.ax2.set_xlim((0.95, 3.5))


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

    @staticmethod
    def set_spectral_features(header):
        """Creates dictionary that defines the instrument configuration

        Both Blue and Red Camera produce slightly different FITS headers being
        the red camera the one that provides more precise and better
        information. This method will recognize the camera and create the
        dictionary accordingly.

        Notes:
            As of August 2017 both headers are FITS compliant and contain the
            same keywords.

        Args:
            header (object): Instance of astropy.io.fits.header.Header.

        Returns:
            A dictionary that contains key information regarding the kind of
            spectroscopic data, basically related to the instrument
            configuration.
            The keywords it reads are: INSTCONF, GRATING, ROI, FILTER, FILTER2,
            SLIT, WAVMODE, CAM_ANG, GRT_ANG.

        """
        # TODO (simon): Use CAM_TARG and GRT_TARG instead of CAM_ANG and GRT_ANG
        if header is None:
            log.error('Header has not been parsed')
        else:
            try:
                log.debug('{:s} Camera'.format(header['INSTCONF']))
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

                # for key in dict.keys():
                # print(key, dict[key])
                return spectral_dict
            except KeyError:
                log.debug('(Old) Blue Camera')
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

                # for key in dict.keys():
                # print(key, dict[key])
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
            header (object): FITS header instance from astropy.io.fits.

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

                        log.info('Keyword: {:s} does not Match'.format(
                            key.upper()))

                        log.info('{:s} - Solution: {:s} - New '
                                 'Data: {:s}'.format(key.upper(),
                                                     self.spectral_dict[key],
                                                     new_dict[key]))

                        return False

                    elif key in ['cam_ang',  'grt_ang'] and \
                                    abs(new_dict[key] -
                                                self.spectral_dict[key]) > 1:

                        log.debug('Keyword: {:s} Lamp: {:s} Data: '
                                  '{:s}'.format(key,
                                                self.spectral_dict[key],
                                                new_dict[key]))

                        log.info('Solution belong to a different Instrument'
                                 'Configuration.')

                        return False
                    # else:
                    #     return True
                elif self.spectral_dict['camera'] == 'blue':
                    if key in ['grating',
                               'ccdsum',
                               'serial_bin',
                               'parallel_bin'] and \
                                    new_dict[key] != self.spectral_dict[key]:

                        log.debug('Keyword: {:s} does not Match', key.upper())

                        log.info('{:s} - Solution: {:s} - New Data: {:s}',
                                 key.upper(),
                                 self.spectral_dict[key],
                                 new_dict[key])

                        return False

                    elif key in ['cam_ang',
                                 'grt_ang'] and abs(
                                float(new_dict[key]) -
                                float(self.spectral_dict[key]))> 1:

                        log.debug('Keyword: {:s} Lamp: {:s} Data: {:s}',
                                  key,
                                  self.spectral_dict[key],
                                  new_dict[key])

                        log.info('Solution belong to a different Instrument '
                                 'Configuration.')

                        return False
                    # else:
                    #     return True
            return True
        else:
            log.error('Header has not been parsed')
            return False

    @staticmethod
    def set_solution_name(header):
        """Defines a name for the solution

        Using the header's information define a string that could be used as a
        keyword to identify a particular solution in the event that multiple
        solutions or instances of this class are stored somewhere/somehow.

        Args:
            header (object): FITS header instance from astropy.io.fits.

        Returns:
            Wavelength solution name.

        """
        name_text = ''
        grating = None
        # Grating part of the flat name
        if header['grating'] == '<NO GRATING>':
            try:
                if header['wavmode'] != 'Imaging':
                    name_text += '_nogrt'
            except KeyError:
                log.error('KeyError: Blue Camera')
        else:
            grating = header['grating'].split('_')[1]
            name_text += '_' + grating
            # Mode for the grating part of the flat name
            try:
                mode = header['wavmode'].split(' ')[1]
                name_text += '_' + mode.upper()
            except KeyError:
                log.error('KeyError: Blue Camera')
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

                    log.debug(center_wavelength)

                    name_text += '_' + \
                                 mode.upper() + \
                                 '_{:d}nm'.format(int(round(center_wavelength)))
                else:
                    # print(mode)
                    log.error('WAVMODE: {:s} not supported'.format(mode))

        # First filter wheel part of the flat name
        if header['filter'] != '<NO FILTER>':
            name_text += '_' + header['filter']
        # Second filter wheel part of the flat name
        if header['filter2'] != '<NO FILTER>':
            name_text += '_' + header['filter2']
        name_text += '_' + re.sub('[a-zA-Z" ]', '', header['slit'])

        wsolution_name = 'ws' + name_text

        log.debug(wsolution_name)
        return wsolution_name


if __name__ == '__main__':
    log.error('This can not be run on its own.')
