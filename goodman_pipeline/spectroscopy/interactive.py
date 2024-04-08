from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import glob
import logging
import os
import re
import sys

import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import scipy.interpolate
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip
from ccdproc import CCDData
from matplotlib.backends.backend_pdf import PdfPages
from scipy import signal

from ..wcs.wcs import WCS
from ..core import (add_linear_wavelength_solution,
                    bin_reference_data,
                    cross_correlation,
                    evaluate_wavelength_solution,
                    get_lines_in_lamp,
                    get_spectral_characteristics,
                    linearize_spectrum,
                    record_wavelength_solution_evaluation,
                    write_fits)

from ..core import (ReferenceData, NoMatchFound)

FORMAT = '%(levelname)s:%(filename)s:%(module)s: 	%(message)s'
logging.basicConfig(level=logging.DEBUG, format=FORMAT)
log = logging.getLogger(__name__)



SHOW_PLOTS = False


class InteractiveWavelengthCalibration(object):
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
            args (object): Runtime arguments.

        """

        # TODO - Documentation missing
        self.poly_order = 2
        self.wcs = WCS()
        self.wsolution = None
        self.rms_error = None
        self.n_rejections = None
        self.n_points = None
        # print(self.args.reference_dir)
        self.reference_data = None
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
        self.line_list = []
        self.line_pixels = []
        self.line_angstroms = []
        self.ref_filling_value = 1000
        self.raw_filling_value = 1000
        self.events = True
        self.first = True
        self.evaluation_comment = None
        # self.binning = self.lamp.header[]
        self.pixel_center = []
        """this data must come parsed"""
        # self.science_pack = sci_pack
        # self.sci_filename = self.science_object.file_name
        # self.history_of_lamps_solutions = {}
        self.reference_solution = None
        self.linearize = False

    def __call__(self,
                 ccd,
                 comp_list,
                 save_data_to,
                 reference_data,
                 object_number=None,
                 output_prefix='w',
                 wsolution_obj=None,
                 linearize=False,
                 plot_results=False,
                 save_plots=False,
                 plots=False):
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

        self.linearize = linearize

        if os.path.isdir(reference_data):
            self.reference_data = ReferenceData(reference_dir=reference_data)
        else:
            log.error(f"Reference Data Directory provided  :{reference_data} does not exist.")

        self.i_fig = None

        log.info('Processing Science Target: '
                 '{:s}'.format(ccd.header['OBJECT']))
        if comp_list is not None:
            for self.lamp in comp_list:

                self.line_list = self.reference_data.get_line_list_by_lamp(ccd=self.lamp)
                print(self.line_list)
                try:
                    self.calibration_lamp = self.lamp.header['GSP_FNAM']
                except KeyError:
                    self.calibration_lamp = ''

                self.raw_pixel_axis = range(self.lamp.shape[0])
                # self.raw_pixel_axis = range(len(self.lamp.data))
                # self.lamp.header = lamp_ccd.header.copy()
                self.lamp_name = self.lamp.header['OBJECT']

                log.info('Processing Comparison Lamp: '
                         '{:s}'.format(self.lamp_name))

                # self.data1 = self.interpolate(self.lamp.data)
                # self.lines_limits = self.get_line_limits()
                # self.lines_center = self.get_line_centers(self.lines_limits)
                self.lines_center = get_lines_in_lamp(ccd=self.lamp)
                self.spectral = get_spectral_characteristics(
                    ccd=self.lamp,
                    pixel_size=self.pixel_size,
                    instrument_focal_length=self.goodman_focal_length)
                object_name = ccd.header['OBJECT']
                # try:
                self.interactive_wavelength_solution(
                    object_name=object_name)
                # except TypeError as error:
                #     log.error(error)

                if self.wsolution is not None:
                    # TODO (simon): plug in a record system
                    # record = '{:s} {:.3f} {:.3f}'.format(
                    #     self.lamp.header['GRATING'],
                    #     self.lamp.header['GRT_TARG'],
                    #     self.lamp.header['CAM_TARG'])
                    #
                    # for par in self.wsolution.parameters:
                    #     record += ' {:.5f}'.format(par)
                    #
                    # os.system("echo \'{:s}\' >> parametros.txt
                    # ".format(record))
                    # print(self.wsolution)
                    # self.lamp.write(os.path.join(self.args.destiny, 'non-linear-raw-1d.fits'))
                    # linear = self.lamp.copy()
                    self.save()

                    print(self.wsolution)

                    if plots or plot_results or save_plots:

                        plt.close(1)
                        if not plots:
                            plt.ion()
                            # plt.show()
                        else:
                            plt.ioff()

                        wavelength_axis = self.wsolution(range(ccd.data.size))

                        object_name = ccd.header['OBJECT']
                        grating = ccd.header['GRATING']

                        fig_title = 'Wavelength Calibrated Data : ' \
                                    '{:s}\n{:s}'.format(object_name, grating)

                        fig, ax1 = plt.subplots(1)
                        fig.canvas.set_window_title(ccd.header['GSP_FNAM'])
                        # ax1 = fig.add_subplot(111)
                        manager = plt.get_current_fig_manager()
                        if plt.get_backend() == u'GTK3Agg':
                            manager.window.maximize()
                        elif plt.get_backend() == u'Qt5Agg':
                            manager.window.showMaximized()

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
                            plots_dir = os.path.join(save_data_to, 'plots')
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

                        if plot_results:
                            plt.show()
                        else:
                            plt.draw()
                            plt.pause(2)
                            plt.ioff()

                    return self.wsolution
                else:
                    log.error('It was not possible to get a wavelength '
                              'solution from this lamp.')
                    return None
        else:
            print('Data should be saved anyways')

    def save(self, format='fits'):

        # self.lamp.header = self._a
        rms_error, n_points, n_rejections = self.evaluate_solution()

        header.set('GSP_WRMS', value=rms_error)
        header.set('GSP_WPOI', value=n_points)
        header.set('GSP_WREJ', value=n_rejections)

        if evaluation_comment is None:
            self.evaluation_comment = 'Lamp Solution RMSE = {:.3f} ' \
                                      'Npoints = {:d}, ' \
                                      'NRej = {:d}'.format(rms_error,
                                                           n_points,
                                                           n_rejections)
        if self.linearize:
            self.lamp = record_wavelength_solution_evaluation(ccd=self.lamp,
                                                              rms_error=self.rms_error,
                                                              n_points=self.n_points,
                                                              n_rejections=self.n_rejections)
            linearized_lamp_x_axis, self.lamp.data = linearize_spectrum(
                data=self.lamp.data,
                wavelength_solution=self.wsolution)

            self.lamp.header = self.add_linear_wavelength_solution(
                new_header=self.lamp.header,
                spectrum=self.linear_lamp[0],
                original_filename=self.calibration_lamp,
                index=object_number)

            ccd = record_wavelength_solution_evaluation(ccd=ccd,
                                                        rms_error=self.rms_error,
                                                        n_points=self.n_points,
                                                        n_rejections=self.n_rejections)
            linearized_ccd_x_axis, ccd.data = linearize_spectrum(
                data=ccd.data,
                wavelength_solution=self.wsolution)

            self.header = self.add_linear_wavelength_solution(
                new_header=ccd.header,
                spectrum=self.linearized_sci[0],
                original_filename=ccd.header['GSP_FNAM'],
                index=object_number)

        else:
            pass


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

        new_filename = save_data_to + \
                       output_prefix + \
                       original_filename.replace('.fits', '') + \
                       f_end

        new_header.set('GSP_FNAM', value=os.path.basename(new_filename))

        #  print('spectrum[0]')
        # print(spectrum[0])
        # print('spectrum[1]')
        # print(spectrum[1])
        # print(len(spectrum))

        ccd = CCDData(data=spectrum[1], header=new_header, unit=u.adu)
        ccd.write(new_filename, overwrite=True)
        # print(ccd.header['GSP_FNAM'])

        # fits.writeto(new_filename, spectrum[1], new_header, overwrite=True)
        log.info('Created new file: {:s}'.format(new_filename))
        # print new_header
        return new_header


    def create_reference_lamp(self):
        log.info("Saving as template")
        self.lamp = record_wavelength_solution_evaluation(ccd=self.lamp,
                                                          rms_error=self.rms_error,
                                                          n_points=self.n_points,
                                                          n_rejections=self.n_rejections)

        self.reference_data.create_reference_lamp(
            ccd=self.lamp,
            wavelength_solution=self.wsolution,
            lines_pixel=self.line_pixels,
            lines_angstrom=self.line_angstroms)

    def get_best_filling_value(self, data):
        """Find the best y-value to locate marks

        The autmatically added points will be placed at a fixed location in the
        y-axis. This value is calculated by doing a 2-sigma clipping with 5
        iterations. Then the masked out values are removed and the median is
        calculated.

        Args:
            data (array): Array of 1D data

        Returns:
            Median value of clipped data.

        """
        clipped_data = sigma_clip(data, sigma=2, maxiters=5)
        clean_data = clipped_data[~clipped_data.mask]
        log.debug("Found best filling value"
                       " at {:f}".format(np.median(clean_data)))
        return np.median(clean_data)

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
                new_center.append(line)
            else:
                new_center.append(centroid)
        if plots:
            fig = plt.figure(9)
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

    @staticmethod
    def recenter_broad_lines(lamp_data, lines, order):
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
                abs(self.line_list - x_data))

            reference_line_value = self.line_list[reference_line_index]

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
            sub_y = self.lamp.data[pseudo_center - 10: pseudo_center + 10]
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
            This method uses the Qt5Agg backend, in theory it could also work
            with GTK3Agg but it is being forced to use Qt5Agg.

        """
        log.debug("Starting interactive wavelength calibration")
        plt.switch_backend('Qt5Agg')

        # disable full screen to allow the use of f for fitting the solution

        plt.rcParams['keymap.fullscreen'] = [u'ctrl+f']

        try:
            reference_lamp = self.reference_data.get_reference_lamp(
                header=self.lamp.header)
        except NotImplementedError:
            reference_lamp = self.reference_data.get_reference_lamps_by_name(
                lamp_name=self.lamp.header['OBJECT'])
            log.warning('Could not find a perfect match for reference '
                             'data')
            # reference_file = None
            # log.critical('Could not find a comparison lamp in the
            # reference.')

        # reference_file = self.reference_data.get_reference_lamps_by_name(
        #     self.lamp_name)

        if reference_lamp is not None:
            log.info('Using reference file: {:s}'.format(reference_lamp.header['GSP_FNAM']))
            reference_plots_enabled = True

            self.reference_solution = self.wcs.read_gsp_wcs(ccd=reference_lamp)
            print(self.reference_solution)
        else:
            reference_plots_enabled = False
            log.error('Please Check the OBJECT Keyword of your reference '
                           'data')

        # update filling value
        self.raw_filling_value = self.get_best_filling_value(
            data=self.lamp.data)

        if reference_lamp is not None:
            self.ref_filling_value = self.get_best_filling_value(data=reference_lamp.data)

        # ------- Plots -------
        self.i_fig, ((self.ax1, self.ax2), (self.ax3, self.ax4)) = \
            plt.subplots(2,
                         2,
                         gridspec_kw={'width_ratios': [4, 1]})

        self.i_fig.canvas.setWindowTitle('Science Target: {:s}'.format(
            object_name))

        manager = plt.get_current_fig_manager()
        if plt.get_backend() == u'GTK3Agg':
            manager.window.maximize()
        elif plt.get_backend() == u'Qt5Agg':
            manager.window.showMaximized()

        self.ax1.set_title('Raw Data - {:s}\n{:s} - {:s}'.format(
            self.lamp_name,
            self.lamp.header['GRATING'],
            self.lamp.header['SLIT']))

        self.ax1.set_xlabel('Pixels')
        self.ax1.set_ylabel('Intensity (counts)')
        self.ax1.plot([], linestyle='--', color='r', label='Detected Lines')
        for idline in self.lines_center:
            self.ax1.axvline(idline, linestyle='--', color='r')

        self.ax1.plot(self.raw_pixel_axis,
                      self.lamp.data,
                      color='k',
                      label='Raw Data')

        self.ax1.set_xlim((0, len(self.lamp.data)))
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
        self.ax3.set_xlim((self.spectral['blue'].value, self.spectral['red'].value))
        self.ax3.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

        self.ax3.plot([],
                      linestyle=':',
                      color='r',
                      label='Reference Line Values')

        for rline in self.line_list:
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

    def on_click(self, event):
        """Handles Click events for Interactive Mode

        Calls the method register_mark

        Args:
            event (object): Click event
        """
        if event.button == 3:
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
            print("ctrl+t: Save as template")
            print("ctrl+z: Go back to previous solution "
                  "(deletes automatic added points")
            print('Right Button Click or p: records data location.')
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
            else:
                log.debug('Wavelength solution is None')
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
                self.linearize_spectrum(self.lamp.data, plots=True)

        elif event.key == 'ctrl+z':
            log.info('Deleting automatic added points. If exist.')

            if self.raw_data_marks_x is not [] and \
                            self.reference_marks_x is not []:

                to_remove = []
                for i in range(len(self.raw_data_marks_x)):
                    # print self.raw_data_marks[i], self.filling_value
                    if self.raw_data_marks_y[i] == self.raw_filling_value:
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
        elif event.key == 'ctrl+t':
            self.create_reference_lamp()

        elif event.key == 'p':
            self.register_mark(event)

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
                    abs(self.line_list - wlines[i]))

                rline = self.line_list[closer_index]

                rw_difference = wlines[i] - rline
                # print('Difference w - r ', rw_difference, rline)
                square_differences.append(rw_difference ** 2)
                new_physical.append(self.lines_center[i])
                new_wavelength.append(rline)
            clipped_differences = sigma_clip(square_differences,
                                             sigma=2,
                                             maxiters=3)

            if len(new_wavelength) == len(new_physical) == \
                    len(clipped_differences):

                for i in range(len(new_wavelength)):
                    if clipped_differences[i] is not \
                            np.ma.masked and new_wavelength[i] not in \
                            self.reference_marks_x:

                        self.reference_marks_x.append(new_wavelength[i])
                        self.reference_marks_y.append(self.ref_filling_value)
                        self.raw_data_marks_x.append(new_physical[i])
                        self.raw_data_marks_y.append(self.raw_filling_value)
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
                    log.debug('Removing reference marks')
                except:
                    log.debug('Reference points is None')
                    pass
            log.debug("Plot new marks")
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
                    self.lamp.data,
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
                    abs(self.line_list - wline))

                rline = self.line_list[closer_index]

                rw_difference = wline - rline
                # print 'Difference w - r ', rw_difference, rline
                differences = np.append(differences, rw_difference)

            clipping_sigma = 2.
            # print(differences)
            # clipped_differences = sigma_clip(differences,
            #                                  sigma=clipping_sigma,
            #                                  iters=5,
            #                                  cenfunc=np.ma.median)
            #
            # once_clipped_differences = sigma_clip(differences,
            #                                       sigma=clipping_sigma,
            #                                       iters=1,
            #                                       cenfunc=np.ma.median)
            clipped_differences = differences
            once_clipped_differences = differences

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
                self.line_pixels = []
                self.line_angstroms = []
                for i in range(len(self.reference_marks_x)):
                    self.line_pixels.append(self.raw_data_marks_x[i])
                    self.line_angstroms.append(self.reference_marks_x[i])

                self.wsolution = self.wcs.fit(physical=self.line_pixels,
                                              wavelength=self.line_angstroms,
                                              model_name='chebyshev',
                                              degree=self.poly_order)

                self.evaluate_solution(plots=True)

        else:
            log.error('Clicks record is empty')
            self.display_onscreen_message(message='Clicks record is empty')
            if self.wsolution is not None:
                self.wsolution = None


    def add_gsp_wcs(self):
        pass

    def add_linear_wavelength_solution(self,
                                header,
                                spectrum,
                                original_filename,
                                evaluation_comment=None):
        """Add wavelength solution to the new FITS header

        Defines FITS header keyword values that will represent the wavelength
        solution in the header so that the image can be read in any other
        astronomical tool. (e.g. IRAF)

        Notes:
            This method also saves the data to a new FITS file, This should be
            in separated methods to have more control on either process.

        Args:
            header (object): An Astropy header object
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

        new_crpix = 1
        new_crval = spectrum[new_crpix - 1]
        new_cdelt = spectrum[new_crpix] - spectrum[new_crpix - 1]

        header['BANDID1'] = 'spectrum - background none, weights none, ' \
                                'clean no'
        header['WCSDIM'] = 1
        header['CTYPE1'] = 'LINEAR  '
        header['CRVAL1'] = new_crval
        header['CRPIX1'] = new_crpix
        header['CDELT1'] = new_cdelt
        header['CD1_1'] = new_cdelt
        header['LTM1_1'] = 1.
        header['WAT0_001'] = 'system=equispec'
        header['WAT1_001'] = 'wtype=linear label=Wavelength units=angstroms'
        header['DC-FLAG'] = 0
        header['DCLOG1'] = 'REFSPEC1 = {:s}'.format(self.calibration_lamp)
        return header

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

        self.ax2.text(1.46, 7, 'Delete all recorded marks.',
                      fontsize=13)

        self.ax2.text(1, 6.5, 'Ctrl+t:',
                      fontsize=13)

        self.ax2.text(1.46, 6.5, 'Save as reference lamp template.',
                      fontsize=13)

        self.ax2.text(1, 6, 'Ctrl+z:',
                      fontsize=13)

        self.ax2.text(1.46, 6, 'Remove all automatic added points.',
                      fontsize=13)

        self.ax2.text(1.46, 5.5, 'Undo what F3 does.',
                      fontsize=13)

        self.ax2.text(1, 5, 'Middle Button or p:',
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
