import argparse
import os
import sys

import astropy.units as u
import logging

from importlib.metadata import version
import matplotlib as mpl
import numpy as np
from astropy.nddata import CCDData
from matplotlib import pyplot as plt
from pandas import DataFrame

from goodman_pipeline.core import (get_lines_in_lamp, get_spectral_characteristics, NoMatchFound)
from goodman_pipeline.core import ReferenceData
from goodman_pipeline.wcs import WCS

mpl.use('QtAgg')
# mpl.rcParams['savefig.dpi'] = 300
__version__ = version('goodman_pipeline')

FIGURE_SIZES_FOR_SCREEN = {
    "small" : (10, 5),
    "medium" : (16, 8),
    "large" : (20, 10),
}

def get_args(arguments=None):
    log = logging.getLogger()

    parser = argparse.ArgumentParser(
        description="Creates reference lamp.\nPipeline Version: {:s}".format(__version__))
    parser.add_argument("--comparison-lamp", action="store", help="Comparison lamp file name")
    parser.add_argument("--reference-lamp", action="store", default=None, help="Already calibrated comparison lamp file name.")
    parser.add_argument("--plots-theme", action="store", default="dark", choices=["light", "dark"], help="Choose a theme for plotting, default is dark.")
    parser.add_argument("--screen-size", action="store", default='large', choices=['small', 'medium', 'large'], help="Choose a screen size for sizing the plots, default is large.")
    parser.add_argument("--debug", action="store_true", default=False, help="Enable debug mode.")
    parser.add_argument("-v", "--version", action="version", version=__version__)
    args = parser.parse_args(args=arguments)

    if not args.comparison_lamp:
        log.error("Comparison lamp file name not specified.")
        parser.print_help()
        sys.exit("Please specify a comparison lamp file name.")

    return args

class CreateReferenceLamp:

    def __init__(self):
        self.log = logging.getLogger(__name__)
        self.args = None
        self.pixel_size = 15 * u.micrometer
        self.instrument_focal_length = 377.3 * u.mm
        self.wcs = WCS()
        self.pixel = []
        self.angstrom = []
        self.reference_data = ReferenceData(reference_dir=os.path.join(os.path.dirname(__file__), '../data/ref_comp'))
        self.wavelength_solution = None
        self.comparison_lamp = None
        self.comparison_lines = []
        self.comp_xmin = None
        self.comp_xmax = None
        self.comp_ymin = None
        self.comp_ymax = None
        self.comp_spectral_characteristics = None
        self.reference_lamp = None
        self.reference_lines = []
        self.ref_xmin = None
        self.ref_xmax = None
        self.ref_ymin = None
        self.ref_ymax = None
        self.ref_wavelength = None
        self.ref_intensity = None
        self.fig = None
        self.ax_ref = None
        self.ax_comp = None
        self.comp_over_ref = None
        self.ref_markers = None
        self.comp_markers = None

        # recenter
        self.recenter_fig = None
        self.recenter_ax = None
        self.line_center = None
        self.recenter_plot = None
        self.recenter_center = None
        self.recenter_callback = None

    def __call__(self, args=None):
        if self.args is None:
            self.args = get_args(arguments=args)
        else:
            self.args = args

        if self.args.debug:
            reference_data_string = str(self.reference_data).split("\n")
            for line in reference_data_string:
                self.log.debug(line)

        if not os.path.exists(self.args.comparison_lamp):
            self.log.error("Comparison lamp file not found.")
            sys.exit("Please specify a comparison lamp file name.")

        self.comparison_lamp = CCDData.read(self.args.comparison_lamp, unit=u.adu)
        self.comparison_lines = get_lines_in_lamp(ccd=self.comparison_lamp, peak_percent_for_threshold=3)
        self.comp_spectral_characteristics = get_spectral_characteristics(
            ccd=self.comparison_lamp,
            pixel_size=self.pixel_size,
            instrument_focal_length=self.instrument_focal_length)

        self.ref_xmin = self.comp_spectral_characteristics['blue'].value
        self.ref_xmax = self.comp_spectral_characteristics['red'].value

        if self.args.reference_lamp and os.path.exists(self.args.reference_lamp):
            self.reference_lamp = CCDData.read(self.args.reference_lamp, unit=u.adu)
        else:
            reference_lamp_full_path = self._identify_best_reference_lamp()
            self.reference_lamp = CCDData.read(reference_lamp_full_path, unit=u.adu)

        self.ref_wavelength, self.ref_intensity = self.wcs.read_gsp_wcs(ccd=self.reference_lamp)

        if self.args.plots_theme == "dark":
            plt.style.use('dark_background')
        elif self.args.plots_theme == "light":
            plt.style.use('default')
        self.log.debug("Disabling matplotlib shortcut for full screen.")
        mpl.rcParams['keymap.fullscreen'] = []
        mpl.rcParams['keymap.home'] = []
        mpl.rcParams['keymap.yscale'] = []

        self.__make_main_plot()



    def __make_main_plot(self):
        fig, (ax1, ax2) = plt.subplots(
            nrows=2,
            ncols=1,
            figsize=FIGURE_SIZES_FOR_SCREEN[self.args.screen_size])

        self.fig = fig
        self.ax_ref = ax1
        self.ax_comp = ax2

        # Comparison lamp plot
        comp_min = self.comparison_lamp.data.min()
        comp_max = self.comparison_lamp.data.max()
        comp_range = comp_max - comp_min
        self.comp_ymin = comp_min - 0.05 * comp_range
        self.comp_ymax = comp_max + 0.3 * comp_range
        comp_x_edge = 10
        self.comp_xmin = -comp_x_edge
        comp_xmax = self.comparison_lamp.data.shape[0] + comp_x_edge
        self.ax_comp.set_ylim(self.comp_ymin, self.comp_ymax)
        self.ax_comp.set_xlim(self.comp_xmin, comp_xmax)
        self.ax_comp.plot(self.comparison_lamp.data, color='C0')
        self.ax_comp.set(xlabel='Pixels', ylabel='Intensity (ADU)', title='Comparison lamp')
        if self.comparison_lines:
            for line in self.comparison_lines:
                text = f"{line:.3f}"
                text_y_position = self.comp_ymin + 0.98 * (self.comp_ymax - self.comp_ymin)
                line_intensity = np.max(self.comparison_lamp.data[int(line) - 1: int(line) + 1])
                line_ymin = (line_intensity - self.comp_ymin) / (self.comp_ymax - self.comp_ymin)
                self.ax_comp.axvline(x=line, ymin=line_ymin + 0.04, ymax=np.max([line_ymin, 0.8]), alpha=0.5, linestyle=':')
                self.ax_comp.text(line, text_y_position, text, rotation=90, verticalalignment='top',
                                  horizontalalignment='center', clip_on=True)

        # Reference lamp plot
        ref_min_index = np.abs(self.ref_wavelength - self.ref_xmin).argmin()
        ref_max_index = np.abs(self.ref_wavelength - self.ref_xmax).argmin()
        ref_intensity_subsample = self.ref_intensity[ref_min_index:ref_max_index]

        ref_min = np.min(ref_intensity_subsample)
        ref_max = np.max(ref_intensity_subsample)
        ref_range = ref_max - ref_min
        self.ref_ymin = ref_min - 0.05 * ref_range
        self.ref_ymax = ref_max + 0.3 * ref_range

        self.ax_ref.plot(self.ref_wavelength, self.ref_intensity, label=f"Reference lamp {self.reference_lamp.header['OBJECT']}", color='C0')
        self.ax_ref.set_ylim(self.ref_ymin, self.ref_ymax)
        self.ax_ref.set_xlim(self.ref_xmin, self.ref_xmax)
        self.ax_ref.set(xlabel='Wavelength (Angstrom)', ylabel='Intensity (ADU)', title=f"Reference lamp - {self.reference_lamp.header['OBJECT']}")
        for angstrom_key in self.reference_lamp.header['GSP_A*']:
            if int(float(self.reference_lamp.header[angstrom_key])) != 0:
                self.reference_lines.append(float(self.reference_lamp.header[angstrom_key]))
        if len(self.reference_lines) > 0:
            for line in self.reference_lines:
                text = f"{line:.3f}"
                text_y_position = self.ref_ymin + 0.98 * (self.ref_ymax - self.ref_ymin)
                line_index = np.abs(self.ref_wavelength - line).argmin()
                line_intensity = np.max(self.reference_lamp.data[int(line_index) - 1: int(line_index) + 1])
                line_ymin = (line_intensity - self.ref_ymin) / (self.ref_ymax - self.ref_ymin)
                self.ax_ref.axvline(x=line, ymin=line_ymin, ymax=np.max([line_ymin, 0.8]), alpha=0.5, linestyle=':')
                self.ax_ref.text(line, text_y_position, text, rotation=90, verticalalignment='top', horizontalalignment='center', clip_on=True)

        self.fig.tight_layout()
        self.fig.canvas.mpl_connect('button_press_event', self._on_click)
        self.fig.canvas.mpl_connect('key_press_event', self._on_key_pressed)
        print(f"\nPress 'Control' and click the line you want to select. Press 'h' for help.\n")
        plt.show()



    def __print_selected_points(self):
        if len(self.pixel) > 0 or len(self.angstrom) > 0:
            print(f"\nSelected points:")
            print(f"{'Pixel':8}\tAngstrom")
            for i in range(max([len(self.pixel), len(self.angstrom)])):
                pixel = f"{self.pixel[i]:.3f}" if len(self.pixel) > i else '-' * 8
                angstrom = f"{self.angstrom[i]:.3f}" if len(self.angstrom) > i else '-' * 8
                print(f"{pixel:8}\t{angstrom}")

    @staticmethod
    def __print_main_help():
        print("Commands available:\n")
        print("\tl : Print list of recorded values.")
        print("\td : Remove nearby data point.")
        print("\th : Print this help message.")

    @staticmethod
    def __print_recenter_help():
        print("""\nPress 'left' or 'right' arrow to adjust the line center. Then press 'Enter' to confirm.\n""")

    @staticmethod
    def __report_click_position(click_position, units):
        print(f"\nThe clicked position is at: {click_position:.3f} {units}.\nPress 'Control' + click to mark a selection.\n")

    def __delete_data_point(self, event):

        def get_index_of_element_to_remove_from_array(point: float, input_array: list, tolerance:float, units: str):
            closes_point_index = np.argmin(input_array - point)
            closes_point = input_array[closes_point_index]
            if closes_point - event.xdata <= tolerance:
                return closes_point_index
            else:
                self.log.error(f"Unable to find a data point within {tolerance} {units}")
                return None

        if event.inaxes == self.ax_comp and len(self.pixel) > 0:
            idx = get_index_of_element_to_remove_from_array(point=event.xdata, input_array=self.pixel, tolerance=1, units='Pixels')
            if idx is not None:
                removed = self.pixel.pop(idx)
                self.log.info(f"Removed point {removed:.3f} ")
        elif event.inaxes == self.ax_ref and len(self.angstrom) > 0:
            idx = get_index_of_element_to_remove_from_array(point=event.xdata, input_array=self.angstrom, tolerance=1, units='Angstrom')
            if idx is not None:
                removed = self.angstrom.pop(idx)
                self.log.info(f"Removed point {removed:.3f} ")

    def _refine_line_center(self, center, xaxis, data):
        self.line_center = center
        center_index = np.abs(xaxis - center).argmin()
        self.recenter_fig, self.recenter_ax = plt.subplots()
        offset = 10
        xaxis_sample = xaxis[int(center_index - offset):int(center_index + offset)]
        intensity_sample = data[int(center_index - offset):int(center_index + offset)]
        self.recenter_step = 0.01 * (xaxis_sample[-1] - xaxis_sample[0])
        self.recenter_ax.set(title=f"Line center at {self.line_center}")
        self.recenter_plot, = self.recenter_ax.plot(xaxis_sample, intensity_sample, color='C0')
        self.recenter_center = self.recenter_ax.axvline(self.line_center, color='C3', linestyle='--')
        self.recenter_callback = self.recenter_fig.canvas.mpl_connect('key_press_event', self._on_key_pressed_for_recenter)
        plt.show(block=False)
        self.__print_recenter_help()
        self.recenter_fig.canvas.start_event_loop()
        return self.line_center

    def _on_click(self, event):
        if event.inaxes in [self.ax_comp, self.ax_ref]:
            if event.inaxes == self.ax_comp:
                if event.button in [1, 2, 3] and event.key == 'control':
                    self._refine_line_center(center=event.xdata, xaxis=range(self.comparison_lamp.shape[0]), data=self.comparison_lamp.data)
                    self.pixel.append(self.line_center)
                    print(f"Register data point at {self.line_center:.3f} pixels.")
                    if len(self.pixel) > len(self.angstrom):
                        print(f"Now find the corresponding line in the reference lamp.")
                elif event.button == 1:
                    self.__report_click_position(click_position=event.xdata, units="Pixels")
            if event.inaxes == self.ax_ref:
                if event.button in [1, 2, 3] and event.key == 'control':
                    self.line_center = self._refine_line_center(center=event.xdata, xaxis=self.ref_wavelength, data=self.ref_intensity)
                    self.angstrom.append(self.line_center)
                    print(f"Register data point at {self.line_center:.3f} Angstrom.")
                    if len(self.angstrom) > len(self.pixel):
                        print(f"Now find the corresponding line in the comparison lamp.")
                elif event.button == 1:
                    self.__report_click_position(click_position=event.xdata, units="Angstroms")

            if len(self.pixel) > 0 or len(self.angstrom) > 0:
                self._draw_markers()

            if len(self.angstrom) == len(self.pixel) and len(self.angstrom) > 4:
                self._fit_wavelength_solution()
                self._overplot_comparison_lamp_with_solution()

    def _on_key_pressed(self, event):
        if event.key == 'f':
            self.log.info("Fitting data to model")
            self._fit_wavelength_solution()
            self._overplot_comparison_lamp_with_solution()
        elif event.key == 'l':
            self.__print_selected_points()
        elif event.key == 'h':
            self.__print_main_help()
        elif event.key == 'd':
            self.__delete_data_point(event=event)
            self._draw_markers(delete=True)

    def _on_key_pressed_for_recenter(self, event):
        replot = False
        if event.key == 'left':
            self.line_center -= self.recenter_step
            replot = True
        if event.key == 'right':
            self.line_center += self.recenter_step
            replot = True
        if event.key == 'enter':
            self.recenter_fig.canvas.mpl_disconnect(self.recenter_callback)
            self.recenter_fig.canvas.stop_event_loop()
            plt.close(self.recenter_fig)
            return
        if replot:
            if self.recenter_center is not None:
                self.recenter_center.remove()
                self.recenter_ax.relim()
            self.recenter_ax.set(title=f"Line center at {self.line_center}")
            self.recenter_center = self.recenter_ax.axvline(self.line_center, color='C3', linestyle='--')
            self.recenter_fig.canvas.draw()

    def _fit_wavelength_solution(self):
        if len(self.angstrom) == len(self.pixel) and (len(self.angstrom) >= 4 or len(self.pixel) >= 4):
            self.wavelength_solution = self.wcs.fit(physical=self.pixel,
                                                    wavelength=self.angstrom,
                                                    model_name='chebyshev',
                                                    degree=3)
            # print(self.wavelength_solution)

    def _overplot_comparison_lamp_with_solution(self):
        if self.wavelength_solution is not None:
            if self.comp_over_ref is not None:
                self.comp_over_ref.remove()
                self.ax_ref.relim()
            data_normalized_to_reference_lamp = (self.comparison_lamp.data / np.max(self.comparison_lamp.data)) * np.max(self.reference_lamp.data)
            self.comp_over_ref, = self.ax_ref.plot(self.wavelength_solution(range(self.comparison_lamp.data.shape[0])), data_normalized_to_reference_lamp, c='C3')
            self.fig.canvas.draw()

    def _draw_markers(self):
        if len(self.angstrom) > 0:
            if self.ref_markers is not None:
                self.ref_markers.remove()
                self.ax_ref.relim()
            ref_markers_yaxis = [self.reference_lamp.data.min()] * len(self.angstrom)
            self.ref_markers, = self.ax_ref.plot(self.angstrom, ref_markers_yaxis, marker='^', markersize=5, color='C6', linestyle='None')
        if len(self.pixel) > 0:
            if self.comp_markers is not None:
                self.comp_markers.remove()
                self.ax_comp.relim()
            comp_markers_yaxis = [self.comparison_lamp.data.min()] * len(self.pixel)
            self.comp_markers, = self.ax_comp.plot(self.pixel, comp_markers_yaxis, marker='^', markersize=5, color='C6', linestyle='None')
        self.fig.canvas.draw()

    def __estimate_spectral_features_of_reference_lamps(self, reference_lamps: DataFrame):
        spectral_data = []
        for _file in reference_lamps['file'].to_list():
            lamp_full_path = os.path.normpath(os.path.join(self.reference_data.reference_dir, _file))
            self.log.debug(f"Lamp full file path: {lamp_full_path}")
            ccd = CCDData.read(lamp_full_path, units='adu')
            wavelength, _ = self.wcs.read_gsp_wcs(ccd=ccd)
            spectral_features = get_spectral_characteristics(ccd=ccd, pixel_size=self.pixel_size, instrument_focal_length=self.instrument_focal_length)
            spectral_data.append([_file, wavelength[0], wavelength[-1], spectral_features['blue'].value, spectral_features['red'].value, spectral_features['center'].to(u.angstrom).value])
        spectral_df = DataFrame(data=spectral_data, columns=['file', 'wavelength_start', 'wavelength_end', 'blue', 'red', 'center'])
        return spectral_df


    def _identify_best_reference_lamp(self):
        try:
            reference_lamps_df = self.reference_data.get_reference_lamps_by_lamp_status_keyword(header=self.comparison_lamp.header)
            if self.args.debug:
                print(reference_lamps_df[['file', 'lamp_hga', 'lamp_ne', 'lamp_ar', 'lamp_fe' , 'lamp_cu']].to_string(index=False))
            ref_spectral = self.__estimate_spectral_features_of_reference_lamps(reference_lamps=reference_lamps_df)
            comp_blue = self.comp_spectral_characteristics['blue'].value
            comp_red = self.comp_spectral_characteristics['red'].value
            comp_center = self.comp_spectral_characteristics['center'].to(u.angstrom).value

            compatible_reference_lamps = ref_spectral[
                (
                    (ref_spectral['wavelength_start'] <= comp_center) &
                    (comp_center <= ref_spectral['wavelength_end'])
                ) | (
                    (ref_spectral['blue'] <= comp_center) &
                    (comp_center <= ref_spectral['red'])
                )]
            if compatible_reference_lamps.empty:
                raise NoMatchFound("Unable to find a reference lamp so that the comparison lamp's spectral center fits within.")
            elif len(compatible_reference_lamps) == 1:
                reference_lamp_filename = compatible_reference_lamps.file.to_string(index=False).strip()
                return os.path.normpath(os.path.join(self.reference_data.reference_dir, reference_lamp_filename))
            else:
                self.log.info(f"Found {len(compatible_reference_lamps)} compatible reference lamps. Comparison lamp theoretical wavelength starts at {comp_blue:.3f} Angstrom and ends at {comp_red:.3f} Angstrom.")
                max_file_length = np.max([len(f) for f in compatible_reference_lamps['file'].to_list()])
                comp_range = comp_red - comp_blue
                compatible_reference_lamps['wavelength_range'] = (compatible_reference_lamps['wavelength_end'] - compatible_reference_lamps['wavelength_start'])
                compatible_reference_lamps['range_difference'] = (compatible_reference_lamps['wavelength_range'] - comp_range).abs()
                compatible_reference_lamps = compatible_reference_lamps.sort_values(by=['range_difference'], ascending=True)
                for index, row in compatible_reference_lamps.iterrows():
                    self.log.info(f"File: {row['file']:>{max_file_length}}\tWavelength Solution Start: {row['wavelength_start']:.3f} Angstrom\tWavelength End: {row['wavelength_end']:.3f} Angstrom.")
                idx = compatible_reference_lamps['range_difference'].idxmin()
                row = compatible_reference_lamps.loc[idx]
                self.log.info(f"Selected {row['file']} as best reference lamp. Wavelength Start: {row['wavelength_start']:.3f} Angstrom and Wavelength End: {row['wavelength_end']:.3f} Angstrom.")
                self.log.warning(f"If the matched lamp is not good try any other from the list above. Use --reference-lamp {self.reference_data.reference_dir}/<file name>")
                reference_lamp_filename = row['file']
                return os.path.normpath(os.path.join(self.reference_data.reference_dir, reference_lamp_filename))

        except NoMatchFound:
            self.log.info("Here is a detailed view of all available reference lamps that match at least one of the lamps:")
            reference_data_with_some_compatibility = self.reference_data.get_reference_lamps_with_some_lamps_matching(header=self.comparison_lamp.header)
            print(reference_data_with_some_compatibility[['file', 'lamp_hga', 'lamp_ne', 'lamp_ar', 'lamp_fe' , 'lamp_cu']].to_string())
            sys.exit(f"Please specify a valid reference lamp with --reference-lamp {self.reference_data.reference_dir}/<file name>")
