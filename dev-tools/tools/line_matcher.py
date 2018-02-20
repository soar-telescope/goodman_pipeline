import glob
from ccdproc import CCDData
import astropy.units as u
import re
import matplotlib.pyplot as plt
import time
from threading import Thread
import numpy as np
import sys

# from pipeline.spectroscopy.new_wavelength import WavelengthCalibration

"""Tool to record angstrom values for lines already recorded in pixels.

Usage:
  >>> ./line_matcher.py ll_*.fits
  
  A plot will show up. Clicking on the plot will trigger a function (which runs
  as a separated thread). Then you should go to the terminal an enter the
  corresponding Angstrom value for the line requested. For all the lines.
  
Notes:
  1. You can only enter values that can be converted to float.
  2. If you run the tool again on the same image will update the values.
  

"""


def get_spectral_characteristics(ccd):
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
    pixel_size = 15 * u.micrometer
    pixel_scale = 0.15 * u.arcsec
    goodman_focal_length = 377.3 * u.mm
    # TODO (simon): find a definite solution for this, this only work
    # TODO (simon): (a little) for one configuration
    blue_correction_factor = -50 * u.angstrom
    red_correction_factor = -37 * u.angstrom

    grating_frequency = float(
        re.sub('[A-Za-z_-]',
               '',
               ccd.header['GRATING'])) / u.mm

    # print('Grating Frequency ' +
    # '{:d}'.format(int(grating_frequency)))
    grating_angle = float(ccd.header['GRT_ANG']) * u.deg
    camera_angle = float(ccd.header['CAM_ANG']) * u.deg

    serial_binning, parallel_binning = [
        int(x) for x in ccd.header['CCDSUM'].split()]

    pixel_count = len(ccd.data)
    # Calculations
    # TODO (simon): Check whether is necessary to remove the
    # TODO (simon): slit_offset variable
    alpha = grating_angle.to(u.rad)
    beta = camera_angle.to(u.rad) - grating_angle.to(u.rad)

    center_wavelength = (np.sin(alpha) +
                              np.sin(beta)) / grating_frequency

    limit_angle = np.arctan(
        pixel_count *
        (pixel_size / goodman_focal_length) / 2)

    blue_limit = ((
                           np.sin(alpha) +
                           np.sin(beta - limit_angle.to(u.rad))) /
                       grating_frequency).to(u.angstrom) + \
                      blue_correction_factor

    red_limit = ((
                          np.sin(alpha) +
                          np.sin(beta +
                                 limit_angle.to(u.rad))) /
                      grating_frequency).to(u.angstrom) + \
                     red_correction_factor

    pixel_one = 0
    pixel_two = 0
    # log.debug(
    #     'Center Wavelength : {:.3f} Blue Limit : '
    #     '{:.3f} Red Limit : {:.3f} '.format(center_wavelength,
    #                                         blue_limit,
    #                                         red_limit))

    spectral_characteristics = {'center': center_wavelength,
                                'blue': blue_limit,
                                'red': red_limit,
                                'alpha': alpha,
                                'beta': beta,
                                'pix1': pixel_one,
                                'pix2': pixel_two}
    return spectral_characteristics


class LineMatcher(object):

    def __init__(self, search_pattern='ll_*fits'):
        """Tool created to record Angstrom values to matching pixel values

        Detected emission lines are recorded in the header of comparison lamps
        Also there is a placeholder for recording angstrom values.

        This tool shows a plot of the comparison lamp with the line values and
        ask you what are the angstrom values for each line recorded.

        :param search_pattern: Pattern search for `glob.glob()` to filter data
        :type search_pattern: str
        """
        self.file_list = glob.glob(search_pattern)
        self.fig = None
        self.ax = None
        self.arrow = None
        self.lines = []
        self.ccd = None
        self.lock_identify = False
        self.threads = []

    def __call__(self):
        """Run the tool for all the images matching `search_pattern`"""
        for fits_file in self.file_list:
            print(fits_file)
            self.ccd = CCDData.read(fits_file, unit=u.adu)
            if not self.threads:
                # plot_thread = Thread(target=self._create_plot)
                # plot_thread.start()
                id_thread = Thread(target=self.identify_matching_line)
                id_thread.start()
                # self.identify_matching_line()
                # self.threads.
                self._create_plot()
                id_thread.join()
            self.ccd.write(fits_file, overwrite=True)

    def identify_matching_line(self):
        """Interactive recording lines

        This function runs in an independent thread and is triggered by a click
        event on the plot.

        There is a very rudimentary locking system defined by a boolean
        `self.lock_identify`. This lock is activated at the beginning and
        deactivated at the end of the execution of this function.

        :return:
        """
        self.lock_identify = True
        line_key_list = self.ccd.header[r'GSP_P*']
        for key in line_key_list:
            if re.match(r'GSP_P\d{3}', key) is not None:
                angstrom_key = re.sub('_P', '_A', key)
                while True:
                    try:
                        angstrom_value = float(input("Enter the angstrom value "
                                                     "for the line at "
                                                     "{:f}({:f}):".format(
                            self.ccd.header[key],
                            self.ccd.header[angstrom_key])))
                    except ValueError:
                        print("Please enter a valid Angstrom Value")
                        continue
                    else:

                        self.ccd.header[angstrom_key] = angstrom_value
                        print(angstrom_key, angstrom_value)
                        break

        # self._print_header()
        self.lock_identify = False
        plt.close('all')

    def _print_header(self):
        """Prints keywords and values for a header, without the comments"""
        for key in self.ccd.header.keys():
            print("{:8s} = {:s}".format(key, str(self.ccd.header[key])))

    def _on_click(self, event):
        """Creates the thread for the line matching routine.

        The thread will run `self.identify_matching_lines` which in turn will
        activate the `self.lock_identify` lock. Therefore it can't create
        thread for an image already in process.

        :param event:
        :return:
        """
        pass
        # if not self.lock_identify:
        #     print("Threading")
        #     thread = Thread(target=self.identify_matching_line)
        #     thread.start()

    def _key_pressed(self, event):
        """Closes all plots when the `enter` key is pressed"""
        if event.key == 'enter':
            plt.close('all')

    # def _add_arrow(self, x_loc, y_loc, value):
    #     if self.arrow is not None:
    #         try:
    #             self.arrow.remove()
    #             self.ax.relim()
    #         except:
    #             pass
    #
    #     self.arrow, = self.ax.plot(x_loc, y_loc,
    #                                ha='top',
    #                                va='top',
    #                                linestyle='None',
    #                                marker='o',
    #                                color='r')
    #     # plt.draw()
    #     print("drawing")
    #     self.fig.canvas.draw()

    def _create_plot(self):
        """Creates reference plot of lamp and all registered lines.

        :return:
        """
        if self.fig is not None:
            self.ax.cla()
            self.ax.relim()
        self.fig = plt.figure(11)
        self.fig.canvas.set_window_title(self.ccd.header['GSP_FNAM'])
        self.ax = self.fig.add_subplot(111)
        self.ax.set_title("Lamp: {:s}\nGrating: {:s}\nSlit: {:s}\n\n\n".format(
            self.ccd.header['OBJECT'],
            self.ccd.header['WAVMODE'],
            self.ccd.header['SLIT']), color='b')

        sec_ax = self.ax.twiny()
        wave_char = get_spectral_characteristics(ccd=self.ccd)
        # print(wave_char)
        theoretical_angstrom_ax = np.linspace(wave_char['blue'],
                                              wave_char['red'],
                                              len(self.ccd.data))
        sec_ax.plot(theoretical_angstrom_ax, self.ccd.data, alpha=0)
        sec_ax.set_xlabel("Theoretical Angstrom Values")

        line_key_list = self.ccd.header[r'GSP_P*']
        for key in line_key_list:
            if re.match(r'GSP_P\d{3}', key) is not None:
                self.ax.axvline(self.ccd.header[key], color='r', alpha=0.4)
                x_loc = self.ccd.header[key]
                y_loc = self.ccd.data[int(round(x_loc))]
                self.ax.text(x_loc, y_loc, x_loc, rotation=90, size=13)
        self.ax.plot(self.ccd.data)
        self.ax.set_xlabel("Pixels")
        self.ax.set_ylabel("Intensity (ADU)")
        self.fig.tight_layout()
        fig_manager = plt.get_current_fig_manager()
        fig_manager.window.showMaximized()
        # self.fig.canvas.mpl_connect('button_press_event', self._on_click)
        # self.fig.canvas.mpl_connect('key_press_event', self._key_pressed)
        plt.show()


if __name__ == '__main__':
    try:
        line_matcher = LineMatcher(search_pattern=sys.argv[1])
    except IndexError:
        line_matcher = LineMatcher()
    line_matcher()



