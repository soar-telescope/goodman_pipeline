import glob
from ccdproc import CCDData
import astropy.units as u
import re
import matplotlib.pyplot as plt
import time
from threading import Thread


class LineMatcher(object):

    def __init__(self, search_pattern='ll_*fits'):
        self.file_list = glob.glob(search_pattern)
        self.fig = None
        self.ax = None
        self.arrow = None
        self.lines = []
        self.ccd = None
        self.lock_identify = False

    def __call__(self):
        for fits_file in self.file_list:
            print(fits_file)
            self.ccd = CCDData.read(fits_file, unit=u.adu)
            self._create_plot()
            self.ccd.write(fits_file, overwrite=True)

    def identify_matching_line(self):
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

        self._print_header()
        self.lock_identify = False
        plt.close('all')

    def _print_header(self):
        for key in self.ccd.header.keys():
            print("{:8s} = {:s}".format(key, str(self.ccd.header[key])))

    def _on_click(self, event):
        if not self.lock_identify:
            print("Threading")
            thread = Thread(target=self.identify_matching_line)
            thread.start()


    def _key_pressed(self, event):
        if event.key == 'enter':
            plt.close('all')

    def _add_arrow(self, x_loc, y_loc, value):
        if self.arrow is not None:
            try:
                self.arrow.remove()
                self.ax.relim()
            except:
                pass
        bbox_props = dict(boxstyle="rarrow,pad=0.3", fc="cyan", ec="b", lw=2)
        # self.arrow = self.ax.text(x_loc, y_loc, value,
        #                           ha="center",
        #                           va="center",
        #                           rotation=45,
        #                           size=15,
        #                           bbox=bbox_props)
        self.arrow, = self.ax.plot(x_loc, y_loc, ha='top', va='top', linestyle='None', marker='o', color='r')
        # plt.draw()
        print("drawing")
        self.fig.canvas.draw()

    def _create_plot(self):
        if self.fig is not None:
            self.ax.cla()
            self.ax.relim()
        self.fig = plt.figure(11)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_title("Lamp: {:s}\nGrating: {:s}\nSlit: {:s}".format(
            self.ccd.header['OBJECT'],
            self.ccd.header['WAVMODE'],
            self.ccd.header['SLIT']))
        line_key_list = self.ccd.header[r'GSP_P*']
        for key in line_key_list:
            if re.match(r'GSP_P\d{3}', key) is not None:
                self.ax.axvline(self.ccd.header[key], color='r', alpha=0.4)
                x_loc = self.ccd.header[key]
                y_loc = self.ccd.data[int(round(x_loc))]
                self.ax.text(x_loc, y_loc, x_loc, rotation=90, size=15)
        self.ax.plot(self.ccd.data)
        self.ax.set_xlabel("Pixels")
        self.ax.set_ylabel("Intensity (ADU)")
        self.fig.tight_layout()
        fig_manager = plt.get_current_fig_manager()
        fig_manager.window.showMaximized()
        self.fig.canvas.mpl_connect('button_press_event', self._on_click)
        # self.fig.canvas.mpl_connect('key_press_event', self._key_pressed)
        plt.show()



if __name__ == '__main__':
    line_matcher = LineMatcher()
    line_matcher()



