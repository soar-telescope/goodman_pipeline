from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ccdproc import CCDData
import argparse
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
# disables the s key event for saving.
# plt.rcParams['keymap.save'] = ''
import astropy.units as u
import glob
import sys
import re

from goodman_spec.wsbuilder import ReadWavelengthSolution


def get_args(arguments=None):
    parser = argparse.ArgumentParser(
        description="Plots image or spectrum")

    parser.add_argument('file',
                        action='store',
                        nargs='+',
                        help="File containing fits data.")
    args = parser.parse_args(args=arguments)

    return args


class DataPlotter(object):

    def __init__(self, args):
        self.args = args
        self.fig = None
        self.ax = None
        self.file = None

    def __call__(self, in_file, save=False):

        self.file = in_file
        self.fig, self.ax = plt.subplots()

        # read data and get its wavelength solution
        ccd = CCDData.read(self.file, unit=u.adu)
        wcs_reader = ReadWavelengthSolution(header=ccd.header,
                                            data=ccd.data)
        wavelength, intensity = wcs_reader()


        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
        plt.title('{:s}\n{:s}'.format(self.file, ccd.header['OBJECT']))

        self.ax.plot(wavelength, intensity, color='k', label='Data')
        self.ax.axvline(6562.8, color='r')
        self.ax.set_xlim((wavelength[0], wavelength[-1]))
        self.ax.set_ylabel('Intensity (ADU)')
        self.ax.set_xlabel('Wavelength (Angstrom)')

        plt.legend(loc='best')
        plt.subplots_adjust(left=0.05,
                            right=0.99,
                            top=0.96,
                            bottom=0.04,
                            hspace=0.17,
                            wspace=0.11)
        # plt.tight_layout()

        if not save:
            self.fig.canvas.mpl_connect('key_press_event', self.key_pressed)
            plt.show()
        # else:
        #     output = re.sub('.fits', '.png', self.file)
        #     plt.savefig(output, dpi=600)

    def key_pressed(self, event):
        # print(event.key)
        if event.key == 'q':
            plt.close(self.fig)
            sys.exit()
        # elif event.key == 's':
        #     self.__call__(in_file=self.file, save=True)
        elif event.key == 'n':
            plt.close(self.fig)





if __name__ == '__main__':
    args = get_args()

    if type(args.file) == list and len(args.file) > 1:
        file_list = args.file
    elif len(args.file) == 1:
        # print(args.file)
        file_list = glob.glob(args.file[0])
    #print(file_list)


    plotter = DataPlotter(args=args)
    for image in file_list:
        plotter(in_file=image)