from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import glob
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from ccdproc import CCDData
from astropy import units as u
import numpy as np
import re
import cPickle as pickle

from goodman_spec.wsbuilder import ReadWavelengthSolution
from goodman_spec.linelist import ReferenceData


class CleanLineList(object):

    def __init__(self, reference_dir):
        self.cleaned_list = []


        reference_data = ReferenceData(reference_dir)

        file_list = glob.glob(os.path.join(reference_dir, '*fits'))

        for lfile in file_list:
            self.fig, self.ax = plt.subplots()
            ccd = CCDData.read(lfile, unit=u.adu)
            read_wavelength = ReadWavelengthSolution(ccd.header, ccd.data)
            wavelength, intensity = read_wavelength()

            self.line_list = reference_data.get_line_list_by_name(
                ccd.header['OBJECT'])

            file_name = lfile.split('/')[-1]
            pickle_file_name = re.sub('.fits', '_list.pkl', file_name)

            manager = plt.get_current_fig_manager()
            manager.window.showMaximized()
            if os.path.isfile(pickle_file_name):
                with open(pickle_file_name, 'rb') as pickled_file:
                    line_list = pickle.load(pickled_file)
                    for ref_line in line_list:
                        self.ax.axvline(ref_line, color='r', alpha=1)
            else:
                for ref_line in self.line_list:
                    self.ax.axvline(ref_line, color='r', alpha=.4)


            self.ax.set_title(file_name)
            self.ax.plot(wavelength, intensity, color='k')
            self.ax.set_xlabel('Wavelength (Angstrom)')
            self.ax.set_ylabel('Intensity (ADU)')
            self.ax.set_xlim((wavelength[0], wavelength[-1]))

            self.fig.canvas.mpl_connect('button_press_event', self.on_click)
            plt.show()

            if self.cleaned_list != []:

                with open(pickle_file_name, 'wb') as list_file:
                    pickle.dump(self.cleaned_list,
                                list_file,
                                protocol=pickle.HIGHEST_PROTOCOL)

    def on_click(self, event):
        print(event)
        if event.button == 2:
            print(event.xdata)
            differences = np.abs(self.line_list - event.xdata)
            min_index = np.argmin(differences)
            ref_value = self.line_list[min_index]
            if ref_value is not None:
                self.cleaned_list.append(ref_value)
                print(ref_value)
                self.ax.plot(ref_value, event.ydata, marker='o', color='b')
                self.fig.canvas.draw()




if __name__ == '__main__':
    this_file_path = os.path.dirname(os.path.abspath(__file__)).split('/')[:-1]
    this_file_path.append('spectroscopy')
    this_file_path.append('ref_comp')

    reference_dir = '/'.join(this_file_path)

    CleanLineList(reference_dir=reference_dir)




