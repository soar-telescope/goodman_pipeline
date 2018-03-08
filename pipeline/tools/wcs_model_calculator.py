from glob import glob
from ccdproc import CCDData
import os
import re
import sys
import numpy as np
import pandas
from pipeline.wcs import WCS
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


class GSPWcsCalculator(object):

    def __init__(self):
        self.path = None
        self.files = None
        self.ccd = None
        self.pixel = None
        self.angstrom = None
        self.spectrum = None
        self.wcs = WCS()
        self.nist = None

    def __call__(self, path, *args, **kwargs):
        print("Working on: {:s}".format(path))
        pattern = kwargs.get('pattern', 'll*.fits')
        self.path = path
        self.files = glob(pathname=os.path.join(path, pattern))
        for image_file in self.files:
            print("Processing file: {:s}".format(image_file))
            self.ccd = CCDData.read(image_file, unit='adu')
            self._recover_lines()
            if not self._validate_lines():
                print("Please check lines for lamp {:s}".format(image_file))
            else:
                pixel = np.asarray(self.pixel, dtype=float)
                angstrom = np.asarray(self.angstrom, dtype=float)
                wcs_model = self.wcs.fit(physical=pixel,
                                         wavelength=angstrom)
                self.ccd = self.wcs.write_gsp_wcs(self.ccd, wcs_model)
                image_name = os.path.basename(image_file)
                new_name = os.path.join(self.path, image_name[24:])
                # print(new_name, image_name)
                self.ccd.write(new_name, overwrite=True)

                self._pdf_generator(model=wcs_model)

    def _recover_lines(self):
        self.pixel = []
        self.angstrom = []
        pixel_keys = self.ccd.header['GSP_P*']
        for pixel_key in pixel_keys:
            if re.match(r'GSP_P\d{3}', pixel_key) is not None:
                angstrom_key = re.sub('GSP_P', 'GSP_A', pixel_key)
                assert pixel_key[-3:] == angstrom_key[-3:]
                assert angstrom_key in self.ccd.header
                if int(self.ccd.header[angstrom_key]) != 0:
                    self.pixel.append(float(self.ccd.header[pixel_key]))
                    self.angstrom.append(float(self.ccd.header[angstrom_key]))
                else:
                    print("File: {:s}".format(self.ccd.header['GSP_FNAM']))
                    print("Ignoring keywords: {:s}={:f}, {:s}={:f}".format(pixel_key, self.ccd.header[pixel_key],
                          angstrom_key, self.ccd.header[angstrom_key]))
        # self.pixel = np.asarray(self.pixel, dtype=float)
        # self.angstrom = np.asarray(self.angstrom, dtype=float)

    def _validate_lines(self):
        """Calls all available validation methods

        Returns:
            True if none of the validation fails.
        """
        assert len(self.pixel) == len(self.angstrom)
        if not self._order_validation(self.pixel):
            return False
        if not self._order_validation(self.angstrom):
            return False
        # if not self._validate_line_existence():
        #     return False
        self._validate_line_existence()
        return True

    @staticmethod
    def _order_validation(lines_array):
        """Checks that the array of lines only increases."""
        previous = None
        for line_value in lines_array:
            # print(line_value)
            if previous is not None:
                try:
                    assert line_value > previous
                    previous = line_value
                except AssertionError:
                    print("Error: Line {:f} is not larger "
                          "than {:f}".format(line_value, previous))
                    return False
            else:
                previous = line_value
        return True

    def _load_nist_list(self, **kwargs):
        """Load all csv files from strong lines in nist."""
        nist_path = kwargs.get(
            'path',
            os.path.join(os.path.dirname(sys.modules['pipeline'].__file__),
                         'data/nist_list'))
        assert os.path.isdir(nist_path)
        nist_files = glob(os.path.join(nist_path, "*.txt"))
        for nist_file in nist_files:
            key = os.path.basename(nist_file)[22:-4]
            nist_data = pandas.read_csv(nist_file, names=['intensity',
                                                          'air_wavelength',
                                                          'spectrum',
                                                          'reference'])
            self.nist[key] = nist_data

    def _validate_line_existence(self):
        """Check if a line actually exists in any table of NIST
        Notes:
            It does not actually check NIST, it loads six csv tables
            from NIST's strong lines for the elements used in lamps.
            Ar,  Cu, Fe, He, Hg, Ne. It does not work perfect so far
            so the it is not checking existence actually but if it finds it
            it will get the value at "spectrum" column in NIST tables which
            correspond to the source of the line for instance Ne I, Ar II.
            """

        lamp_elements = []
        lamp_name = self.ccd.header['OBJECT']
        if len(lamp_name) % 2 == 0:
            for element_index in range(0, len(lamp_name), 2):
                element = lamp_name[element_index:element_index + 2]
                lamp_elements.append(element)

        if self.nist is None:
            self.nist = {}
            self._load_nist_list()

        self.spectrum = list(self.angstrom)
        for i in range(len(self.angstrom)):
            for element in lamp_elements:
                line_info = self.nist[element][
                    self.nist[element].air_wavelength == self.angstrom[i]]
                if line_info.empty:
                    # print(self.angstrom[i], 'no-info')
                    self.spectrum[i] = ''
                else:
                    self.spectrum[i] = line_info['spectrum'].to_string(index=False)
                    # print(self.angstrom[i], line_info['spectrum'].to_string(index=False))
            print(self.angstrom[i], self.spectrum[i])

    def _pdf_generator(self, model):
        """Creates a pdf file Using Reference lines."""
        file_name = self.ccd.header['GSP_FNAM']
        pdf_file_name = os.path.join(
            os.path.join(os.path.dirname(sys.modules['pipeline'].__file__),
                         '../docs/goodman_comp_pdf'),
            re.sub('.fits', '.pdf', file_name[24:]))

        with PdfPages(pdf_file_name) as pdf:
            plt.figure(1, (20, 10))
            plt.plot(model(range(len(self.ccd.data))), self.ccd.data, color='k')
            line_spacer = 0

            top_lim = 1.1 * self.ccd.data.max()
            bottom_lim = self.ccd.data.min() - 0.05 * self.ccd.data.max()
            plt.ylim((bottom_lim, top_lim))
            plt.xlim((model(0), model(len(self.ccd.data))))

            for i in range(len(self.angstrom)):
                line_str = "{:.4f} {:s}".format(self.angstrom[i], str(self.spectrum[i]))
                try:
                    if self.angstrom[i+1] - self.angstrom[i] < 25:
                        line_spacer = (25 - (self.angstrom[i+1] - self.angstrom[i]))
                    elif line_spacer > 0:
                        line_spacer *= -1
                    else:
                        line_spacer = 0
                    # print(line_spacer)
                except IndexError:
                    line_spacer = 0
                x_pos = self.angstrom[i]
                y_pos = np.max((self.ccd.data[int(np.floor(self.pixel[i]))],
                                self.ccd.data[int(np.ceil(self.pixel[i]))]))
                y_offset = 0.05 * self.ccd.data.max()
                # plt.plot(x_pos, y_pos, x_pos - line_spacer, y_pos + y_offset, color='r')
                # y_min = y_pos / self.ccd.data.max()
                # y_max = (y_pos + y_offset) / self.ccd.data.max()
                # plt.axvline(x=x_pos, ymin=y_min, ymax=y_max, color='r')
                plt.text(x_pos - line_spacer, y_pos + y_offset, line_str, rotation=90,
                         verticalalignment='bottom',
                         horizontalalignment='center')
                # plt.axvline(aline, color='r', alpha=0.5)

            title = "{:s} - {:s} - {:s}".format(self.ccd.header['OBJECT'],
                                                self.ccd.header['WAVMODE'],
                                                self.ccd.header['SLIT'])
            plt.title(title)
            plt.xlabel('Wavelength (Angstrom)')
            plt.ylabel('Intensity (ADU)')
            plt.tight_layout()
            pdf.savefig()
            # pdf.close()
            # plt.show()
            plt.close()



if __name__ == '__main__':
    calculator = GSPWcsCalculator()
    path = '/data/simon/data/soar/comp_lamp_lib/work/comparison_lamp_library/completed'
    path_david = '/user/simon/data/soar/work/2018-02-08_David/comp/completed'
    calculator(path=path_david)