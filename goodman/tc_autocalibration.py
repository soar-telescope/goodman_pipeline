import matplotlib.pyplot as plt
from astropy.modeling import models
import numpy as np
import pandas
import logging as log
import wsbuilder
from astropy.io import fits


class Experiment(object):

    def __init__(self):
        self.line_list_files = {'cu': 'Cu_3000A-10000A_clean.csv',
                                'he': 'He_3000A-10000A_clean.csv',
                                'ne': 'Ne_3000A-10000A_clean.csv',
                                'ar': 'Ar_3000A-10000A_clean.csv',
                                'hg': 'Hg_3000A-10000A_clean.csv'}
        self.args_reference_dir = '/user/simon/development/soar/goodman/refdata/'
        # self.file_name = '/user/simon/data/soar/work/goodman/test/etest2/gfzh.0043.SO2016A-019_0320.fits'
        self.file_name = '/user/simon/data/soar/work/goodman/test/extraction-tests/CuHeAr_600.fits'
        self.data = fits.getdata(self.file_name)
        self.header = fits.getheader(self.file_name)

    def get_ref_spectrum_from_linelist(self, blue, red, name):
        """Experimental

        Builds a unidimensional spectrum to be used as a template for finding an automatic wavelength solution
        """
        wsc = wsbuilder.ReadWavelengthSolution(self.header, self.data)
        sample_data = wsc.get_wavelength_solution()
        # wav_axis = np.linspace(blue, red, 5000)
        wav_axis = sample_data[0]
        blue = sample_data[0][0]
        red = sample_data[0][-1]
        print 'blue', blue, 'red', red
        int_axis = np.zeros(len(wav_axis))
        if len(name) % 2 == 0:
            elements = [name[i:i + 2].lower() for i in range(0, len(name), 2)]
            for element in elements:
                linelist_file = self.line_list_files[element]
                print linelist_file
                pdf = pandas.read_csv(self.args_reference_dir + linelist_file)
                filtered_pdf = pdf[(pdf.ion == '%s I' % element.title())
                                   & pandas.notnull(pdf.airwave)
                                   & pandas.notnull(pdf.relint)
                                   & (pdf.airwave > blue)
                                   & (pdf.airwave < red)]
                slice = filtered_pdf[['airwave', 'relint']]
                slice = slice.values.tolist()
                for value in slice:
                    wavelength = value[0]
                    rel_int = value[1]
                    print wavelength, rel_int
                    gaussian = models.Gaussian1D(amplitude=rel_int, mean=wavelength, stddev=4)
                    int_axis += gaussian(wav_axis)

            plt.title(self.header['OBJECT'])
            plt.ylim((np.min(sample_data[1]), np.max(sample_data[1]) + 1000))
            plt.xlim(sample_data[0][0], sample_data[0][-1])
            plt.plot(sample_data[0], sample_data[1], label='Real Spectrum')
            plt.plot(wav_axis, int_axis, label='Sintetic Spectrum (NIST)')
            plt.legend(loc='best')
            plt.savefig('/user/simon/data/soar/work/goodman/test/extraction-tests/automatic_attempt_%s.png' % self.header['OBJECT'], dpi=300)
            plt.show()


        else:
            log.error('Error in the calibration lamp name: %s', name)
            return None


if __name__ == '__main__':
    exp = Experiment()
    exp.get_ref_spectrum_from_linelist(3000., 6270., exp.header['OBJECT'])
    """
    (
                                   | (pdf.ion == '%s II' % element.title()))

    """