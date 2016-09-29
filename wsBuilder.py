from specutils.wcs import specwcs
from astropy.io import fits
import numpy as np
import shlex
import sys
import logging as log
from astropy.modeling import models, fitting
# from specutils.io import read_fits
import tc_read_fits_copy
import astropy.units as u
import matplotlib.pyplot as plt

from specutils.wcs import specwcs
from specutils import Spectrum1D

log.basicConfig(level=log.DEBUG)

class WavelengthFitter:
    def __init__(self, model='chebyshev', degree=3):
        self.model_name = model
        self.degree = degree
        self.model = None
        self.model_fit = None
        self.model_constructor()

    def model_constructor(self):
        if self.model_name == 'chebyshev':
            self.model = models.Chebyshev1D(degree=self.degree)
            self.model_fit = fitting.LevMarLSQFitter()
        elif self.model_name == 'linear':
            self.model = models.Linear1D()
            self.model_fit = fitting.LinearLSQFitter()
        # return model

    def ws_fit(self, physical, wavelength):
        if self.model and self.model_fit is not None:
            return self.model_fit(self.model, physical, wavelength)
        else:
            log.error('Either model or model fitter were not constructed')



class ReadWavelengthSolution:

    def __init__(self, header, data, reference=''):
        self.header = header
        self.data = data
        self.reference_file_name = reference
        self.wat_wcs_dict = dict()
        self.wcs_dict = dict()

        self.wave_intens = []

    def non_linear_solution(self, dimension):
        header = self.header
        ctypen = header['CTYPE%s' % dimension]
        if ctypen == 'MULTISPE':
            wat_head = header['WAT%s*' % dimension]
            if len(wat_head) == 1:
                print('Get units')
                wat_array = wat_head[0].split(' ')
                for pair in wat_array:
                    split_pair = pair.split('=')
                    self.wat_wcs_dict[split_pair[0]] = split_pair[1]
                    # print(wat_head[0].split(' '))
            elif len(wat_head) > 1:
                wat_string = ''
                for key in wat_head:
                    wat_string += header[key]
                wat_array = shlex.split(wat_string.replace('=', ' '))
                if len(wat_array) % 2 == 0:
                    for i in range(0, len(wat_array), 2):
                        # if wat_array[i] not in self.wcs_dict.keys():
                        self.wat_wcs_dict[wat_array[i]] = wat_array[i + 1]
                        print wat_array[i], wat_array[i + 1]

        for key in self.wat_wcs_dict.keys():
            log.debug("%s -%s- %s", dimension, key, self.wat_wcs_dict[key])
        if 'spec1' in self.wat_wcs_dict.keys():
            print self.wat_wcs_dict['spec1']
            spec = self.wat_wcs_dict['spec1'].split()
            aperture = int(spec[0])
            beam = int(spec[1])
            disp_type = int(spec[2])
            disp_start = float(spec[3])
            disp_del_av = float(spec[4])
            pix_num = int(spec[5])
            dopp_fact = float(spec[6])
            aper_low = int(float(spec[7]))
            aper_high = int(float(spec[8]))
            weight = float(spec[9])
            zeropoint = float(spec[10])
            function_type = spec[11]
            order = int(float(spec[12]))
            min_pix_val = int(float(spec[13]))
            max_pix_val = int(float(spec[14]))
            # resto = list(spec[9:])
            params = [float(i) for i in spec[15:]]
            self.wcs_dict = {'aperture': aperture,
                             'beam': beam,
                             'dtype': disp_type,
                             'dstart': disp_start,
                             'avdelt': disp_del_av,
                             'pnum': pix_num,
                             'z': dopp_fact,
                             'alow': aper_low,
                             'ahigh': aper_high,
                             'weight': weight,
                             'zeropoint': zeropoint,
                             'ftype': function_type,
                             'order': order,
                             'pmin': min_pix_val,
                             'pmax': max_pix_val,
                             'fpar': params}

            function = ReadMathFunctions(self.wcs_dict)
            solution = function.get_solution()
            print solution
            print "%s"%params, len(params)
            wav1 = [4545.0519,
                    4579.3495,
                    4589.8978,
                    4609.5673,
                    4726.8683,
                    4735.9058,
                    4764.8646,
                    4806.0205,
                    4847.8095,
                    5495.8738,
                    5506.1128,
                    5558.702,
                    5572.5413,
                    5606.733,
                    5650.7043]
            # data = fits.getdata('/data/simon/data/soar/work/goodman/test/extraction-tests/cuhear600nonlinearli.fits')
            x = range(1, len(self.data) + 1)
            # x0 = range(len(self.data))
            print('x data', x[0], x[-1], len(self.data))
            plt.title(self.header['OBJECT'])
            plt.xlabel("%s (%s)" % (self.wat_wcs_dict['label'], self.wat_wcs_dict['units']))
            plt.plot(solution(x), self.data)
            # plt.plot(solution(x0), self.data, color='g')
            for line in wav1:
                plt.axvline(line, color='r')
            plt.show()
            print(spec)

    def linear_solution(self):
        crval = float(self.header['CRVAL1'])
        crpix = int(self.header['CRPIX1'])
        cdelt = float(self.header['CDELT1'])
        self.wcs_dict = {'crval': crval,
                         'crpix': crpix,
                         'cdelt': cdelt,
                         'dtype': 0}
        function = ReadMathFunctions(self.wcs_dict)
        solution = function.get_solution()
        # data = fits.getdata('/data/simon/data/soar/work/goodman/test/extraction-tests/CuHeAr_600.fits')
        x = range(1, len(self.data) + 1)

        # plt.xlabel("%s (%s)" % (self.wat_wcs_dict['label'], self.wat_wcs_dict['units']))
        self.wave_intens = [solution(x), self.data]
        """
        plt.title(self.header['OBJECT'])
        plt.plot(solution(x), self.data)
        print solution(crpix), crval

        wav1 = [4545.0519,
                4579.3495,
                4589.8978,
                4609.5673,
                4726.8683,
                4735.9058,
                4764.8646,
                4806.0205,
                4847.8095,
                5495.8738,
                5506.1128,
                5558.702,
                5572.5413,
                5606.733,
                5650.7043]
        for line in wav1:
            plt.axvline(line, color='r')
        plt.show()
        """
        return solution



    def get_wavelength_solution(self):
        header = self.header
        wcsdim = int(header['WCSDIM'])
        for dim in range(1, wcsdim + 1):
            ctypen = header['CTYPE%s' % dim]
            if ctypen == 'LINEAR':
                print('Linear Solution')
                # self.wcs_dict = {'dtype': 0}
                self.linear_solution()
            elif ctypen == 'MULTISPE':
                self.non_linear_solution(dim)
        return self.wave_intens


class ReadMathFunctions:

    def __init__(self, wcs_dict):
        self.wcs = wcs_dict
        self.solution = None
        if self.wcs['dtype'] == -1:
            self.none()
        elif self.wcs['dtype'] == 0:
            self.solution = self.linear_solution()
        elif self.wcs['dtype'] == 1:
            self.log_linear()
        elif self.wcs['dtype'] == 2:
            if self.wcs['ftype'] == '1':
                self.solution = self.chebyshev()
            elif self.wcs['ftype'] == '2':
                self.non_linear_legendre()
            elif self.wcs['ftype'] == '3':
                self.non_linear_cspline()
            elif self.wcs['ftype'] == '4':
                self.non_linear_lspline()
            elif self.wcs['ftype'] == '5':
                #pixel coordinates
                pass
            elif self.wcs['ftype'] == '6':
                # sampled coordinate array
                pass
            else:
                log.error('Not Implemented')
        else:
            log.error('Not Implemented')

    def none(self):
        return 0

    def linear_solution(self):
        intercept = self.wcs['crval'] - (self.wcs['crpix'] - 1) * self.wcs['cdelt']
        linear = models.Linear1D(slope=self.wcs['cdelt'], intercept=intercept)
        return linear

    def log_linear(self):
        return False

    def chebyshev(self):
        cheb = models.Chebyshev1D(degree=self.wcs['order'], domain=[self.wcs['pmin'], self.wcs['pmax']], )
        for p in range(self.wcs['order']):
            cheb.parameters[p] = self.wcs['fpar'][p]
        return cheb

    def non_linear_legendre(self):
        return False

    def non_linear_lspline(self):
        return False

    def non_linear_cspline(self):
        """Cubic Spline"""
        # cubic_spline = models.
        return True

    def get_solution(self):
        if self.solution != None:
            return self.solution
        else:
            log.error("The solution hasn't been found")
"""
def onclick(event):
    print 'click ', event.xdata, ' ', event.ydata, ' ', event.button

def interactive_solution(raw_file, reference_file):
    # ------- Reference -------
    ref_data = fits.getdata(reference_file)
    ref_header = fits.getheader(reference_file)
    reference = ReadWavelengthSolution(ref_header, ref_data)
    reference_solution = reference.get_wavelength_solution()
    # ------- RAW -------
    raw_data = fits.getdata(raw_file)
    raw_header = fits.getheader(raw_file)
    raw_pixel_axis = range(1, len(raw_data) + 1, 1)
    raw = ReadWavelengthSolution(raw_data, raw_header)
    # raw_solution = raw.get_wavelength_solution()
    # ------- Plots -------
    fig = plt.figure(1)
    manager = plt.get_current_fig_manager()
    manager.window.maximize()
    plt.subplot(211)
    plt.title('Raw Data')
    plt.xlabel('Pixels')
    plt.ylabel('Intensity (counts)')
    plt.plot(raw_pixel_axis, raw_data)
    plt.xlim((0, 4096))

    plt.subplot(212)
    # plt.xlim((3589, 4930))
    plt.title('Reference Data')
    plt.xlabel('Wavelength (Angstrom)')
    plt.ylabel('Intensity (counts)')
    plt.plot(reference_solution[0], reference_solution[1])
    plt.xlim((reference_solution[0][0], reference_solution[0][-1]))
    plt.tight_layout()
    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
"""


if __name__ == '__main__':
    # filename = '/data/simon/data/soar/work/goodman/test/extraction-tests/CuHeAr_600_nonlinear_2.fits'
    # filename = '/data/simon/data/soar/work/goodman/test/extraction-tests/CuHeAr_600.fits'
    # filename = '/data/simon/data/soar/work/goodman/test/extraction-tests/cuhear600nonlinearli.fits'
    raw_filename = '/data/simon/data/soar/work/goodman/test/extraction-tests/exfc_0047.SO2016A-019_0320.fits'
    ref_lamp = '/data/simon/data/soar/work/goodman/test/extraction-tests/cuhear_reference_noao.fits'
    # header = fits.getheader(raw_filename)
    # data = fits.getdata(raw_filename)
    # c = ReadWavelengthSolution(header, data, reference=ref_lamp)
    # c.interactive_solution()
    interactive_solution(raw_filename, ref_lamp)
    sys.exit(0)


# print(wat)
# spectrum = tc_read_fits_copy.read_goodman_non_linear_spectrum(filename=filename, dispersion_unit=u.angstrom)
# print(spectrum[0].wavelength)
# print(dir(spectrum[0]))
# print(dir(spectrum[1]))
# plt.plot(spectrum[0].wavelength, spectrum[0].data)
# plt.show()
