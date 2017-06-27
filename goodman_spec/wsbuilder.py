# -*- coding: utf8 -*-
"""Set of class and method to handle wavelength solutions


"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging
import shlex

from astropy.modeling import models, fitting
from matplotlib import pyplot as plt


# log.basicConfig(level=log.DEBUG)
log = logging.getLogger('redspec.wsbuilder')


class WavelengthFitter(object):
    """Contains methods to do pixel to angstrom fit"""

    def __init__(self, model='chebyshev', degree=3):
        """Initializes the class

        Args:
            model (str): Name of the model to fit, Chebyshev (default) or Linear
            degree (int): Degree of the model. Only needed by Chebyshev model.
        """
        self.model_name = model
        self.degree = degree
        self.model = None
        self.model_fit = None
        self.model_constructor()

    def model_constructor(self):
        """Generates callable mathematical model

        It can do chebyshev and linear model only but is easy to implement others.
        Chebyshev 3rd degree is by default since provided the best results for Goodman data.
        """
        if self.model_name == 'chebyshev':
            self.model = models.Chebyshev1D(degree=self.degree)
            self.model_fit = fitting.LevMarLSQFitter()
        elif self.model_name == 'linear':
            self.model = models.Linear1D()
            self.model_fit = fitting.LinearLSQFitter()
        # return model

    def ws_fit(self, physical, wavelength):
        """Wavelength solution fit

        Takes a list of pixel values and its respective wavelength values to do a fit to
        the mathematical model defined with the class.

        Returns:
            fitted_model (class): Fitted model
        """
        if self.model and self.model_fit is not None:
            try:
                # print(physical)
                # print(wavelength)
                fitted_model = self.model_fit(self.model, physical, wavelength)
                return fitted_model
            except TypeError as error:
                log.info('Please add more data points.')
                log.error('TypeError: %s', error)
        else:
            log.error('Either model or model fitter were not constructed')


class ReadWavelengthSolution(object):
    """Read wavelength solutions from a fits header"""

    def __init__(self, header, data, reference=''):
        self.header = header
        self.data = data
        self.reference_file_name = reference
        self.wat_wcs_dict = dict()
        self.wcs_dict = dict()
        self.wave_intens = []
        self.math_model = None

    def __call__(self):
        """call method

        Discriminates from header's keywords what kind of solution is present
        and call the appropriate method for linear or non-linear solutions.

        Returns:
            self.wave_intens (list): Returns the class attribute which is a two dimension list
            the first element is the x_axis in angstrom and the second is the y_axis in intensity
        """
        wcsdim = int(self.header['WCSDIM'])
        for dim in range(1, wcsdim + 1):
            ctypen = self.header['CTYPE%s' % dim]
            if ctypen == 'LINEAR':
                log.info('Reading Linear Solution')
                # self.wcs_dict = {'dtype': 0}
                self.math_model = self.linear_solution()
            elif ctypen == 'MULTISPE':
                self.non_linear_solution(dim)
        return self.wave_intens

    def non_linear_solution(self, dimension):
        """Non linear solutions reader

        Notes:
            Not all kind of non-linear solutions are implemented

        Args:
            dimension (int): Solutions can be multi-dimensionals, this method is called for each one of them.

        Returns:

        """
        header = self.header
        ctypen = header['CTYPE%s' % dimension]
        if ctypen == 'MULTISPE':
            wat_head = header['WAT%s*' % dimension]
            if len(wat_head) == 1:
                log.debug('Get units')
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
                        # print(wat_array[i], wat_array[i + 1])

        for key in self.wat_wcs_dict.keys():
            log.debug("%s -%s- %s", dimension, key, self.wat_wcs_dict[key])
        if 'spec1' in self.wat_wcs_dict.keys():
            # print(self.wat_wcs_dict['spec1'])
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
            # print(solution)
            # print("%s %s" % (params, len(params)))
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
            x_axis = range(1, len(self.data) + 1)
            # x0 = range(len(self.data))
            # print('x data', x_axis[0], x_axis[-1], len(self.data))
            plt.title(self.header['OBJECT'])
            plt.xlabel("%s (%s)" % (self.wat_wcs_dict['label'], self.wat_wcs_dict['units']))
            plt.plot(solution(x_axis), self.data)
            # plt.plot(solution(x0), self.data, color='g')
            for line in wav1:
                plt.axvline(line, color='r')
            plt.show()
            # print(spec)

    def linear_solution(self):
        """Linear solution reader

        This method read the apropriate keywords and defines a linear wavelength solution

        Returns:
            solution (class): Callable wavelength solution model
        """
        crval = float(self.header['CRVAL1'])
        crpix = int(self.header['CRPIX1'])
        cdelt = float(self.header['CDELT1'])
        self.wcs_dict = {'crval': crval,
                         'crpix': crpix,
                         'cdelt': cdelt,
                         'dtype': 0}
        self.function = ReadMathFunctions(self.wcs_dict)
        solution = self.function.get_solution()
        # data = fits.getdata('/data/simon/data/soar/work/goodman/test/extraction-tests/CuHeAr_600.fits')
        # x_axis = range(1, len(self.data) + 1)
        x_axis = range(len(self.data))

        # plt.xlabel("%s (%s)" % (self.wat_wcs_dict['label'], self.wat_wcs_dict['units']))
        self.wave_intens = [solution(x_axis), self.data]
        return solution



class ReadMathFunctions(object):

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
                # pixel coordinates
                pass
            elif self.wcs['ftype'] == '6':
                # sampled coordinate array
                pass
            else:
                log.error('Not Implemented')
        else:
            log.error('Not Implemented')

    @staticmethod
    def none():
        return 0

    def linear_solution(self):
        intercept = self.wcs['crval'] - (self.wcs['crpix'] - 1) * self.wcs['cdelt']
        # intercept = self.wcs['crval'] - self.wcs['crpix'] * self.wcs['cdelt']
        linear = models.Linear1D(slope=self.wcs['cdelt'], intercept=intercept)
        return linear

    @staticmethod
    def log_linear():
        return False

    def chebyshev(self):
        cheb = models.Chebyshev1D(degree=self.wcs['order'], domain=[self.wcs['pmin'], self.wcs['pmax']], )
        for param_index in range(self.wcs['order']):
            cheb.parameters[param_index] = self.wcs['fpar'][param_index]
        return cheb

    def non_linear_legendre(self):
        raise NotImplementedError

    def non_linear_lspline(self):
        raise NotImplementedError

    def non_linear_cspline(self):
        """Cubic Spline"""
        # cubic_spline = models.
        raise NotImplementedError

    def get_solution(self):
        if self.solution is not None:
            return self.solution
        else:
            log.error("The solution hasn't been found")

if __name__ == '__main__':
    pass
