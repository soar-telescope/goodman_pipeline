# -*- coding: utf8 -*-
"""Set of class and method to handle wavelength solutions


"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import logging
import matplotlib.pyplot as plt
import shlex

from astropy.modeling import models, fitting, Model
from ccdproc import CCDData

# self.log.basicConfig(level=self.log.DEBUG)


class WCS(object):

    def __init__(self):
        # wavelength solution fitter variables
        self.log = logging.getLogger(__name__)
        self.model_name = None
        self.degree = None
        self.model = None
        self.model_fitter = None
        self.fitted_model = None

        # wavelength solution reader from header variables
        self.ccd = None
        self.data = None
        self.header = None
        self.wcs = None
        self.wat_wcs_dict = dict()
        self.wcs_dict = dict()
        self.wavelength_and_intensity = []

    def __call__(self, *args, **kwargs):
        self.log.error("Please use the method read(), fit() or get_model().")

    def fit(self, physical, wavelength, model_name='chebyshev', degree=3):
        """Fits a mathematical model

        Args:
            physical (list): List of line centers in pixel
            wavelength (list): List of line centers in Angstroms
            model_name (str): Name of the Mathematical model that needs to be
              created
            degree (int): Degree or order of the mathematical model (usually is
              some kind of polynomial).
        """
        self.model_name = model_name
        self.degree = degree
        self._model_constructor()
        self.fitted_model = self._fitter(physical=physical,
                                         wavelength=wavelength)
        return self.fitted_model

    # def read(self, ccd=None, header=None):
    def read(self, ccd=None):
        assert isinstance(ccd, CCDData)
        self.ccd = ccd
        self.header = ccd.header
        self.data = ccd.data
        self.wcs = ccd.wcs.wcs

        wcsdim = self.wcs.naxis
        for dim in range(1, wcsdim + 1):
            ctypen = self.wcs.ctype[0]
            if ctypen == 'LINEAR':
                self.log.info('Reading Linear Solution')
                # self.wcs_dict = {'dtype': 0}
                self.model = self._read_linear()
            elif ctypen == 'MULTISPE':
                self._read_non_linear(dim)
            else:
                raise NotImplementedError("CTYPE {:s} is "
                                          "not recognized".format(ctypen))
        return self.wavelength_and_intensity

    @staticmethod
    def write_gsp_wcs(ccd, model):
        """Writes a GSP-specific wavelength solution

        In an effort to easily write non-linear wavelength solutions into a fits
        header this method add a set of keywords that describes a pixel to
        angstrom relationship by means of using the astropy's modeling tools.

        GSP stands for Goodman Spectroscopic Pipeline.

        Notes:
            A limited amount of mathematical models are implemented on the read
            side. So you have to be careful what you write.

        Args:
            ccd (object): CCDData instance. Its header attribute will be
              modified
            model (object): astropy.modeling.Model instance.

        Returns:

        """
        assert isinstance(ccd, CCDData)
        assert isinstance(model, Model)

        ccd.header.set('GSP_FUNC', value=model.__class__.name,
                       comment="Mathematical model")
        ccd.header.set('GSP_ORDR', value=model.degree,
                       comment="Mathematical model order")
        ccd.header.set('GSP_NPIX', value=ccd.size,
                       comment="Number of Pixels")
        for i in range(model.degree + 1):
            ccd.header.set('GSP_C{:03d}'.format(i),
                           value=model.__getattr__('c{:d}'.format(i)).value,
                           comment="Value of parameter c{:d}".format(i))
        return ccd

    def read_gsp_wcs(self, ccd):
        """Read a GSP-specific wavelength solution

        Args:
            ccd (object): CCDData instance

        Returns:
            astropy.modeling.Model instance

        """
        assert isinstance(ccd, CCDData)
        self.model_name = ccd.header['GSP_FUNC']
        self.degree = ccd.header['GSP_ORDR']

        if self.model_name == 'Chebyshev1D':
            self.model = models.Chebyshev1D(degree=self.degree)
            for i in range(ccd.header['GSP_ORDR'] + 1):
                self.model.__getattr__('c{:d}'.format(i)).value = ccd.header[
                    'GSP_C{:03d}'.format(i)]
            self.wavelength_and_intensity = [
                self.model(range(ccd.header['GSP_NPIX'])), ccd.data]

            return self.model

    def _model_constructor(self):
        """Generates callable mathematical model

        It can do chebyshev and linear model only but is easy to implement
        others. Chebyshev 3rd degree is by default since provided the best
        results for Goodman data.
        """
        if self.model_name == 'chebyshev':
            self.model = models.Chebyshev1D(degree=self.degree)
            self.model_fitter = fitting.LevMarLSQFitter()
        elif self.model_name == 'linear':
            self.model = models.Linear1D()
            self.model_fitter = fitting.LinearLSQFitter()
        else:
            raise NotImplementedError("The model {:s} is "
                                      "not implemented".format(self.model_name))

    def _fitter(self, physical, wavelength):
        """Wavelength solution fit

        Takes a list of pixel values and its respective wavelength values to do
        a fit to the mathematical model defined with the class.

        Args:
            physical (list): Position of lines in pixel values.
            wavelength (list): Position of lines in angstrom values.

        Returns:
            fitted_model i.e. wavelength solution.

        """
        if self.model and self.model_fitter is not None:
            try:

                fitted_model = self.model_fitter(self.model,
                                                 physical,
                                                 wavelength)
                return fitted_model
            except TypeError as error:
                self.log.info('Unable to do fit, please add more data points.')
                self.log.error('TypeError: %s', error)
                return None
        else:
            self.log.error('Either model or model fitter were not constructed')
            return None

    # wavelength solution reader private methods.
    def _read_non_linear(self, dimension):
        """Non linear solutions reader

        Notes:
            Not all kind of non-linear solutions are implemented. Apparently is
            not fully implemented either.

        Args:
            dimension (int): Solutions can be multi-dimensionals, this method is
                called for each one of them.

        Returns:

        """
        # TODO (simon): Complete implementation.
        header = self.header
        ctypen = header['CTYPE{:d}'.format(dimension)]
        if ctypen == 'MULTISPE':
            # TODO (simon): What is the * (asterisc) doing here?.
            wat_head = header['WAT{:d}*'.format(dimension)]
            if len(wat_head) == 1:
                self.log.debug('Get units')
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

            self.log.debug("{:d} -{:s}- {:s}".format(dimension,
                                                key,
                                                self.wat_wcs_dict[key]))

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

            # This section of the code only shows a plot to see if the code
            # above actually worked which means this methods has not been fully
            # developed neither tested

            self._set_math_model()

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
            x_axis = range(1, len(self.data) + 1)
            # x0 = range(len(self.data))
            # print('x data', x_axis[0], x_axis[-1], len(self.data))
            plt.title(self.header['OBJECT'])

            plt.xlabel("%s (%s)" % (self.wat_wcs_dict['label'],
                                    self.wat_wcs_dict['units']))

            plt.plot(self.model(x_axis), self.data)
            # plt.plot(solution(x0), self.data, color='g')
            for line in wav1:
                plt.axvline(line, color='r')
            plt.show()
            # print(spec)

    def _read_linear(self):
        """Linear solution reader

        This method read the apropriate keywords and defines a linear wavelength
        solution

        Returns:
            Callable wavelength solution model. Instance of
                astropy.modeling.Model
        """
        self.wcs_dict = {'crval': self.wcs.crval[0],
                         'crpix': self.wcs.crpix[0],
                         'cdelt': self.wcs.cd[0],
                         'dtype': 0}

        self._set_math_model()

        x_axis = range(len(self.data))

        self.wavelength_and_intensity = [self.model(x_axis), self.data]

        return self.model

    def _set_math_model(self):

        if self.wcs_dict['dtype'] == -1:
            self._none()
        elif self.wcs_dict['dtype'] == 0:
            self._linear_solution()
        elif self.wcs_dict['dtype'] == 1:
            self._log_linear()
        elif self.wcs_dict['dtype'] == 2:
            if self.wcs_dict['ftype'] == '1':
                self._chebyshev()
            elif self.wcs_dict['ftype'] == '2':
                self._non_linear_legendre()
            elif self.wcs_dict['ftype'] == '3':
                self._non_linear_cspline()
            elif self.wcs_dict['ftype'] == '4':
                self._non_linear_lspline()
            elif self.wcs_dict['ftype'] == '5':
                # pixel coordinates
                raise NotImplementedError
            elif self.wcs_dict['ftype'] == '6':
                # sampled coordinate array
                raise NotImplementedError
            else:
                self.log.error('Not Implemented')
        else:
            self.log.error('Not Implemented')

    @staticmethod
    def _none():
        # TODO (simon): Document why this is this way
        return 0

    def _linear_solution(self):
        """Returns a linear 1D model"""
        intercept = self.wcs_dict['crval'] -\
                    (self.wcs_dict['crpix'] - 1) *\
                    self.wcs_dict['cdelt']

        self.model = models.Linear1D(slope=self.wcs_dict['cdelt'],
                                     intercept=intercept)

    @staticmethod
    def _log_linear():
        """Not implemented, returns False"""
        return False

    def _chebyshev(self):
        """Returns a chebyshev model"""
        self.model = models.Chebyshev1D(degree=self.wcs_dict['order'],
                                        domain=[self.wcs_dict['pmin'],
                                                self.wcs_dict['pmax']], )
        for param_index in range(self.wcs_dict['order']):
            self.model.parameters[param_index] = self.wcs_dict['fpar'][param_index]

    def _non_linear_legendre(self):
        """Not implemented

        Raises:
            NotImplementedError

        """
        raise NotImplementedError

    def _non_linear_lspline(self):
        """Not implemented

        Raises:
            NotImplementedError

        """
        raise NotImplementedError

    def _non_linear_cspline(self):
        """Not implemented

        Raises:
            NotImplementedError

        """
        raise NotImplementedError

    def get_model(self):
        """Returns the wavelength solution model if exists."""
        if self.model is not None:
            return self.model
        else:
            self.log.error("The solution hasn't been found")


if __name__ == '__main__':
    pass


