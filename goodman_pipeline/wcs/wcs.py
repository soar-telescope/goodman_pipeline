# -*- coding: utf8 -*-
"""Set of class and method to handle wavelength solutions


"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import logging
import shlex
import sys

from astropy.modeling import models, fitting, Model
from ccdproc import CCDData

# self.log.basicConfig(level=self.log.DEBUG)


class WCS(object):
    """World Coordinate System class for Spectroscopy

    This class is intended to contain methods for all operations regarding
    WCS for spectroscopy or wavelength solution operations. Starting on the
    fitting to writing (and reading) from a FITS header.
    """

    def __init__(self):
        # wavelength solution fitter variables
        self.log = logging.getLogger(__name__)
        self.model_name = None
        self.degree = None
        self.model = None
        self.model_fitter = None
        self.fitted_model = None
        # dispersion direction binning
        self._binning = 1

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
        sys.exit(1)

    @property
    def binning(self):
        return self._binning

    @binning.setter
    def binning(self, value):
        if value > self._binning:
            self._binning = value
            self._set_binning_in_model()
        else:
            raise NotImplementedError("Impossible to perform")

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
        """Read WCS from FITS header

        Notes:
            The mathematical model stays as an attribute of the class as `model`

        Args:
            ccd (CCDData) Instance of `:class:`~astropy.nddata.CCDData``
              with FITS's wavelength solution.

        Returns:
            A list with an array representing the wavelength axis and another
            representing the intensity (ccd.data).

        """
        assert isinstance(ccd, CCDData)
        self.ccd = ccd
        self.header = ccd.header
        self.data = ccd.data
        self.wcs = ccd.wcs.wcs

        wcsdim = self.header['WCSDIM']
        for dim in range(1, wcsdim + 1):
            ctypen = self.wcs.ctype[dim - 1]
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

    def write_fits_wcs(self, ccd, model):
        """Write FITS WCS to the header

        Notes:
            This method is not implemented, in the current version the
            equivalent method resides within
            `goodman.pipeline.spectroscopy.wavelength.py`

        Args:
            ccd (CCDData) Instance of `:class:`~astropy.nddata.CCDData``
            model (object): Instance of `astropy.modeling.Model` that should be
              the mathematical representation of the wavelength solution of
              `ccd`

        Raises:
            NotImplementedError
        """
        raise NotImplementedError

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
            ccd (CCDData) :class:`~astropy.nddata.CCDData` instance. Its header
              attribute will be modified
            model (object): astropy.modeling.Model instance.

        Returns:

        """
        assert isinstance(ccd, CCDData)
        assert isinstance(model, Model)

        ccd.header.set('GSP_FUNC',
                       value=model.__class__.name,
                       comment="Mathematical model of non-linearized data",
                       after='GSP_WREJ')
        ccd.header.set('GSP_ORDR', value=model.degree,
                       comment="Mathematical model order",
                       after='GSP_FUNC')
        ccd.header.set('GSP_NPIX', value=ccd.size,
                       comment="Number of Pixels",
                       after='GSP_ORDR')
        last_keyword = 'GSP_NPIX'
        for i in range(model.degree + 1):
            ccd.header.set('GSP_C{:03d}'.format(i),
                           value=model.__getattribute__('c{:d}'.format(i)).value,
                           comment="Value of parameter c{:d}".format(i),
                           after=last_keyword)

            last_keyword = 'GSP_C{:03d}'.format(i)

        return ccd

    def read_gsp_wcs(self, ccd):
        """Read a GSP-specific wavelength solution

        Args:
            ccd (CCDData) :class:`~astropy.nddata.CCDData` instance

        Returns:
            astropy.modeling.Model instance

        """
        assert isinstance(ccd, CCDData)
        self.model_name = ccd.header['GSP_FUNC']
        self.degree = ccd.header['GSP_ORDR']
        self._binning = int(ccd.header['CCDSUM'].split()[0])

        if self.model_name == 'Chebyshev1D':
            self.model = models.Chebyshev1D(degree=self.degree)
            for i in range(ccd.header['GSP_ORDR'] + 1):
                self.model.__getattribute__('c{:d}'.format(i)).value = ccd.header[
                    'GSP_C{:03d}'.format(i)]
            self.wavelength_and_intensity = [
                self.model(range(ccd.header['GSP_NPIX'])), ccd.data]

            return self.wavelength_and_intensity

    def _model_constructor(self):
        """Generates callable mathematical model

        It can do Chebyshev and Linear model only but is easy to implement
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
        """Wavelength solution fitter

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
                # print(physical, wavelength, self.model)
                self.log.info('Unable to do fit, please add more data points.')
                self.log.error('TypeError: %s', error)
                return None
        else:
            self.log.error('Either model or model fitter were not constructed')
            raise RuntimeError("Undefined model and fitter")

    # wavelength solution reader private methods.
    def _read_non_linear(self, dimension):
        """Non linear solutions reader

        Notes:
            Not all kind of non-linear solutions are implemented. Apparently is
            not fully implemented either.

        Args:
            dimension (int): Solutions can be multi-dimensionals, this method is
                called for each one of them.

        Raises:
            NotImplementedError: This method is not fully implemented.

        """
        # TODO (simon): Complete implementation.
        ctypen = self.wcs.ctype[dimension - 1]
        if ctypen == 'MULTISPE':
            # TODO (simon): What is the * (asterisc) doing here?.
            wat_head = self.header['WAT{:d}*'.format(dimension)]
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
                    wat_string += self.header[key]
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
            function_type = int(spec[11])
            order = int(float(spec[12]))
            min_pix_val = int(float(spec[13]))
            max_pix_val = int(float(spec[14]))
            # resto = list(spec[9:])
            params = [float(i) for i in spec[15:]]
            # TODO (simon): Document self.wcs_dict (what is what).
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
            self.wavelength_and_intensity = [
                self.model(range(self.wcs_dict['pnum'])),
                self.ccd.data]

    def _read_linear(self):
        """Linear solution reader

        This method read the appropriate keywords and defines a linear
        wavelength solution

        Returns:
            Callable wavelength solution model. Instance of
              :class:`~astropy.modeling.Model`
        """
        self.wcs_dict = {'crval': self.wcs.crval[0],
                         'crpix': self.wcs.crpix[0],
                         'cdelt': self.wcs.cd[0],
                         'dtype': self.ccd.header['DC-FLAG'],
                         'pnum': self.ccd.header['NAXIS1']}

        self._set_math_model()

        x_axis = range(len(self.data))

        self.wavelength_and_intensity = [self.model(x_axis), self.data]

        return self.model

    def _set_binning_in_model(self):
        """Modified the model so it matches binned data"""

        self.log.debug('"Binning" model')
        if self.model.__class__.__name__ == 'Linear1D':
            self.model.slope.value *= self._binning
        elif self.model.__class__.__name__ in ['Chebyshev1D',
                                               'Legendre1D']:
            for par_index in range(len(self.model.parameters)):
                parameter_name = self.model.param_names[par_index]
                self.model.__getattribute__(parameter_name).value *= \
                    self.binning ** par_index
        else:
            raise NotImplementedError("Model {:s} not implemented"
                                      "".format(self.model.__class__.__name__))

    def _set_math_model(self):

        if self.wcs_dict['dtype'] == -1:
            self._none()
        elif self.wcs_dict['dtype'] == 0:
            self._linear_solution()
        elif self.wcs_dict['dtype'] == 1:
            self._log_linear()
        elif self.wcs_dict['dtype'] == 2:
            if self.wcs_dict['ftype'] == 1:
                self._chebyshev()
            elif self.wcs_dict['ftype'] == 2:
                self._non_linear_legendre()
            elif self.wcs_dict['ftype'] == 3:
                self._non_linear_cspline()
            elif self.wcs_dict['ftype'] == 4:
                self._non_linear_lspline()
            elif self.wcs_dict['ftype'] == 5:
                # pixel coordinates
                raise NotImplementedError
            elif self.wcs_dict['ftype'] == 6:
                # sampled coordinate array
                raise NotImplementedError
            else:
                raise SyntaxError('ftype {:d} is not defined in the '
                                  'standard'.format(self.wcs_dict['ftype']))
        else:
            raise SyntaxError('dtype {:d} is not defined in the '
                              'standard'.format(self.wcs_dict['dtype']))

    @staticmethod
    def _none():
        """Required to handle No-wavelength solution

        No wavelength solution is considered in the FITS standard (dtype = -1)
        This method is placed here for completness even though is not
        implemented.

        Raises:
            NotImplementedError
        """
        raise NotImplementedError

    def _linear_solution(self):
        """Constructs a linear 1D model based on the WCS information obtained
        from the header.
        """
        intercept = self.wcs_dict['crval'] -\
            (self.wcs_dict['crpix'] - 1) *\
            self.wcs_dict['cdelt']

        self.model = models.Linear1D(slope=self.wcs_dict['cdelt'],
                                     intercept=intercept)

    def _log_linear(self):
        """Not implemented

        Raises:
            NotImplementedError
        """
        # intercept = np.power(10, self.wcs_dict['crval'] +
        #                          self.wcs_dict['cdelt'] *
        #                         (self.wcs_dict['crpix'] - 1))
        #
        # slope = np.power(10, self.wcs_dict['cdelt'])
        #
        # self.model = models.Linear1D(slope=slope,
        #                              intercept=intercept)
        # print(self.model)
        raise NotImplementedError

    def _chebyshev(self):
        """Returns a chebyshev model"""
        self.model = models.Chebyshev1D(degree=self.wcs_dict['order'] - 1,
                                        domain=[self.wcs_dict['pmin'],
                                                self.wcs_dict['pmax']], )
        # self.model.parameters[0] = self.wcs_dict['pmin']
        for param_index in range(self.wcs_dict['order']):
            self.model.parameters[param_index] = self.wcs_dict['fpar'][
                param_index]

    def _non_linear_legendre(self):
        """Set model to Legendre1D
        """
        """Returns a legendre model"""
        self.model = models.Legendre1D(degree=self.wcs_dict['order'] - 1,
                                       domain=[self.wcs_dict['pmin'],
                                               self.wcs_dict['pmax']], )
        # self.model.parameters[0] = self.wcs_dict['pmin']
        for param_index in range(self.wcs_dict['order']):
            self.model.parameters[param_index] = self.wcs_dict['fpar'][
                param_index]
        # raise NotImplementedError

    def _non_linear_lspline(self):
        """Not implemented

        Raises:
            NotImplementedError

        """
        raise NotImplementedError('Linear spline is not implemented')

    def _non_linear_cspline(self):
        """Not implemented

        Raises:
            NotImplementedError

        """
        raise NotImplementedError('Cubic spline is not implemented')

    def get_model(self):
        """Returns the wavelength solution model if exists."""
        if self.model is not None:
            return self.model
        else:
            self.log.error("The solution hasn't been found")


if __name__ == '__main__':
    pass
