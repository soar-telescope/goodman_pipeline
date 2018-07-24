.. _flat-normalization:

Flat Normalization
******************

There are three possible ``<method>`` (s) to do the normalization of master flats.
For the method using a model the default model's order is ``15``. It can be set
using ``--flat-norm-order <order>``.

``mean``:
  Calculates the mean of the image using numpy's :func:`~numpy.mean` and divide
  the image by it.

``simple`` (default):
  Collapses the master flat across the spatial direction, fits a
  :class:`~astropy.modeling.polynomial.Chebyshev1D` model of order ``15`` and
  divide the full image by this fitted model.

``full``:
  Fits a :class:`~astropy.modeling.polynomial.Chebyshev1D` model to every
  line/column (dispersion axis) and divides it by the fitted model.
  This method takes too much to process and it has been left in the code for
  experimentation purposes only.