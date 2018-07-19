.. _flat-normalization:

Flat Normalization
******************

There are three possible ``<method>``(s) to do the normalization of master flats.
For the method using a model the default model's order is 15. It can be set using
``--flat-norm-order <order>``.

``mean``:
  Calculates the mean of the image using ``numpy.mean`` and divide the image by it.

``simple``:
  Collapses the master flat across the spatial direction, fits a model and divide
  the full image by this fitted model.

``full``:
  Fits a model to every line/column (dispersion axis) and divide by itself by the model.
  This method takes too much to process and it has been left in the code for
  experimentation purposes.