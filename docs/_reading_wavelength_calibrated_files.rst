.. _reading-wavelength-calibrated-files:

Read Wavelength Calibrated Files
********************************

.. important::

  This change was introduced in :ref:`v3.0.0` as a consequence of several requirements received
  requesting to eliminate resampling of the wavelength calibrated data. It introduces breaking changes
  for instance the wavelength calibrated file is no longer stored as a linear wavelength solution,
  instead it creates a FITS Binary Table.

Using astropy.io.fits
^^^^^^^^^^^^^^^^^^^^^
There are several ways of reading the spectrum, it is no longer possible to read it using ``ccdproc.CCDData`` anymore.

Now we read directly using Astropy's fits.

.. code-block:: python

    from astropy.io import fits

    hdul = fits.open("/full/path/to/file.fits")

And if we do ``hdul.info()`` we get the file's structure.

.. code-block:: shell

    Filename: /full/path/to/file.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU     324   ()
      1  SPECTRUM      1 BinTableHDU     15   2030R x 2C   [D, D]

So in order to get the header we can do:

.. code-block:: python

    hdul['PRIMARY'].header


Or in this case its equivalent:

.. code-block:: python

    hdul[0].header


And to get access the spectrum, for a plot for instance:

.. code-block:: python

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    ax.plot(hdul['SPECTRUM'].data['WAVELENGTH'], hdul['SPECTRUM'].data['INTENSITY'])

    plt.show()


Using astropy.table.Table
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from astropy.table import Table

    tb = Table.read("/full/path/to/file.fits", hdu="SPECTRUM")

    wavelength = tb["WAVELENGTH"]
    intensity = tb["INTENSITY"]



Recreate Mathematical Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mathematical model is described in the primary header. Here is an example of the the parameters for
a *Chebyshev* of order 3:

.. code-block:: text

    GSP_FUNC= 'Chebyshev1D'        / Mathematical model of non-linearized data
    GSP_ORDR=                    3 / Mathematical model order
    GSP_NPIX=                 2030 / Number of Pixels
    GSP_C000=     2999.82491781623 / Value of parameter c0
    GSP_C001=   1.9959306167200839 / Value of parameter c1
    GSP_C002= 2.10872155121733E-06 / Value of parameter c2
    GSP_C003= -1.2192825594282E-09 / Value of parameter c3

With this information it is possible to recreate the exact mathematical model used to create the wavelength axis
and it is based on ``astropy.modeling.Model``

So, here is a full example:

.. code-block:: python

    from astropy.io import fits
    from astropy.modeling import models

    header = fits.getheader("/full/path/to/file.fits", extname="PRIMARY")

    model_name = header['GSP_FUNC'] # This is where the model name is defined, in this example is not used
    degree = header['GSP_ORDR]

    model = models.Chebyshev1D(degree=degree)

    for i in range(degree + 1):
        model.__getattribute__(f"c{i:d}").value = header[f"GSP_C{i:03d}"]

Then you can use model as you wish, for instance you could plot and compare, for simplicity here we
present the example of plotting:

.. code-block:: python

    from astropy.table import Table

    tb = Table.read("/full/path/to/file.fits", hdu="SPECTRUM")

    # wavelength = tb["WAVELENGTH"]
    intensity = tb["INTENSITY"]

    # use GSP_NPIX to get the number of pixels in the spectrum
    wavelength = model(range(header['GSP_NPIX']))

    # reutilize the intensity obtained before

    fig, ax = plt.subplots()
    ax.plot(wavelength, intensity)

    plt.show()
