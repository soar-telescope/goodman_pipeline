.. _reading-wavelength-calibrated-files:

Read Wavelength Calibrated Files
********************************

.. important::

  This change was introduced in :ref:`v3.0.0` in response to several requests to eliminate the
  resampling of wavelength calibrated data. It introduces breaking changes: for instance, the
  wavelength calibrated file is no longer stored with a linear wavelength solution. Instead, it is
  stored as a FITS binary table.

Using astropy.io.fits
^^^^^^^^^^^^^^^^^^^^^

Wavelength calibrated spectra can no longer be read with ``ccdproc.CCDData``. Instead, read them
directly with Astropy's ``fits`` module:

.. code-block:: python

    from astropy.io import fits

    hdul = fits.open("/full/path/to/file.fits")

Calling ``hdul.info()`` shows the file's structure:

.. code-block:: text

    Filename: /full/path/to/file.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU     324   ()
      1  SPECTRUM      1 BinTableHDU     15   2030R x 2C   [D, D]

To get the header, use the HDU name:

.. code-block:: python

    hdul['PRIMARY'].header

or, equivalently, the HDU index:

.. code-block:: python

    hdul[0].header

To access the spectrum, for example to plot it:

.. code-block:: python

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    ax.plot(hdul['SPECTRUM'].data['WAVELENGTH'], hdul['SPECTRUM'].data['INTENSITY'])
    ax.set_xlabel("Wavelength")
    ax.set_ylabel("Intensity")

    plt.show()

Using astropy.table.Table
^^^^^^^^^^^^^^^^^^^^^^^^^

Alternatively, read the ``SPECTRUM`` extension as a table:

.. code-block:: python

    from astropy.table import Table

    tb = Table.read("/full/path/to/file.fits", hdu="SPECTRUM")

    wavelength = tb["WAVELENGTH"]
    intensity = tb["INTENSITY"]

Recreate the mathematical model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mathematical model used to compute the wavelength axis is described in the primary header.
For example, these are the parameters for a *Chebyshev* model of order 3:

.. code-block:: text

    GSP_FUNC= 'Chebyshev1D'        / Mathematical model of non-linearized data
    GSP_ORDR=                    3 / Mathematical model order
    GSP_NPIX=                 2030 / Number of Pixels
    GSP_C000=     2999.82491781623 / Value of parameter c0
    GSP_C001=   1.9959306167200839 / Value of parameter c1
    GSP_C002= 2.10872155121733E-06 / Value of parameter c2
    GSP_C003= -1.2192825594282E-09 / Value of parameter c3

With this information you can recreate the exact model used to build the wavelength axis, using
``astropy.modeling``. Here is a full example:

.. code-block:: python

    from astropy.io import fits
    from astropy.modeling import models

    header = fits.getheader("/full/path/to/file.fits", extname="PRIMARY")

    # GSP_FUNC holds the model name. This example assumes Chebyshev1D and does not use it.
    degree = header['GSP_ORDR']

    model = models.Chebyshev1D(degree=degree)

    for i in range(degree + 1):
        getattr(model, f"c{i}").value = header[f"GSP_C{i:03d}"]

You can then evaluate the model as needed. For example, here is how to plot the spectrum using the
recreated wavelength axis:

.. code-block:: python

    import matplotlib.pyplot as plt
    from astropy.table import Table

    tb = Table.read("/full/path/to/file.fits", hdu="SPECTRUM")
    intensity = tb["INTENSITY"]

    # GSP_NPIX is the number of pixels in the spectrum
    wavelength = model(range(header['GSP_NPIX']))

    fig, ax = plt.subplots()
    ax.plot(wavelength, intensity)
    ax.set_xlabel("Wavelength")
    ax.set_ylabel("Intensity")

    plt.show()

Using the integrated convenience method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``goodman_pipeline`` package includes a convenience method that combines
the steps described above into a single call. See
:meth:`~goodman_pipeline.wcs.wcs.WCS.read_wcs_from_binary_table` for more
details.

.. code-block:: python

    from astropy.io import fits
    from goodman_pipeline.wcs import WCS

    gsp_wcs = WCS()

    full_file_path = "/full/path/to/file.fits"

    hdulist = fits.open(full_file_path)

    wavelength, intensity, model = gsp_wcs.read_wcs_from_binary_table(hdulist=hdulist)

The method returns the wavelength axis (``wavelength``), the spectrum
(``intensity``), and ``model``, the recovered mathematical model of the
wavelength solution. Because the model is stored in the file, you can reuse it
without re-fitting the wavelength axis.

As an example, here is how to plot the spectrum:

.. code-block:: python

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot(wavelength, intensity, label="Spectrum")
    ax.set_xlabel("Wavelength (Angstrom)")
    ax.set_ylabel("Intensity (ADU)")
    ax.legend()

Alternatively, you can recreate the wavelength axis by evaluating the model:

.. code-block:: python

    reconstructed_wavelength = model(range(hdulist['PRIMARY'].header['GSP_NPIX']))

    fig, ax = plt.subplots()
    ax.plot(reconstructed_wavelength, intensity, label="Spectrum with reconstructed wavelength axis")
    ax.set_xlabel("Reconstructed Wavelength (Angstrom)")
    ax.set_ylabel("Intensity (ADU)")
    ax.legend()