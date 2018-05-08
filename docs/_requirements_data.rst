.. _`data requirements`:
Data Requirements
*****************
The |goodman HTS|'s data has seen some evolution in the past years, in shape and
most importantly in its headers. The |pipeline full name| relies heavily on the data's
header so this is in fact very important.

The headers must be `FITS Compliant <https://fits.gsfc.nasa.gov/fits_standard.html>`_,
otherwise the software exits with errors.

Bear in mind that the |goodman hts| has `two cameras <http://www.ctio.noao.edu/soar/content/goodman-spectrograph-overview>`_, *Blue* and *Red*.

Recent data obtained with the Red Camera already meets all the requirements in
terms of format and header cards. Data obtained with the Blue Camera before
|headers change| is expected to have several format issues:

    1. There were non fits-compliant characters in some comments. To solve that, you can edit the header using the most recent version of AstroPy, IRAF or WCSTOOLS to remove the following keywords: ``PARAM0``, ``PARAM61``, ``PARAM62`` and ``PARAM63``.

    2. The data was defined as 3D, just like a single frame of a data cube. To solve this, you will have to read the data and rewrite it with only two dimensions using AstroPy or IRAF.

    3. Some keywords were added with time.

        * INSTCONF: contains the name of the Goodman Camera used, e.g., "Blue" or "Red".

        * WAVMODE: contains the ruling number of the used grating and the mode name, e.g., "400 m1" or "400 m2".

        * ROI: the name of the region of interest, e.g., "Spectroscopic 1x1", "user-defined", etc.

    4. Duplicated keywords. Make sure that your data does not contain duplicated keywords.

Reference Lamp Files
^^^^^^^^^^^^^^^^^^^^
Having an *automatic wavelength calibration method* relies on having previously calibrated
reference lamps obtained in the same configuration or mode. It is also important
that the lamp names are correct, for instance ``HgAr`` is quite different than
``HgArNe``. Spaces between lamps are not allowed. And the name is case
insensitive: you may write "HgAr", "HGAR" or "hgar".
The list of current lamps is the following.


.. _`Table Supported Modes`:

.. table:: List of |goodman hts| supported modes

   ========= ====== ======== ======================================
    Grating   Mode   Filter    Lamp   
   ========= ====== ======== ======================================
      400      M1    None     HgAr, HgArNe
      400      M2    GG455    Ar, Ne, HgAr, HgArNe, CuHeAr, FeHeAr
   ========= ====== ======== ======================================


.. important::

    More lamps will be made public shortly.

.. raw:: pdf

    PageBreak

.. _`Header Requirements`:

Headers Requirements
~~~~~~~~~~~~~~~~~~~~

Goodman HTS spectra have small non-linearities on their wavelength solutions.
They are small but must be taken into account.

It was necessary to  implement a custom way of storing non-linear wavelength
solutions that at the same time allowed for keeping data as *untouched* as
possible. The main reason is that linearizing the reference lamps made it
harder to track down those non-linearities on the new data being calibrated. Also,
The documentation on how to write non-linear solution to a FITS header is
far too complex for our use case and there is no apparent effort on trying to
simplify it. Below we compile a list of required keywords for
comparison lamps if they want to be used as reference lamps. The full list of
keywords is listed under `New Keywords`_.

General Custom Keywords:
  Every image processed with the *Goodman Spectroscopic Pipeline* will have the
  `general keywords`_. The one required for a reference lamp is the following:

    ``GSP_FNAM = file-name.fits / Current file name``


Record of `detected lines`_ in Pixel and Angstrom:
  Reference lamps have a record of usable lines in its header. Initially the lamp
  is run through a tool that identifies the lines and records its pixel value.
  The root string is ``GSP_P`` followed by a zero-padded three digit sequential number
  (001, 002, etc). For instance.

    ``GSP_P001= 499.5377036976768  / Line location in pixel value``

    ``GSP_P002= 810.5548319623747  / Line location in pixel value``

    ``GSP_P003= 831.6984711087946  / Line location in pixel value``

  Later, the equivalent values in angstrom are then recorded with the root string
  ``GSP_A`` and the same numerical pattern as before.

    ``GSP_A001= 5460.75            / Line location in angstrom value``

    ``GSP_A002= 5769.61            / Line location in angstrom value``

    ``GSP_A003= 5790.67            / Line location in angstrom value``


  ``GSP_P001`` and ``GSP_A001`` are a match. If any of the angstrom value entries
  have a value of ``0`` (default value) the equivalent pair pixel/angstrom entry is ignored.
  Also they must be organized in an *always increasing* way, if they are not, they
  will be ignored too.

.. important::

  Those keywords are used to calculate the mathematical fit of the
  wavelength solution and are not used on normal operation. Our philosophy here
  is that the line identification has to be done only once and then the
  model can be fitted several times, actually you can try several models
  if you want. (On your own)

.. raw:: pdf

    PageBreak

`Non-linear wavelength solution`_:
  The method for recording the non-linear wavelength solution is actually
  very simple. It requires: ``GSP_FUNC`` which stores a string with the name of
  the mathematical model from ``astropy.modeling.models``. ``GSP_ORDR`` stores
  the order or degree of the model. ``GSP_NPIX`` stores the number of pixels in
  the spectral axis. Then there is N+1 parameter keywords where N is the order
  of the model defined by ``GSP_ORDR``. The root string of the keyword is ``GSP_C``
  and the rest is a zero-padded three digit number starting on zero to N.
  See the example below.

    ``GSP_FUNC= Chebyshev1D          / Mathematical model of non-linearized data``

    ``GSP_ORDR= 3                    / Mathematical model order``

    ``GSP_NPIX= 4060                 / Number of Pixels``

    ``GSP_C000= 4963.910057577853    / Value of parameter c0``

    ``GSP_C001= 0.9943952599223119   / Value of parameter c1``

    ``GSP_C002= 5.59241584012648e-08 / Value of parameter c2``

    ``GSP_C003= -1.2283411678846e-10 / Value of parameter c3``

.. warning::

    This method has been developed and tested to write correctly polynomial-like
    models. And ONLY reads ``Chebyshev1D`` models.
    Other models will just be ignored. More development will be done based on
    request, suggestions or needs.

File organization
^^^^^^^^^^^^^^^^^
redccd and redspec will look for all FITS files inside the current working
directory or inside the path provided with the ``--raw-path`` (redccd)/``--data-path`` (redspec)
flag non-recursively. Make sure to have only data that contains relevant signal.
Data obtained during the focusing process, saturated flats, etc, must be removed.

Also, we recommend you follow these good practices:

- Delete all unnecessary files (focus,  test, acquisition, unwanted exposures, etc)
- Don't mix different ROI (Region Of Interest), Gain and Readout Noises.
- Make sure all the required file types are present: BIAS, FLAT, COMP, OBJECT.


