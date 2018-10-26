.. _`New Keywords`:

Description of custom keywords
******************************

The pipeline adds several keywords to keep track of the process and in general
for keeping important information available. The following table gives a description
of all the keywords added by |pipeline name|, though not all of them are
added to all the images.

.. _`general keywords`:

General Purpose Keywords
^^^^^^^^^^^^^^^^^^^^^^^^

These keywords are used for record purpose, except for ``GSP_FNAM`` which is
used to keep track of the file name.

.. _`table general keywords`:

.. table:: General purpose keywords, added to all images at the moment of the first read.

    ========== =============================================================
     Keyword    Purpose
    ========== =============================================================
     GSP_VERS   Pipeline version.
     GSP_ONAM   Original file name, first read.
     GSP_PNAM   Parent file name.
     GSP_FNAM   Current file name.
     GSP_PATH   Path from where the file was read.
     GSP_TECH   Observing technique. Imaging or Spectroscopy.
     GSP_DATE   Date of processing.
     GSP_OVER   Overscan region.
     GSP_TRIM   Trim section.
     GSP_SLIT   Slit trim section. From slit-illuminated area.
     GSP_BIAS   Master bias file used.
     GSP_FLAT   Master flat file used.
     GSP_SCTR   Science target file name (for lamps only)
     GSP_LAMP   Reference lamp used to obtain wavelength solution
     GSP_NORM   Master flat normalization method.
     GSP_COSM   Cosmic ray rejection method.
     GSP_TERR   RMS error of target trace
     GSP_EXTR   Extraction window at first column
     GSP_BKG1   First background extraction zone
     GSP_BKG2   Second background extraction zone
     GSP_WRMS   Wavelength solution RMS Error.
     GSP_WPOI   Number of points used to calculate RMS Error.
     GSP_WREJ   Number of points rejected from RMS Error Calculation.
     GSP_DCRR   Reference paper for DCR software (cosmic ray rejection).
    ========== =============================================================

.. _`target-trace-model`:

Target Trace Model
^^^^^^^^^^^^^^^^^^

.. _`table trace model keywords`:

.. table:: Keywords used to describe the model used to fit the target's trace.

     ========== ====================================================================
      Keyword    Purpose
     ========== ====================================================================
      GSP_TMOD   Name of mathematical model from astropy's :mod:`~astropy.modeling`
      GSP_TORD   Order of the model used.
      GSP_TC00   Value of parameter ``c0``.
      GSP_TC01   Value of parameter ``c1``.
      GSP_TC02   Value of parameter ``c2``. This goes on depending the order.
     ========== ====================================================================


.. _`non-linear wavelength solutions`:

Non-linear wavelength solution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since writing non-linear wavelength solutions to the headers using the FITS
standard (reference) is extremely complex and not necessarily well documented,
we came up with the solution of simply describing the mathematical model
from astropy's :mod:`~astropy.modeling`. This allows for maintaining the data
*untouched* while keeping a reliable description of the wavelength solution.

The current implementation will work for writting any polinomial model. Reading is implemented only for :class:`~astropy.modeling.polynomial.Chebyshev1D` which is the
model by default.

.. _`table non-linear keywords`:

.. table:: Keywords used to describe a non-linear wavelength solution.

     ========== ====================================================================
      Keyword    Purpose
     ========== ====================================================================
      GSP_FUNC   Name of mathematical model from astropy's :mod:`~astropy.modeling`
      GSP_ORDR   Order of the model used.
      GSP_NPIX   Number of pixels.
      GSP_C000   Value of parameter ``c0``.
      GSP_C001   Value of parameter ``c1``.
      GSP_C002   Value of parameter ``c2``. This goes on depending the order.
     ========== ====================================================================


.. _`combined images`:

Combined Images
^^^^^^^^^^^^^^^

Every image used in a combination of images is recorded in the header of the
resulting one. The order does not have importance but most likely the header
of the first one will be used.

The combination is made using the :func:`~ccdproc.combine` method with the following parameters

- ``method='median'``
- ``sigma_clip=True``
- ``sigma_clip_low_thresh=1.0``
- ``sigma_clip_high_thresh=1.0``

At this moment these parameters are not user-configurable.

.. _`table combined images key`:

.. table:: Keywords that list all the images used to produce a combined image.

    ========== =============================================================
     Keyword    Purpose
    ========== =============================================================
     GSP_IC01   First image used to create combined.
     GSP_IC02   Second image used to create combined.
    ========== =============================================================

.. _`detected lines`:
Detected lines
^^^^^^^^^^^^^^

The *reference lamp library* maintains the lamps non-linearized and also they
get a record of the pixel value and its equivalent in angstrom. In the following
table a three-line lamp is shown.

.. _`table line list`:

.. table:: Description of all the keywords used to list lines in lamps in Pixel and Angstrom.

     ========== =============================================================
      Keyword    Purpose                                                     
     ========== =============================================================
      GSP_P001   Pixel value for the first line detected.
      GSP_P002   Pixel value for the second line detected.
      GSP_P003   Pixel value for the third line detected.
      GSP_A001   Angstrom value for the first line detected.
      GSP_A002   Angstrom value for the second line detected.
      GSP_A003   Angstrom value for the third line detected.
     ========== =============================================================
