Obtaining Astrometric Solution
******************************

The ``redastrometry`` command performs astrometric calibration using offline Astrometry.net software to determine world coordinate system (WCS) solutions with SIP projections.

Prerequisites
^^^^^^^^^^^^^

Before using ``redastrometry``, ensure you have:

- Downloaded the appropriate astrometry.net index files (see :ref:`install` section)
- A science image to process
- Optionally, a master flat field matching the filter and size of your image

Basic Usage
^^^^^^^^^^^

The simplest way to run astrometric calibration:

.. code-block:: bash

   redastrometry image.fits

This will process ``image.fits`` and create ``image_wcs.fits`` with the astrometric solution.

Recommended Usage with Flat Field
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning::

   It is **highly recommended** to use a master flat image to obtain a precise mask of the useful image area:

.. code-block:: bash

   redastrometry image.fits --flat flat_V.fits

The master flat must match the input file's filter, otherwise a default circular mask will be used.

Key Parameters
^^^^^^^^^^^^^^

Source Detection
"""""""""""""""""

- ``--initial-fwhm FWHM`` - Initial FWHM in pixels for source detection (default: 3.0)
- ``--detection-threshold THRESHOLD`` - Detection threshold above noise factor (default: 5.0)

Astrometry Configuration
""""""""""""""""""""""""

- ``--pixel-scale SCALE`` - Expected pixel scale in arcsec/pixel (default: 0.15)
- ``--pixel-scale-tolerance TOL`` - Tolerance for pixel scale matching (default: 0.02)
- ``--downsample-factor FACTOR`` - Downsample factor for astrometry.net (default: 2)

Processing Options
""""""""""""""""""

- ``--flat FLAT`` or ``--flat-image FLAT`` - Path to master flat image
- ``--index-directory DIR`` - Custom directory for astrometry.net index files
- ``--ignore-goodman-vignetting`` - Skip pre-source detection masking for non-Goodman images
- ``--overwrite`` - Allow overwriting existing astrometry results

Output and Debugging
"""""""""""""""""""""

- ``--plots`` - Show diagnostic plots
- ``--debug`` - Enable debug mode
- ``--verbose`` - Enable verbose astrometry.net logging

Example Commands
^^^^^^^^^^^^^^^^

Basic processing with plots:

.. code-block:: bash

   redastrometry image.fits --plots

With custom parameters and flat field:

.. code-block:: bash

   redastrometry image.fits --flat flat_V.fits --pixel-scale 0.12 --detection-threshold 4.0 --plots

Using custom index directory:

.. code-block:: bash

   redastrometry image.fits --index-directory /path/to/astrometry/index --overwrite

Output Files
^^^^^^^^^^^^

``redastrometry`` produces several output files:

- ``image_wcs.fits`` - Main output with astrometric solution
- Additional diagnostic files and plots (if ``--plots`` is enabled)

The ``_wcs.fits`` file contains the original image data with updated FITS headers including:

- World Coordinate System (WCS) keywords
- SIP distortion coefficients
- Astrometric quality metrics

Next Steps
^^^^^^^^^^

The ``image_wcs.fits`` file can be used as input for photometric processing with ``redphotometry``:

.. code-block:: bash

   redphotometry image_wcs.fits

Troubleshooting
^^^^^^^^^^^^^^^

- If astrometry fails, try adjusting ``--pixel-scale`` and ``--pixel-scale-tolerance``
- For crowded fields, consider increasing ``--detection-threshold``
- Use ``--debug`` mode to get detailed processing information
- Ensure you have the appropriate index files for your image's field of view
- Verify that your flat field matches the filter of the input image
