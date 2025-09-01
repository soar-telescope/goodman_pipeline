Performing Photometry
*********************

The ``redphotometry`` command performs aperture photometry on astrometrically calibrated images, providing magnitude measurements cross-matched with Gaia catalog sources for photometric calibration.

Prerequisites
^^^^^^^^^^^^^

Before using ``redphotometry``, you must have:

- An astrometrically calibrated image (``_wcs.fits`` file from ``redastrometry``)
- Optionally, a master flat field matching the filter and size of your image
- Internet connection for Gaia catalog queries (if not cached locally)

Basic Usage
^^^^^^^^^^^

The typical workflow starts with an astrometrically calibrated image:

.. code-block:: bash

   redphotometry image_wcs.fits

This will perform photometry on ``image_wcs.fits`` and produce photometric results.

Recommended Usage with Flat Field
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning::

   It is **highly recommended** to use a master flat image to obtain a precise mask of the useful image area:

.. code-block:: bash

   redphotometry image_wcs.fits --flat flat_V.fits

The master flat must match the input file's filter, otherwise a default circular mask (90% of shorter axis diameter) will be used.

Key Parameters
^^^^^^^^^^^^^^

Photometry Configuration
""""""""""""""""""""""""

- ``--aperture-radius RADIUS`` - Aperture radius in pixels (default: 4.0)
- ``--initial-fwhm FWHM`` - Initial FWHM in pixels for source detection (default: 3.0)
- ``--detection-threshold THRESHOLD`` - Detection threshold above noise factor (default: 5.0)

Gaia Integration
""""""""""""""""

- ``--gaia-sources-limit LIMIT`` - Maximum number of Gaia sources to process (default: 5000)
- ``--gaia-photometry-column COLUMN`` - Gaia filter corresponding column name

Processing Options
""""""""""""""""""

- ``--flat FLAT`` or ``--flat-image FLAT`` - Path to master flat image
- ``--disable-mask-creation`` - Disable mask creation (default: False)
- ``--overwrite`` - Allow overwriting existing photometry results
- ``--aperture-curve-of-growth`` - Run aperture curve of growth diagnostics and exit

FITS Header Keywords
""""""""""""""""""""

- ``--imaging-filter-keyword KEYWORD`` - FITS header keyword for imaging filter (default: 'FILTER')

Output and Debugging
""""""""""""""""""""

- ``--plots`` - Show diagnostic plots
- ``--debug`` - Enable debug mode

Example Commands
^^^^^^^^^^^^^^^^

Basic photometry with plots:

.. code-block:: bash

   redphotometry image_wcs.fits --plots

With custom aperture and Gaia settings:

.. code-block:: bash

   redphotometry image_wcs.fits --aperture-radius 5.0 --gaia-sources-limit 3000 --plots

Using flat field and overwriting existing results:

.. code-block:: bash

   redphotometry image_wcs.fits --flat flat_V.fits --overwrite --plots

Aperture curve of growth analysis:

.. code-block:: bash

   redphotometry image_wcs.fits --aperture-curve-of-growth

Complete Workflow Example
^^^^^^^^^^^^^^^^^^^^^^^^^^

Here's a typical processing workflow from raw image to photometry:

.. code-block:: bash

   # Step 1: Perform astrometric calibration
   redastrometry raw_image.fits --flat flat_V.fits --plots

   # Step 2: Perform photometry on the astrometrically calibrated image
   redphotometry raw_image_wcs.fits --flat flat_V.fits --aperture-radius 4.5 --plots

   # Optional: Run aperture curve of growth analysis
   redphotometry raw_image_wcs.fits --aperture-curve-of-growth

Output Files and Results
^^^^^^^^^^^^^^^^^^^^^^^^

``redphotometry`` produces several outputs:

- Photometric catalog with source positions and magnitudes
- Cross-matched Gaia sources for calibration
- Magnitude measurements with uncertainties
- Diagnostic plots (if ``--plots`` is enabled)

The photometric results include:

- Instrumental magnitudes for detected sources
- Calibrated magnitudes using Gaia reference stars
- Photometric uncertainties and quality flags
- Source positions in pixel and world coordinates

Aperture Photometry Details
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The photometry process involves:

1. **Source Detection** - Using configurable FWHM and threshold parameters
2. **Aperture Photometry** - Measuring flux within specified radius
3. **Background Estimation** - Local background subtraction around each source
4. **Gaia Cross-matching** - Matching detected sources with Gaia catalog
5. **Photometric Calibration** - Using Gaia sources as photometric standards
6. **Quality Assessment** - Providing measurement uncertainties and flags

Troubleshooting
^^^^^^^^^^^^^^^

- Ensure input file is astrometrically calibrated (has WCS solution)
- If few sources detected, try lowering ``--detection-threshold``
- For crowded fields, consider increasing ``--detection-threshold``
- Adjust ``--aperture-radius`` based on seeing conditions and source sizes
- Use ``--debug`` mode for detailed processing information
- Verify flat field matches the filter of the input image
- Check internet connection for Gaia catalog access

Filter Considerations
^^^^^^^^^^^^^^^^^^^^^

The photometric calibration depends on matching the observing filter with appropriate Gaia photometric bands. Ensure that:

- The FITS header contains correct filter information
- The ``--imaging-filter-keyword`` points to the correct header keyword
- The ``--gaia-photometry-column`` matches the appropriate Gaia band for your filter
