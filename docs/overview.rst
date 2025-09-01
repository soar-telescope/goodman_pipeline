Overview
########

The Goodman Data Reduction Pipeline - |pipeline name| - is a
Python-based package for producing science-ready data from both spectroscopic
and imaging observations. The goal of |pipeline name| is to provide SOAR users
with an easy to use, very well documented software for reducing data obtained with the
`Goodman High Throughput Spectrograph <http://www.ctio.noao.edu/soar/content/goodman-high-throughput-spectrograph>`_.
For spectroscopic data, the pipeline produces wavelength-calibrated, 1-D
spectra, while for imaging data it provides calibrated images with astrometric
solutions and photometric measurements.

The pipeline is available for local installation. The |pipeline full name|
project is hosted at GitHub at
`it's GitHub Repository <https://github.com/soar-telescope/goodman_pipeline>`_.
Instructions for running the software are provided in the :ref:`usage` section
of this guide and installation instructions are in :ref:`install`.

Currently the pipeline is separated into several main components. The initial processing is done by ``redccd``, which does the following processes:

- Identifies calibrations and science frames.
- Create master bias.
- Creates master flats and normalizes it.
- Apply overscan correction.
- Trims the image.
- For spectroscopic data find slit edges and trims again.
- Applies bias correction.
- Applies flat correction.
- Applies cosmic ray removal.

The spectroscopic processing is done by ``redspec`` and carries out the following steps:

- Identifies point-source targets.
- Traces the spectra.
- Extracts the spectra.
- Estimates and subtract background.
- Saves extracted (1D) spectra, without wavelength calibration.
- Finds the wavelength solution.
- Linearizes data (resample)
- Writes the wavelength solution to FITS header
- Creates a new file for the wavelength-calibrated 1D spectrum

The astrometric processing is done by ``redastrometry`` and carries out the following steps:

- Detects sources in the image using configurable parameters.
- Creates masks based on flat fields or default circular masks.
- Uses offline Astrometry.net software for plate solving with local index files.
- Determines pixel scale and world coordinate system (WCS) using SIP projections.
- Writes astrometric solution to FITS header.
- Creates a new file with WCS information (_wcs.fits).

The photometric processing is done by ``redphotometry`` and carries out the following steps:

- Uses astrometrically calibrated images as input.
- Performs aperture photometry on detected sources.
- Cross-matches sources with Gaia catalog.
- Applies photometric calibration using Gaia reference stars.
- Provides magnitude measurements and uncertainties.
- Supports configurable aperture sizes and detection parameters.
- Generates diagnostic plots and photometric analysis.

.. include:: _shortcuts.rst
