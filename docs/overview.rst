Overview
########
The Goodman Spectroscopic Data Reduction Pipeline - |pipeline name| - is a
Python-based package for producing science-ready, wavelength-calibrated, 1-D
spectra. The goal of |pipeline name| is to provide SOAR users with an easy to
use, very well documented software for reducing spectra obtained with the
`Goodman High Troughput Spectrograph <http://www.ctio.noao.edu/soar/content/goodman-high-throughput-spectrograph>`_.
Though the current implementation assumes offline data reduction, our aim is to
provide the capability to run it in real time, so 1-D wavelength calibrated
spectra can be produced shortly after the shutter closes.

The pipeline is primarily intended to be run on a data reduction dedicated
computer though it is available for local installation. The |pipeline full name|
project is hosted at GitHub at
`it's GitHub Repository <https://github.com/soar-telescope/goodman_pipeline>`_.

Instructions for running the software are provided in the :ref:`usage` section
of this guide. How to access the the data reduction server is on
:ref:`remote-access` or if you prefer to install a local version instructions
are in :ref:`install`




Currently the pipeline is separated into two main components. The initial
processing is done by ``redccd``, which does the following processess.

- Identifies calibrations and science frames.
- Create master bias.
- Creates master flats and normalizes it.
- Apply overscan correction.
- Trims the image.
- For spectroscopic data find slit edges and trims again.
- Applies bias correction.
- Applies flat correction.
- Applies cosmic ray removal.


The spectroscopic processing is done by ``redspec`` and carries out the
following steps:

- Identifies point-source targets.
- Traces the spectra.
- Extracts the spectra.
- Estimates and subtract background.
- Saves extracted (1D) spectra, without wavelength calibration.
- Finds the wavelength solution.
- Linearizes data (resample)
- Writes the wavelength solution to FITS header
- Creates a new file for the wavelength-calibrated 1D spectrum


.. include:: _shortcuts.rst