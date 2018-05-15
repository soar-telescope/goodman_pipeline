Overview
********
The Goodman Spectroscopic Data Reduction Pipeline - |pipeline name| - is a Python-based
package for producing science-ready, wavelength-calibrated, 1-D spectra. The
goal of |pipeline name| is to provide SOAR users with an easy to use, very well
documented software for reducing spectra obtained with the Goodman spectrograph.
Though the current implementation assumes offline data reduction, our aim is to
provide the capability to run it in real time, so 1-D wavelength calibrated
spectra can be produced shortly after the shutter closes.

The pipeline is primarily intended to be run on a data reduction dedicated
computer. Instructions for running the software are provided in the
`Running Pipeline`_ section of this guide.
The Goodman Spectroscopic Data Reduction Pipeline project is hosted at GitHub at
`it's GitHub Repository <https://github.com/soar-telescope/goodman>`_.


Currently the pipeline is separated into two main components. The initial
processing is done by ``redccd``, which trims the images, and carries out bias
and flat corrections. The spectroscopic processing is done by ``redspec`` and
carries out the following steps:

- Identifies point-source targets.
- Traces the spectra.
- Extracts the spectra.
- Estimates and subtract background.
- Saves extracted (1D) spectra, without wavelength calibration.
- Finds the wavelength solution.
- Linearizes data (resample)
- Writes the wavelength solution to FITS header
- Creates a new file for the wavelength-calibrated 1D spectrum


