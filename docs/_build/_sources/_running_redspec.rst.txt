Run redspec
***********

Is the spectroscopy reduction script. The tasks are the following:

- Classifies data and creates the match of ``OBJECT`` and ``COMP`` if it exists.
- Identifies targets
- Extracts targets
- Saves extracted targets to 1D spectrum
- Finds wavelength solution automatically
- Linearizes data
- Saves wavelength calibrated file

The mathematical model used to define the wavelength solution is recorded
in the header even though the data has been linearized for record purpose.
