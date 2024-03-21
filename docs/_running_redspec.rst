Extracting the spectra
**********************

After you are done :ref:`processing-2d-images` it is time to extract the
spectrum into a wavelength-calibrated 1D file.

The script is called ``redspec``. The tasks performed are the following:

- Classifies data and creates the match of ``OBJECT`` and ``COMP`` if it exists.
- Identifies targets
- Extracts targets
- Saves extracted targets to 1D spectrum
- Finds wavelength solution automatically
- Linearizes data
- Saves wavelength calibrated file

First you have to move into the ``RED`` directory, this is a precautionary method
to avoid unintended deletion of your raw data. Then you can simply do:

  ``redspec``

And the pipeline should work its magic, though this might not be the desired
behavior for every user or science case, we have implemented a set of
*command line arguments* which are listed below.

- ``--data-path <path>`` Folder were data to be processed is located. Default
  is *current working directory*.
- ``--proc-path <path>`` Folder were processed data will be stored. Default
  is *current working directory*.
- ``--search-pattern <pattern>`` Prefix for picking up files. Default is
  ``cfzst-``. See :ref:`file-prefixes`.
- ``--output-prefix <prefix>`` Prefix to be added to calibrated spectrum. Default is
  ``w-``. See :ref:`file-prefixes`.
- ``--extraction <method>`` Select the :ref:`extraction-methods`. The only one
  implemented at the moment is ``fractional`` .
- ``--fit-targets-with {moffat, gaussian}`` Model to fit peaks on spatial profile
  while searching for spectroscopic targets. Default ``moffat``.
- ``--target-min-width <target min width>`` Minimum profile width for fitting the spatial axis of spectroscopic targets.
  If fitting a Moffat it will be reflected as the FWHM attribute of the fitted model and if fitting a Gaussian it will
  be reflected as the STDDEV attribute of the Gaussian model.
- ``--target-max-width <target max width>`` Maximum profile width for fitting the spatial axis of spectroscopic targets.
  If fitting a Moffat it will be reflected as the FWHM attribute of the fitted model and if fitting a Gaussian it will
  be reflected as the STDDEV attribute of the Gaussian model.
- ``--reference-files <path>`` Folder where to find the reference lamps.
- ``--debug`` Shows extended and more messages.
- ``--debug-plot`` Shows plots of intermediate steps.
- ``--max-targets <value>`` Maximum number of targets to detect in a single
  image. Default is 3.
- ``--background-threshold <background threshold>`` Multiplier for background level used to discriminate usable targets.
  Default 3 times background level.
- ``--save-plots`` Save plots.
- ``--plot-results`` Show plots during execution.

The mathematical model used to define the wavelength solution is recorded
in the header even though the data has been linearized for record purpose.
