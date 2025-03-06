Change History
##############

.. _v1.3.8:

V1.3.8 Not Released
^^^^^^^^^^^^^^^^^^^

- Documentation update
- Fix DCR bug that was present when called from CLI.
- Fixed Bug when using relative path for selecting master flat
- Fixed Bug when dcr didn't processed properly when comparison lamps have OBSTYPE = ARC


.. _v1.3.7:

V1.3.7 17-08-2023
^^^^^^^^^^^^^^^^^

- Updated documentation theme to allow dark mode
- Removed python 3.6 and 3.8 from testing and added 3.9 and 3.10
- Unfixed numpy version
- Removed testing with conda on github action because it only added complexity
- Added 3 seconds timeout to request to github to validate version and check if running version is latest
- Added correct exception handler in case the http request timed out.


.. _v1.3.6:

V1.3.6 25-03-2022
^^^^^^^^^^^^^^^^^

- Updated HgArNe reference lamps for 400M1 and 600MID GG395
- Fixed deprecation warning  for numpy.float conversion

.. _v1.3.5:

V1.3.5 24-11-2021
^^^^^^^^^^^^^^^^^

- Created `--target-min-width` and `--target-max-width` flags to allow extraction of extended sources.
- Replaced coveralls with codecov for code coverage.
- Documentation improvements and updates.
- Fixed issue with extracted spectrum that where masked in 2D.
- Code style improvements.


.. _v1.3.4:

V1.3.4 11-05-2021
^^^^^^^^^^^^^^^^^

- Updated default search pattern for `redspec`
- Fixed branch name in Github Actions


.. _v1.3.3:

V1.3.3 19-04-2021
^^^^^^^^^^^^^^^^^

- Fixed bug that caused overscan and bias correction happen in the same image
  resulting in double bias level subtraction.

.. _v1.3.2:

V1.3.2 09-10-2020
^^^^^^^^^^^^^^^^^

- Fixed Github Actions setup
- Removed pandas version constraint and implemented workaround to be able to
  use latest pandas version.
- Modified installation of dcr on travis and Github Actions

.. _v1.3.1:

V1.3.1 23-09-2020
^^^^^^^^^^^^^^^^^

- Documentation improvements:

  + More docstrings.
  + New :ref:`file-suffixes` section.

- Added more test code.
- Cleaned .travis.yml and created special dcr installation script for travis.
- Changed the way `core.combine_data` names new files.
- Fixed version checker due to deprecation of access token as url parameter.
- New `core.identify_technique`. Was developed in the web application context.
- Created `--skip-slit-trim` argument to provide more control for certain use cases.
- Removed python 3.5 because it will not be supported anymore.
- Improved AEON Support:

  + Added values for OBSTYPE and required logic.
  +

- Bugs Fixed:

  + Serial and Parallel binning extraction from header was not working
  + Changed url for astroplan's server

.. _v1.3.0:

V1.3.0 06-03-2020
^^^^^^^^^^^^^^^^^

- Made it compatible with Astropy 4.0
- All versions are free except for Pandas [#314]
- `wavelength.WavelengthCalibration.__call__` can now return a json output.
- `core.setup_logging` can now create a generic logger (same format).
- Modified how master bias are named.
- Removed bias overscan and trimming correction on master bias creation.
- Bugs Fixed:

  + `--max-targets` was not being used, missed connection in `ReduceCCD`.

- Updated keyword values of all reference lamps in the library according to [#292]
- Refactored `wavelength.WavelengthCalibration` class moving several methods to
  `core` [#300, #303]
- Refactored `wavelength.WavelengthCalibration` to be instantiated without
  arguments.
- Improved messages at critical stages of wavelength calibration.
- Moved `setup_logging` call from main package `__init__` to scripts or entry
  points, this allows to re use other master loggers.
- Changed `--background-threshold` to multiply by detection limit instead of
  background level
- Created standard JSON output for :class:`~wavelength.WavelengthCalibration`.


.. _v1.2.1:

V1.2.1 19-08-2019
^^^^^^^^^^^^^^^^^

- Bugs fixed

  + Bias process was not fully ignored when `--ignore-bias` was used [#289].
  + `pandas` version was not specified in `environment.yml` [#288, #290]
  + Target extraction failed for low signal targets because background subtraction
    was being ignored at the step of actually identifying targets.
- Install instructions updated [#290]
- Moved static methods from `ImageProcessor` to `core`.
- Added function to validate ccd regions using regular expressions.
- Using lamps keywords to select reference lamps.
- Replaced `target_stddev` by `target_fwhm` in function `extract` and `extract_fractional`.
- Replaced `nsigmas` by `nfwhm` everywhere.
- Added argument `--background-threshold` with default value `3`.
- Added argument `--fit-targets-with` with options `moffat` and `gaussian`.


.. _v1.2.0:

V1.2.0 26-10-2018
^^^^^^^^^^^^^^^^^
- Bugs removed:

  + If there was more than one lamp for a science target the lamp recorded as used
    was one of them only.
  + A percentage symbols was added to the help of ``--saturation`` argument, this
    caused a crash when ``redccd -h`` or ``redccd --help`` was used.
- Numpy is fixed to the version ``1.15.2`` until further notice.
- Reference lamps now get the extraction window added to the end of the file name.
  This is to avoid overwritting the lamps when they were used for more than one target.
- DCR install script is now more advanced and requires a virtual environment to work on.
- Added SOAR Logo to ReadTheDocs page.
- Changed install instruction with exact steps and commands instead of
  referencing documentation.
- Improved method to detect saturated images. Added a table with the *half full
  well* for all the readout modes possible and created a method to easily
  retrieve the value. This is a big improvement since in earlier versions the
  saturation limit was set to 65000 ADU regardless the input data and the user
  had to set a different one using the argument ``--saturation``.
- Repurposed the ``--saturation`` command line argument, now is used to define
  the percentage of pixels above the saturation level, which for simplicity is
  the value of half full well. A default value of 1 percent was set as default.
- Added record information of target trace into the header and logs.
- Added record of background extraction regions into the header and logs.
- Made all plots full screen and the images using the ``gray`` cmap.
- Trace information is printed in the logs and also is recorded in the image's
  header
- Added sigma clipping to target tracing functions

.. _v1.1.2:

V1.1.2 05-10-2018
^^^^^^^^^^^^^^^^^

- Version 1.1.2 is pip instalable

  ``pip install goodman-pipeline``

- Project and package renamed to ``goodman_pipeline`` this is because the
  previous was too generic. Now we have this structure::

   goodman_pipeline/
      docs/
      goodman_pipeline/
         core/
         images/
         ..etc
      setup.py
      ..etc

- Bugs Fixed:

  + :class:`~pandas.DataFrame` index is unusable when partial parts are eliminated.
    Added ``index_reset(drop=True)``
  + Data conversion from string to integer needed to be converted to float first.

  + For low SNR data there was confusion of noise with targets, added a median
    filter and increased the the ``order`` value of peak detection.

- Created several new keywords:

  ``GSP_EXTR``:
    Extraction window at the first column.

  ``GSP_SCTR``:
    Used for extracted comparison lamps, contains the name of the file of
    science target that the lamp was extracted for.

  ``GSP_LAMP``:
    For science targets, it records the name of the lamp used for the wavelength
    calibration.

- "Sliding" cross correlation window (to trace non-linearity of wavelength
  solution) is set to the maximum value between the length of the lamp spectrum
  in pixels and four times the global cross correlation of the reference lamp to
  the new one.

- Iterations in sigma clipping of differences between obtained wavelength
  values and laboratory values was increased from 1 to 3. This is for removing
  bad fitted lines and also RMS error calculation.

- Gaussian Kernel size for reference lamp convolution is now dependent on slit size and binning

- Added reference lamps for all gratings and their modes except ``1200M0``

- Created script ``install_dcr.sh``

- Increased code coverage

- Eliminated ``None`` elements in list of instances of :class:`goodman_pipeline.core.core.NightDataContainer`

- Improved several logging messages

  + In general, it informs more, when it does an action and when it does not.
    What files are discarded,
  + Debugging plots are more complete for ``identify_targets``.

- Created new argument ``--debug-plot`` dedicated for *graphical debugging*, the
  old ``--debug`` will show additional messages but will not produce any
  graphical output.

- Removed ability to process several folders in sequence, now the pipeline has to
  be run for each folder separately.

.. _v1.1.1:

V1.1.1 23-08-2018
^^^^^^^^^^^^^^^^^

- Bugs Fixed:

  + Added clean exit when pipeline is unable to determine ``instrument`` or
    ``technique`` used.
  + Conversion from string to integer not always works, added intermediate float
    conversion.
  + Abrupt exit when there were non-fits-compliant keywords. Now it attempts to
    fix them all automatically and warns the user. Also, it ends the execution
    and informs the user to try again.

- Removed unused code and tools.
- Relocated module :mod:`goodman_pipeline.core.check_version` to ``pipeline/core``.
- Implemented Authorized GitHub API access and added actual version check
- Moved *command line interface* from ``goodman/bin/`` to ``goodman/pipeline/script/``
- Specified version of :mod:`cython` to be able to build.
- Added reference lamps for all usable modes for the grating 600 l/mm
- Created method to use automatic keyword fix from :mod:`~ccdproc`.
- Improved help information of arguments
- Documentation updates

.. _v1.1.0:

V1.1.0 24-07-2018
^^^^^^^^^^^^^^^^^
- Bugs fixed

  + ``--keep-cosmic-file`` would work for ``dcr`` but not for ``lacosmic``

- Changed organization of ReadTheDocs information

  + New structure
  + Added references to external packages
  + This page is the single place to add changes information. CHANGES.md still
    exist but contains a link here.

- Added ``--version`` argument.
- Implemented `astroscrappy's` LACosmic method
- removed ccdproc's :func:`~ccdproc.cosmicray_lacosmic`.
- created  ``default`` method for cosmic ray rejection.

  + For binning 1x1 default is dcr
  + For binning 2x2 default is lacosmic
  + For binning 3x3 default is lacosmic

methods ``dcr``, ``lacosmic`` or ``none`` can still be forced by using
``--cosmic <method>``

.. _v1.0.3:

V1.0.3 11-07-2018
^^^^^^^^^^^^^^^^^

- Bugs fixed

  + programatically access to the version number did not work because it was
    based purely on ``setup.cfg`` now ``setup.py`` has  a function that creates the
    file :mod:`pipeline.version` which is accessed by ``pipeline/__init__.py``
  + File naming was making some file dissapear by being overwritten for files
    that contained more than one target the next file name would match the
    previous one. A differentiator was added.

.. _v1.0.2:

V1.0.2 10-07-2018
^^^^^^^^^^^^^^^^^

- Removed module ``goodman/pipeline/info.py`` and placed all metadata in ``goodman/setup.cfg``.
- Several updates to documentation

  + Added comment on how to organize data on ``soardata3``.
  + Added link to licence on footer.
  + User manual now is in ReadTheDocs and no longer available as a pdf.
  + Improved information on debug plots

- Bugs Fixed.

  + fixed ``GSP_FNAM``  value for reference lamps
  + Spectral limit calculation by including binning into the equation
  + Included binning in the calculation of the wavelength solution
  + Corrected messages and conditions under which the prefix for cosmic ray rejection is used
  + Image combination call and messages

- Other additions
  + Added lookup table ``dcr.par`` file generator and found optimal parameters for Red camera and binning 2x2

.. _v1.0.1:

V1.0.1 xx-xx-2018
^^^^^^^^^^^^^^^^^

- Moved user manual from external repo to ``goodman/docs/``
- Added version checker
- Centralised metadata (``__version__``, ``__licence__``, etc) in ``goodman/setup.cfg``
- Added ``CHANGES.md``

.. _v1.0.0:

V1.0.0 29-04-2018
^^^^^^^^^^^^^^^^^

- First production ready release
