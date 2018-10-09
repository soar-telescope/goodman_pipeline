Change History
##############

.. _v1.1.3:

V1.1.3 Unreleased
^^^^^^^^^^^^^^^^^

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

- Eliminated ``None`` elements in list of instances of :class:`pipeline.core.core.NightDataContainer`

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
- Relocated module :mod:`pipeline.core.check_version` to ``pipeline/core``.
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