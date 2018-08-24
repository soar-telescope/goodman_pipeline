Change History
##############

.. _v1.1.2:

V1.1.2 Unreleased
^^^^^^^^^^^^^^^^^

- Bugs Fixed:

  + :class:`~pandas.DataFrame` index is unusable when partial parts are eliminated.
    Added ``index_reset(drop=True)``

- Improved some logging messages

  + Informs in debug mode the files that are discarded.

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