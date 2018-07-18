Change History
##############

12-07-2018 V1.1.0
^^^^^^^^^^^^^^^^^

- Implemented astroscrappy's lacosmic method
- removed ccdproc's lacosmic
- created  'default' method for cosmic ray rejection.
  - For binning 1x1 default is dcr
  - For binning 2x2 default is lacosmic
  - For binning 3x3 default is lacosmic

methods 'dcr', 'lacosmic' or 'none' can still be forced by using --cosmic <method>

11-07-2018 V1.0.3
^^^^^^^^^^^^^^^^^

- Bugs fixed
  - programatically access to the version number did not work because it was
    based purely on `setup.cfg` now `setup.py` has  a function that creates the
    file `pipeline/version.py` which is accessed by `pipeline/__init__.py`
  - File naming was making some file dissapear by being overwritten for files
    that contained more than one target the next file name would match the
    previous one. A differentiator was added.

10-07-2018 V1.0.2
^^^^^^^^^^^^^^^^^

- Removed module `goodman/pipeline/info.py` and placed all metadata in `goodman/setup.cfg`.
- Several updates to documentation
  - Added comment on how to organize data on `soardata3`.
  - Added link to licence on footer.
  - User manual now is in ReadTheDocs and no longer available as a pdf.
  - Improved information on debug plots

- Bugs Fixed.
  - fixed `GSP_FNAM`  value for reference lamps
  - Spectral limit calculation by including binning into the equation
  - Included binning in the calculation of the wavelength solution
  - Corrected messages and conditions under which the prefix for cosmic ray rejection is used
  - Image combination call and messages

- Other additions
  - Added lookup table `dcr.par` file generator and found optimal parameters for Red camera and binning 2x2



# xx-xx-2018 V1.0.1

- Moved user manual from external repo to `goodman/docs/`
- Added version checker
- Centralised metadata (`__version__`, `__licence__`, etc) in `goodman/setup.cfg`
- Added `CHANGES.md`

# 29-04-2018 V1.0.0

- First production ready release