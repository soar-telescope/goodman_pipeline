.. _file-suffixes:

File Suffixes
*************

After extraction, suffixes may appear in new files created.
There are two scenarios where this can happen:

 - More than one spectroscopic target for extraction. ``*target_X``.
 - More than one comparison lamp. ``*ws_Y``.
 - Both above


 Let's consider the following scenario: We start with 3 reduced files.


.. table:: Sample files in example.

    ================ ========= ========================================
     File Name        Obstype   Comment
    ================ ========= ========================================
     sci_file.fits    OBJECT    Science file with two spectra.
     lamp_001.fits    COMP      Reference lamp valid for sci_file.fits
     lamp_002.fits    COMP      Another valid reference lamp
    ================ ========= ========================================


Assuming the two targets in `sci_file.fits` are extracted and they are approximately at the position
400 and 600 (pixels in spatial axis), after extraction we'll end up with:

.. code-block:: bash

  esci_file_target_1.fits
  esci_file_target_2.fits
  elamp_001_390-410.fits
  elamp_001_590-610.fits
  elamp_002_390-410.fits
  elamp_002_590-610.fits


The default prefix for extraction is ``e`` and does not have an underscore to separate it from the
file name.

After wavelength calibration, since there are two suitable lamps and due to the fact that the
pipeline does not combine solutions, it will save two wavelength calibrated files with each one
solved by the respective lamp. Then:

.. code-block:: bash

  wesci_file_target_1_ws_1.fits
  wesci_file_target_1_ws_2.fits
  wesci_file_target_2_ws_1.fits
  wesci_file_target_2_ws_2.fits
  welamp_001_390-410.fits
  welamp_001_590-610.fits
  welamp_002_390-410.fits
  welamp_002_590-610.fits
