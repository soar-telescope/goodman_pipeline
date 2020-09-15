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


Assuming the two targets in `sci_file.fits` are extracted we'll end up with