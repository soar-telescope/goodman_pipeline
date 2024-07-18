Prepare Data for Reduction
**************************

If you did a good job preparing and doing the observation this should be an easy
step, either way, keep in mind the following steps.

- Remove all *focus* sequence.
- Remove all *target acquisition* or *test* frames.
- Using your observation's log remove all unwanted files.
- Make sure all data has the same gain (``GAIN``) and readout noise (``RDNOISE``)
- Make sure all data has the same Region Of Interest or ROI (``ROI``).

The pipeline does not modify the original files unless there are problems with
fits compliance, is never a bad idea to keep copies of your original data in
a safe place.

Updating Keywords
^^^^^^^^^^^^^^^^^
Since version :ref:`1.3.0 <v1.3.0>` if your data is older than **August 6, 2019**, you will
need to change the following keywords.

- ``SLIT``: Replace whitespaces with underscore, remove " and all letters are
  uppercase. For instance ``0.45" long slit`` becomes ``0.45_LONG_SLIT``.
- ``GRATING``: Grating's `lines/mm` goes first and then the manufacturer. For instance.
  ``SYZY_400`` becomes ``400_SYZY``.
- ``WAVMODE``: Replace whitespace with underscore and all letters are capitalized.
  For instance. ``400 m1`` becomes ``400_M1``.
- ``INSTRUME``: Instead of using the classical keywords 'Goodman Spectro' and
  'Goodman Imaging', the AEON standard keywords ``ghts_red`` and ``ghts_blue``
  will be used for spectroscopy, and ``ghts_red_imager`` and ``ghts_blue_imager``
  for imaging. This is an exception of the upper case rule.

.. note::

  General rules are: Underscore is the only accepted separator. All letter must
  be upper case. Remove any character that need escaping.