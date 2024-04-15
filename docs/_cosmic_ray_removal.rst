.. _cosmic-ray-removal:

Cosmic Ray Removal
******************

.. warning::

  The parameters for either cosmic ray removal method are not fully understood
  neither tuned but they work for most common instrument configurations. If
  your extracted spectrum shows weird features, specially if you use a custom
  mode, the most likely culprit are the parameters of the method you chose.
  Please let us know.

The argument ``--cosmic <method>`` has four options but there are only two real
methods.

``default`` (default):
  Different methods work different for different binning. So if ``<method>`` is
  set to ``default`` the pipeline will decide as follows:

  ``dcr`` for binning ``1x1``

  ``lacosmic`` for binning ``2x2`` and ``3x3`` though binning ``3x3`` has not
  being tested.

``dcr``:
  It was already said that this method work better for binning ``1x1``. More
  information can be found on :ref:`dcr`. The disadvantages of this method is
  that is a program written in C and it is required to write the file to the
  disk, process it and read it back again. Still is faster than ``lacosmic``.

  The parameters for running ``dcr`` are written in a file called ``dcr.par``
  a lookup table and a file generator have been implemented but you can parse
  custom parameters by placing a ``dcr.par`` file in a different directory and
  point it using ``--dcr-par-file <path>``.

``lacosmic``:
  This is the preferred method for files with binning ``2x2`` and ``3x3``.
  This is the Astroscrappy's implementation and is run with the default
  parameters. Future versions might include some parameter adjustment.


``none``:
  Skips the cosmic ray removal process.

Asymetric binnings have not been tested but the pipeline only takes in
consideration the dispersion axis to decide. This does not mean that the spatial
binning does not impact the performance of any of the methods, we just don't
know it yet.

.. note::

  The prefix ``c`` is added to all the comparison lamps, despite they not being
  affected by cosmic rays.