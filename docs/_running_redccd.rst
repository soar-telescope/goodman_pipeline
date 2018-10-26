.. _processing-2d-images:

Processing your 2D images
*************************

It is the first step in
the reduction process, the main tasks are listed below.

- Create master bias
- Create master flats
- Apply Corrections:

  + Overscan
  + Trim image
  + Detect slit and trim out non-illuminated areas
  + Bias correction
  + Normalized flat field correction
  + Cosmic ray rejection

.. note::

  Some older Goodman HTS data has headers that are not FITS compliant,
  In such cases the headers are fixed and that is the only modification done to
  raw data.

The 2D images are initially reduced using ``redccd``. You can simply move to the
directory where your raw data is located and do:

  ``redccd``

Though you can modify the behavior in several ways.

Running ``redccd`` will create a directory called ``RED`` where it will put your
reduced data. If you want to run it again it will prevent you from accidentally
removing your already reduced data unless you use ``--auto-clean`` this will
tell the pipeline to delete the ``RED`` directory and start over.

  ``redccd --auto-clean``

A summary of the most important *command line arguments* are presented below.

- ``--cosmic <method>`` Let you select the method to do :ref:`cosmic-ray-removal`.
- ``--debug`` Show extended messages and plots of intermediate steps.
- ``--flat-normalize <method>`` Let you select the method to do :ref:`flat-normalization`.
- ``--flat-norm-order <order>`` Set order for the model used to do
  :ref:`flat-normalization`. Default 15.
- ``--ignore-bias`` Ignores the existence or lack of ``BIAS`` data.
- ``--ignore-flats`` Ignores the existence or lack of ``FLAT`` data.
- ``--raw-path <path>`` Set the directory where the raw data is located, can be relative.
- ``--red-path <path>`` Set the directory where the reduced data will be stored. Default ``RED``.
- ``--saturation <saturation>`` Set the saturation threshold in percentage. There
  is a table with all the readout modes and their values at which saturation is
  reached, then all the pixels exceeding that value are counted. If the percentage
  is larger that the threshold defined with this argument the flat is marked as
  saturated. The default value is 1 percent.


This is intended to work with *spectroscopic* and *imaging* data, that it is why
the process is split in two.






