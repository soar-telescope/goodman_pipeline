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

The 2D images are initially reduced using ``redccd``. You can simply move to the
directory where your raw data is located and do:

  ``redccd``

Though you can modify the behavior in several ways.

Running ``redccd`` will create a directory called ``RED`` where it will put your
reduced data. If you want to run it again it will prevent you from accidentally
removing your already reduced data unless you use ``--auto-clean`` this will
tell the pipeline to delete the ``RED`` directory.

  ``redccd --auto-clean``

A summary of the most important arguments are presented below.

- ``--cosmic`` Let you select the method to do :ref:`cosmic-ray-removal`.
- ``--debug`` Show extended messages and plots of intermediate steps.
- ``--flat-normalize`` Let you select the method to do :ref:`flat-normalization`.
- ``--flat-norm-order`` Set order for the model used to do
  :ref:`flat-normalization`. Default 15.
- ``--ignore-bias`` Ignores the existence or lack of ``BIAS`` data.
- ``--ignore-flats`` Ignores the existence or lack of ``FLAT`` data.
- ``--raw-path`` Set the directory where the raw data is located, can be relative.
- ``--red-path`` Set the directory where the reduced data will be stored. Default ``RED``.
- ``--saturation`` Set the saturation level. Flats exceeding the saturation
  level will be discarded. Default 65.000 ADU.








