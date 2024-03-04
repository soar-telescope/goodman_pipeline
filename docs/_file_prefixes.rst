.. _file-prefixes:

File Prefixes
*************

.. note::

  Overscan is no longer performed by default. Since it caused a double bias level substraction.
  Fixed since release :ref:`v1.3.3`.


There are several ways one can do this but we selected adding prefixes to the
file name because is easier to add and also easy to filter using a terminal,
for instance.

  ``ls cfzst*fits``

or in python

.. code-block:: python

  import glob

  file_list = glob.glob('cfzst*fits')

So what does all those letter mean? Here is a table to explain it.

.. _table-prefixes:

.. table:: Characters and meaning of prefixes

    ======== ==================================
     Letter   Meaning
    ======== ==================================
     o        Overscan Correction Applied
     t        Trim Correction Applied
     s        Slit trim correction applied
     z        Bias correction applied
     f        Flat correction applied
     c        Cosmic rays removed
     e        Spectrum extracted to 1D
     w        1D Spectrum wavelength calibrated
    ======== ==================================
   

So, for an original file named ``file.fits``:

  ``eczst_file.fits``

Means the spectrum has been extracted to a 1D  file but the file has not been
flat fielded (``f`` missing).

Ideally after running ``redccd`` the file should be named:

  ``cfzst_file.fits``

And after running ``redspec``:

  ``wecfzst_file.fits``
