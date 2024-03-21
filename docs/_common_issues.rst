.. _common_issues:

Common issues
*************

No comparison lamps were found
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The latest version may introduce changes to wavelength solutions. To
view the current list of available modes and lamps, please refer to the
`GitHub repository <https://github.com/soar-telescope/goodman_pipeline/tree/main/goodman_pipeline/data/ref_comp>`_.

If your lamp, filter, or mode is not included in the repository mentioned
above, ``redspec`` will not work as expected. This is particularly common
with custom modes of the |pipeline full name|. You may encounter the following
error after running ``redspec``:

.. code-block:: console

    $ redspec
    ፧
    [15:02:49][W]: No comparison lamps were provided for file ecfzst_0001_science.fits

This error occurs because the |pipeline full name| relies on specific keywords to
extract spectra, such as:

.. _`table lamp key`:

.. table:: Keywords that are used to produce a wavelength solution.

    ========== =============================================================
     Keyword    Purpose
    ========== =============================================================
     LAMP_HGA   Indicates if HgAr lamp is used.
     LAMP_NE    Indicates if Ne lamp is used.
     LAMP_AR    Indicates if Ar lamp is used.
     LAMP_FE    Indicates if Fe lamp is used.
     LAMP_CU    Indicates if Cu lamp is used.
     WAVMODE    Slit and mode configuration.
    ========== =============================================================

Multiple spectra output
^^^^^^^^^^^^^^^^^^^^^^^^

When taking multiple ``ARC`` images during your observation run, they will be linked
with your science data. This means that if you capture several lamp files,
they will be processed alongside the science images, potentially resulting
in multiple outputs of the same spectrum.
