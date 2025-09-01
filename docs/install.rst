.. _install:

Install
#######

**The preferred method of installation is using pip**: ``pip install goodman-pipeline``

We do not have the resources to provide installation support. However, the installation process is simple, and is described below (we hope that in enough detail) so that users can follow the steps and end up with a running version of the pipeline.

.. note::

  In this tutorial we use Miniconda3 but you can also do it with Anaconda by
  visiting `<https://www.anaconda.com/download/>`_ for downloading
  the Anaconda installer for either operating system.

.. warning::

  **For most cases doing ``pip install goodman-pipeline`` is enough.** You might need to compile DCR separately (see instructions below).

Using pip (Recommended)
***********************

The |pipeline full name| can be installed using pip, which is the preferred and simplest installation method:

**For the latest stable version:**

  ``pip install goodman-pipeline``

**For the latest release candidate (if available):**

  ``pip install --pre goodman-pipeline``

**Note:** This method does not install ``dcr``, so if you need cosmic ray rejection, you'll need to install it separately (see the DCR installation section below).

If you are updating your goodman pipeline installation, you can use:

  ``pip install --upgrade goodman-pipeline``

Installing DCR
**************

If you need cosmic ray rejection functionality, you'll need to install DCR separately:

1. Download the DCR installation script from the repository
2. Run the installation script in your active environment:

  ``sh install_dcr.sh``

To test if DCR installed correctly:

  ``dcr``

You should see output similar to::

    This is a modified version of DCR! for the Goodman Spectroscopic Pipeline
    Please visit the author's site to get the original version:
    Modification Version: 0.0.1

    http://users.camk.edu.pl/pych/DCR/

    USAGE:  dcr  input_file  cleaned_file  cosmicrays_file

File 'dcr.par' must be present in the working directory.

Installing Astrometry.net Index Files
*************************************

If you plan to use the ``redastrometry`` command for astrometric solutions, you'll need to download astrometry.net index files locally. The choice of index files depends on your image's pixel scale, which varies with CCD binning.

Goodman Pixel Scale by Binning
==============================

- **1x1 binning (unbinned)**: 0.15 arcsec/pixel
- **2x2 binning**: 0.30 arcsec/pixel
- **3x3 binning**: 0.45 arcsec/pixel

Selecting Index Files
=====================

Choose index files based on your typical binning mode:

+----------+----------------+---------------------------+------------------------+
| Binning  | Pixel Scale    | Recommended Index Series  | Field of View Coverage |
+==========+================+===========================+========================+
| 1x1      | 0.15"/pixel    | 4200-4204                 | ~0.1° to 2°            |
+----------+----------------+---------------------------+------------------------+
| 2x2      | 0.30"/pixel    | 4200-4206                 | ~0.2° to 4°            |
+----------+----------------+---------------------------+------------------------+
| 3x3      | 0.45"/pixel    | 4200-4207                 | ~0.3° to 6°            |
+----------+----------------+---------------------------+------------------------+

Download Script
===============

Use this script to download the necessary index files::

    #!/bin/bash
    # Target directory for index files
    INDEX_DIR="$HOME/astrometry/index"
    mkdir -p "$INDEX_DIR"

    # Range of index series (4200 to 4207)
    START_INDEX=4200
    END_INDEX=4204  # Adjust based on your binning needs

    # Range of tile numbers (00 to 47)
    TILE_MIN=0
    TILE_MAX=47

    # Base URL
    BASE_URL="http://data.astrometry.net/4200"

    echo "Downloading index files for series ${START_INDEX}–${END_INDEX}..."

    for SERIES in $(seq $START_INDEX $END_INDEX); do
      for TILE in $(seq -w $TILE_MIN $TILE_MAX); do
        FILENAME="index-${SERIES}-${TILE}.fits"
        DEST="${INDEX_DIR}/${FILENAME}"
        URL="${BASE_URL}/${FILENAME}"

        if [ -f "$DEST" ]; then
          echo "✓ $FILENAME already exists, skipping."
        else
          echo "↓ Downloading $FILENAME..."
          curl -L -o "$DEST" "$URL"
          if [ $? -eq 0 ]; then
            echo "✓ Downloaded $FILENAME"
          else
            echo "✗ Failed to download $FILENAME"
            rm -f "$DEST"  # Clean up incomplete file
          fi
        fi
      done
    done

    echo "Download complete! Index files saved to: $INDEX_DIR"

Usage Notes
===========

- **For 1x1 binning**: Use ``END_INDEX=4204`` in the script
- **For 2x2 binning**: Use ``END_INDEX=4206`` in the script
- **For 3x3 binning**: Use ``END_INDEX=4207`` in the script
- **Mixed binning**: Use ``END_INDEX=4207`` to cover all cases

Using Custom Index Directory
============================

If you download index files to a custom location, specify the path when running redastrometry::

    redastrometry image.fits --index-directory $HOME/astrometry/index

Alternative Installation Methods
********************************

Ubuntu
======

This installation process has been tested in a live version
(previously known as *Live CD*) of Ubuntu 18.04.

1. Install ``make`` and ``gcc``

  ``sudo apt install gcc make``

2. Download Miniconda

  ``wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh``

Mac OS
======

This installation was tested on MacOS High Sierra and newer versions.

1. Install `Xcode <https://itunes.apple.com/us/app/xcode/id497799835?mt=12>`_ command Line Tools. You can do it from the *App Store* or from command line.

  ``xcode-select --install``

2. Download Anaconda or Miniconda

  ``curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh --output miniconda.sh``

If you have a different processor like intel, you can search miniconda installer url for your system.

Manual Installation from Source
*******************************

.. warning::

  This installation method is for advanced users only. The pip installation method is recommended for most users.

.. note::

  This installation instructions work with **BASH** only.

1. Install Miniconda (if not already installed)

  ``bash miniconda.sh``

  Answer the questions and reopen the terminal at the end.

2. Configure conda to use the Astroconda Channel

  ``conda config --add channels http://ssb.stsci.edu/astroconda``

3. Get the `latest release <https://github.com/soar-telescope/goodman_pipeline/releases>`_
   of the |pipeline full name| from github. There is a ``.zip`` and ``.tar.gz``
   files, the following steps will continue with the latter.

   Make sure there is a green tag that says Latest release. For instance, for
   the version ``v1.3.11`` the file is called:

   ``goodman_pipeline-1.3.11.tar.gz``


4. Extract the file (you can't copy and paste this line)

  ``tar -xvf goodman_pipeline-<latest tag>.tar.gz``

  where ``<latest tag>`` is the version number of the latest release.

5. Move into the directory containing the package.

  ``cd goodman_pipeline-<latest tag>``

  If you do ``ls`` you should find files such as ``pyproject.toml`` and ``install_dcr.sh``.

6. Create the virtual environment. The ``environment.yml`` file contains a preset
   definition of a virtual environment that Anaconda will understand, and
   ensures that the |pipeline full name| will work. Even the name of the virtual
   environment is set there.

  ``conda env create -f environment.yml``

  This will create a virtual environment called ``goodman_pipeline``. To activate it:

  ``source activate goodman_pipeline``

7. Install ``DCR`` (if needed). This script requires a virtual environment activated.

  ``sh install_dcr.sh``

  To test if it worked you can do:

  ``dcr``

  You should get the DCR usage information as shown above.

8. Install the pipeline using pip

  ``pip install .``

  Or for an editable installation (if you plan to modify the code):

  ``pip install -e .``

Verification
************

To verify that the installation was successful, you can run:

``redccd --help``

and

``redspec --help``

Both commands should display their respective help messages without errors.

Repository Information
**********************

The Goodman Pipeline repository is located at: https://github.com/soar-telescope/goodman_pipeline

For issues, bug reports, or feature requests, please visit the `GitHub Issues page <https://github.com/soar-telescope/goodman_pipeline/issues>`_.

.. include:: _shortcuts.rst
