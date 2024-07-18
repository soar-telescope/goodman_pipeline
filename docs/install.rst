.. _install:

Install
#######

We do not have the resources to provide installation support, thus we provide
a server with the latest and older versions installed that users with access rights
can use. However, the installation process is simple, and is described below (we hope that in enough detail)
so that users can follow the steps and end up with a running version of the pipeline.

.. note::

  In this tutorial we use Miniconda3 as but you can also do it with Anaconda by
  visiting `<https://www.anaconda.com/download/>`_ for downloading
  the Anaconda installer for either operating system.

Ubuntu
******

This installation process has been tested in a live version
(previously known as *Live CD*) of Ubuntu 18.04.

1. Install ``make`` and ``gcc``

  ``sudo apt install gcc make``

2. Download Miniconda

  ``wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh``

Mac OS
******

This installation was tested on MacOS High Sierra.

1. Install `Xcode <https://itunes.apple.com/us/app/xcode/id497799835?mt=12>`_ command Line Tools. You can do it from the *App Store* or from command line.

  ``xcode-select --install``

2. Download Anaconda or Miniconda

  ``curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh --output miniconda.sh``

If you have a different processor like intel, you can search miniconda installer url for your system.

Common Steps
************

.. warning::

  This installation instructions work with **BASH** only.

1. Install Miniconda

  ``bash miniconda.sh``

  Answer the questions and reopen the terminal at the end.

2. Configure conda to use the Astroconda Channel

  ``conda config --add channels http://ssb.stsci.edu/astroconda``

3. Get the `latest release <https://github.com/soar-telescope/goodman_pipeline/releases>`_
   of the |pipeline full name| from github. There is a ``.zip`` and ``.tar.gz``
   files, the following steps will continue with the latter.

   Make sure there is a green tag that says Latest release. For instance, for
   the version ``v1.1.2`` the file is called.

   ``goodman_pipeline-1.2.1.tar.gz``


4. Extract the file (you can't copy and paste this line)

  ``tar -xvf goodman_pipeline-<latest tag>.tar.gz``

  where ``<latest tag>`` is the version number of the latest release.

5. Move into the directory containing the package.

  ``cd goodman_pipeline-<latest tag>``

  If you do ``ls`` you should find some interesting files on it such as:
  ``setup.py`` and ``environment.yml`` and ``install_dcr.sh``.

6. Create the virtual environment. The ``environment.yml`` file contains a preset
   definition of a virtual environment that Anaconda will understand, also
   ensures that the |pipeline full name| will work. Even the name of the virtual
   environment is set there.

  ``conda env create -f environment.yml``

  This will create a virtual environment called ``goodman_pipeline``. To activate it

  ``source activate goodman_pipeline``

7. Install ``DCR``. This script requires a virtual environment activated.

  ``sh install_dcr.sh``

  To test if it worked you can do:


  ``dcr``

  You should get something like this::

      (goodman_pipeline) [user@servername goodman_pipeline]$ dcr

            This is a modified version of DCR! for the Goodman Spectroscopic Pipeline
            Please visit the author's site to get the original version:
            Modification Version: 0.0.1

            http://users.camk.edu.pl/pych/DCR/


            USAGE:  dcr  input_file  cleaned_file  cosmicrays_file

      File 'dcr.par' must be present in the working directory.
          ~~~~~~

8. Run tests.

  ``python setup.py test``

9. Install the pipeline

  ``python setup.py install``

Using pip
*********

The |pipeline full name| Can also be installed using pip, but it does not install
``dcr``, so if you are updating your goodman pipeline only you can use:

  ``pip install goodman-pipeline``


.. include:: _shortcuts.rst
