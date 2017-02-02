Goodman High Throughput Spectrograph
====================================

**WARNING** This is the first release of this pipeline, although we have
thoroughly tested it there might still be bugs. Please let me know by an
e-mail to storres [at] ctio noao edu.

The Goodman High Throughput Spectrograph is an imaging spectrograph, if
you wish to know more about the instrument please check the `SOAR
website <http://www.ctio.noao.edu/soar/content/goodman-high-throughput-spectrograph>`__

To see the GitHub site go `here <https://github.com/simontorres/goodman>`__

What is contained in this package?
----------------------------------

This repository contains tools for spectroscopy, but after the data is
reduced, i.e. bias and flat corrected, for that part please use `David
Sanmartim's <https://github.com/dsanmartim/goodman_ccdreduction>`__
github repository.

This package have the following capabilities.

- Identify targets in images
- Trace the target
- Extract the target with background subtraction
- Find the wavelength Solution (Requires User Input)
- Linearize data
- Write wavelength solution to a FITS header
- Create a new file for the wavelength calibrated 1D spectrum

How to install it?
------------------

Either clone or download this code. If you decide to clone it, just do:

.. code:: shell

    git clone https://github.com/simontorres/goodman.git

Or you can simply go and `click
here <https://github.com/simontorres/goodman/archive/master.zip>`__ for
download a zip file with all the script files.

Requirements
~~~~~~~~~~~~

This software was developed on Python 2.7, use the ``requirements.txt``
file to install all dependencies.

.. code:: shell

    sudo -H pip2.7 install -r requirements.txt

How to use it?
--------------

To get a list of the possible arguments do:

.. code:: shell

    /path/to/this/repo/redspec.py --help

The simplest way of running this pipeline is to go to your data folder,
already processed with
`goodman\_ccdreduction <https://github.com/dsanmartim/goodman_ccdreduction>`__
and execute ``redspec.py``

.. code:: shell

    /path/to/this/repo/redspec.py

Will run the following defaults: 

- Observing Mode **0**: One solution applied to all night 
- Interactive Mode **True**
- Data Path **./** 
- Destination folder for processed data **./**
- Search Pattern **fzh\_** 
- Output prefix **g** 
- Reference Files Path [/path/to/this/repo/]**refdata/** 
- Plots Enabled **False** 
- Reference Lamp **None**
