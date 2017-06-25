Goodman High Throughput Spectrograph
====================================

**Important** This is a Beta Release Candidate version of this pipeline,
although we have thoroughly tested it there might (and will) still be
bugs. Please let us know if you find any. My e-mail is storres [at] ctio
noao edu.

The Goodman High Throughput Spectrograph is an imaging spectrograph, if
you wish to know more about the instrument please check the `SOAR
website <http://www.ctio.noao.edu/soar/content/goodman-high-throughput-spectrograph>`__

To see full documentation please go to the GitHub hosted site for
`Goodman <https://soar-telescope.github.io/goodman/>`__

What is contained in this package?
----------------------------------

This repository contains tools for processing Goodman's spectrograph
data. It is separated into two main components. The **CCD** part is done
by *redccd* originally developed by `David
Sanmartim <https://github.com/dsanmartim>`__ and currently maintained by
our team. It does the standard *ccd image reduction*, i.e. trim, bias
and flat correction. Currently the ccd reduction pipeline has been
integrated into this package so there is no need to look for it
separately. The **Spectroscopic** processing is done by *redspec* and
includes the following features:

-  [x] Identify targets in images
-  [x] Trace the target
-  [x] Extract the target with background subtraction
-  [x] Find the wavelength Solution, interactive and automatically
-  [x] Linearize data (resample)
-  [x] Write wavelength solution to a FITS header
-  [x] Create a new file for the wavelength calibrated 1D spectrum

There is also a library of calibrated spectrum in FITS format. Different
configurations are provided but is very easy to add your own.

How to install it?
------------------

Specific platform instructions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for installing on Ubuntu 16.04 see `this
wiki <https://github.com/simontorres/goodman/wiki/Ubuntu-16.04-Installation-Experience>`__

for installing on Centos 7 see `this
wiki <https://github.com/simontorres/goodman/wiki/Centos-7-Installation>`__

Below you will find the instrucctions anyways, please refer to the wiki
for the specifics. 

Get the code
~~~~~~~~~~~~

Install Git
^^^^^^^^^^^

This step will depend on your platform, below are two examples:

.. code:: shell

    # Centos
    sudo yum install git

    # Ubuntu
    sudo apt-get install git

Clone from GitHub
^^^^^^^^^^^^^^^^^

To clone from GitHub copy and paste the following line in a terminal.

.. code:: shell

    git clone https://github.com/soar-telescope/goodman.git

You are not ready to install the pipeline yet.

Requirements
^^^^^^^^^^^^

This software was developed on Python 2.7, use the ``requirements.txt``
file to install all dependencies.

.. code:: shell

    sudo -H pip2.7 install -r requirements.txt

Remember that if you do not have super user privileges, you can still
install all the requirements by adding the ``--users`` flag. There is
also an option of `installing it using different virtual
environments <http://docs.python-guide.org/en/latest/dev/virtualenvs/>`__
or even using `Anaconda <https://www.continuum.io/downloads>`__.

You may have some issues when using MatPlotLib on ``virtualenvs`` and on
Mac OSX. If so, you can try to follow the instructions on `this
site <http://matplotlib.org/faq/osx_framework.html#osxframework-faq>`__
and, then add the following line on your ``.bashrc`` or
``.bash_profile`` file.

.. code:: shell

    alias python='pythonw' 

Install DCR (Cosmic Ray Rejection)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This pipeline uses `DCR <http://users.camk.edu.pl/pych/DCR/>`__
developed by `Wojtek Pych <mailto:pych@camk.edu.pl>`__ instead of
``ccdproc.cosmicray_lacosmic`` because we got better results with
``DCR``. Unfortunately you will have to compile it, I have successfully
compiled it on Centos 7, Ubuntu 16.04, Linux Mint 18.1, Solaris 11 and
MacOS Sierra.

Follow `this link <http://users.camk.edu.pl/pych/DCR/>`__ and you can
follow the instructions there. The same instructions are reproduced
here.

Download the ``dcr.tar`` file and untar it.

.. code:: shell

    tar -xvf dcr.tar

Compile it

.. code:: shell

    make

If you don't get any errors you can try it without any arguments and you
will get something like this

.. code:: shell

    $ ./dcr

            USAGE:  dcr  input_file  cleaned_file  cosmicrays_file

    File 'dcr.par' must be present in the working directory.
          ~~~~~~

Make it available for the system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now that you have compiled the program you have a file called ``dcr``
you need to put it in the ``$PATH`` variable of your system. I usually
use ``bash`` so if you use something different follow the example below
as a guide.

1. Create a directory to place the executable

.. code:: shell

    $ mkdir ~/.bin

Note that in the example above the
directory .bin will be hidden and the symbol ``~`` denotes your home
directory for instance: ``/home/goodman/``

2. Move ``dcr`` to your new folder.

.. code:: shell

   $ mv dcr ~/.bin

3. Add the directory to the ``PATH`` variable. With your favorite text
   editor, open the file ``~/.bashrc`` 

.. code:: shell

   $ vim ~/.bashrc

At the end add the following line

.. code:: shell

   export PATH=$PATH:/home/user/.bin

If you don't know your home directory do the following 

.. code:: shell

   $ cd
   $ pwd

Whatever the output is there you should replace it for
   ``/home/user/``

4. Reload the environment variables. For this you can simply close and
   reopen the terminal or you can do:
   
.. code:: shell

    $ source ~/.bashrc


How to use it?
--------------

The pipeline is separated in two sub-pipelines. *redccd* and *redspec*.
The ``--help`` argument will print the argument plus some some
description

.. code:: shell 

   $ redccd --help usage: redccd [-h] [-c] [--ignore-bias] [--auto-clean]
                                 [--saturation ] [--raw-path raw_path]
                                 [--red-path red_path] [--debug]
                                 [--log-to-file] [--flat-normalize ]
                                 [--flat-norm-order ] [--dcr-par-dir ]
                                 [--keep-cosmic-files]

    Goodman CCD Reduction - CCDreductions for Goodman spectroscopic data
    
    optional arguments: 
      -h, --help              show this help message and exit -c,
      --cosmic                Clean cosmic rays from science data. 
      --ignore-bias           Ignore bias correction 
      --auto-clean            Automatically clean reduced data directory
      --saturation            Saturation limit. Default to 55.000 ADU (counts)
      --raw-path raw_path     Path to raw data.
      --red-path red_path     Path to reduced data.
      --debug                 Show detailed information of the process.
      --log-to-file           Write log to a file.
      --flat-normalize        Choose a method to normalize the master
                              flat forspectroscoy. Choices are: mean, simple
                              (model) and full (fits model to each line).
      --flat-norm-order       Defines the order of the model to be fitted.
      --dcr-par-dir           Directory of default dcr.par file.
      --keep-cosmic-files     After cleaning cosmic rays with dcr, do not remove
                              the input file and the cosmic rays file.


And for ``redspec``:

.. code:: shell

   $ redspec --help usage: redspec [-h] [--data-path ]
                                   [--proc-path ]
                                   [--search-pattern ]
                                   [--output-prefix ] [--extraction ]
                                   [--reference-files ] [--interactive]
                                   [--debug][--log-to-file]
                                   [--save-plots] [--plot-results]

    Extracts goodman spectra and does wavelength calibration.
    
    optional arguments:
      -h, --help              show this help message and exit
      --data-path             Path for location of raw data. Default <./> 
      --proc-path             Path for destination of processed data. Default <./> 
      --search-pattern        Pattern for matching the goodman's reduced data. 
      --output-prefix         Prefix to add to calibrated spectrum. 
      --extraction            Choose a which extraction to perform. Simple is a
                              sum across the spatial direction after the
                              background has been removed. Optimal is a more
                              advanced method that considers weights and
                              profile fitting. 
      --reference-files       Directory of Reference files location 
      --interactive           Interactive wavelength solution.Disbled by
                              default.
      --debug                 Debugging Mode 
      --log-to-file           Write log to a file
      --save-plots            Save all plots in a directory
      --plot-results          Show wavelength calibrated spectrum at the end.


You should always run ``redccd`` first and then ``redspec``. There are
certain defaults values

redccd Defaults
~~~~~~~~~~~~~~~

.. code:: shell

        --cosmic              False
        --ignore-bias         False
        --auto-clean          False
        --debug               False
        --log-to-file         False
        --keep-cosmic-files   False
        --saturation          55000
        --raw-path            ./
        --red-path            ./RED/
        --flat-normalize      simple
        --dcr-par-dir         files/
        --flat-norm-order     15

redspec Defaults
~~~~~~~~~~~~~~~~

.. code:: shell

        --data-path         ./
        --proc-path         ./
        --search-pattern    cfzsto
        --extraction        simple
        --reference-files   refdata
        --reference-lamp    (empty string)
        --output-prefix     g
        --interactive       False
        --debug             False
        --log-to-file       False
        --save-plots        False
        --plot-results      False

