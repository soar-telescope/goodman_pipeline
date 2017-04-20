Goodman High Throughput Spectrograph
====================================

**WARNING** This is the first release of this pipeline, although we have
thoroughly tested it there might still be bugs. Please let me know by an
e-mail to storres [at] ctio noao edu.

The Goodman High Throughput Spectrograph is an imaging spectrograph, if
you wish to know more about the instrument please check the `SOAR
website <http://www.ctio.noao.edu/soar/content/goodman-high-throughput-spectrograph>`__

To see full documentation please go to the GitHub hosted site for
`Goodman <https://simontorres.github.io/goodman/>`__

What is contained in this package?
----------------------------------

This repository contains tools for processing Goodman's spectrograph
data, from the ccd image on. It is separated into two main components.
The **CCD** part is done by *redccd* originally developed by `David
Sanmartim's <https://github.com/dsanmartim/goodman_ccdreduction>`__ and
currently maintained by our team. It does the standard *ccd image
reduction*, i.e. trim, bias and flat correction. Currently the ccd
reduction pipeline has been integrated into this package so there is no
need to look for it separately. The **Spectroscopic** processing is done
by *redspec* and includes the following features:

-  [x] Identify targets in images
-  [x] Trace the target
-  [x] Extract the target with background subtraction
-  [x] Find the wavelength Solution (Requires User Input)
-  [x] Linearize data
-  [x] Write wavelength solution to a FITS header
-  [x] Create a new file for the wavelength calibrated 1D spectrum

There is also a library of calibrated spectrum in FITS format. Different
configurations are provided but is very easy to add your own.

How to install it?
------------------

Specific platform instructions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*for installing on Ubuntu 16.04 see `this
wiki <https://github.com/simontorres/goodman/wiki/Ubuntu-16.04-Installation-Experience>`__*

*for installing on Centos 7 see `this
wiki <https://github.com/simontorres/goodman/wiki/Centos-7-Installation>`__*

Get the code
~~~~~~~~~~~~

Clone from GitHub
^^^^^^^^^^^^^^^^^

Either clone or download this code. If you decide to clone it, just do:

.. code:: shell

    git clone https://github.com/simontorres/goodman.git

Download package-only tarball
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Go to `this
site <https://github.com/simontorres/goodman/tree/master/dist>`__ and
download the latest (or highest) version and place it wherever you want

Requirements
~~~~~~~~~~~~

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

.. code:: bash

    alias python='pythonw' 

How to use it?
--------------

To get a list of the possible arguments do:

.. code:: shell

    $ redspec --help

The simplest way of running this pipeline is to go to your data folder,
already processed with ``redccd`` (originally
`goodman\_ccdreduction <https://github.com/dsanmartim/goodman_ccdreduction>`__)
and execute ``redspec``

.. code:: shell

    $ redspec

Will run the following defaults: - [ ] Observing Mode **0**: One
solution per configuration applied to all night - [ ] Interactive Mode
**True** - [ ] Data Path **./** - [ ] Destination folder for processed
data **./** - [ ] Search Pattern **fzh\_** - [ ] Output prefix **g** - [
] Reference Files Path [/path/to/this/repo/]**refdata/** - [ ] Plots
Enabled **False** - [ ] Reference Lamp **None**
