Goodman High Throughput Spectrograph
====================================

The Goodman High Throughput Spectrograph (Goodman HTS) Data-Reduction Pipeline
is the SOAR Telescope's official data reduction pipeline for *Goodman HTS*.

It has been fully developed in Python 2.7 and mostly astropy affiliated packages
with the exception of `dcr <http://users.camk.edu.pl/pych/DCR/>`_ an external tool
that does cosmic ray identification and correction. The reason for using it
instead of LACosmic was that it works very well for spectroscopic data and the
results are evidently superior. Some of the negative aspects of using this
external (meaning outside of Python domains) software were: The integration into
the pipeline's workflow and the use of an external ``dcr.par`` parameter file.
Such parameters have to be changed by hand and can't be integrated into the
pipeline's workflow itself. In particular for binning 2x2 and custom ROI those
parameters contained in *dcr.par* has to be specifically tuned.


The pipeline is currently in its first beta release version. We have thoroughly
tested it but we know there are still bugs to be found. If you found one or have
any comment please contact Simón Torres at storres [at] ctio noao edu.

The Goodman High Throughput Spectrograph is an imaging spectrograph, if
you wish to know more about the instrument please check the `SOAR
website <http://www.ctio.noao.edu/soar/content/goodman-high-throughput-spectrograph>`__

To see full documentation please go to the GitHub hosted site for
`Goodman <https://soar-telescope.github.io/goodman/>`__

What is contained in this package?
----------------------------------

The Goodman HTS Pipeline has two packages, *goodman_ccd* dedicated to the basic
ccd reduction process and *goodman_spec* for advanced processing for
spectroscopic data, such as target identification and extraction, for instance.

They are supposed to work in sequence and they are called from two indepentent
scripts, *redccd* and *redspec*.

Other files included are:
- Readme.md (this)
- requirements.txt
- user_manual.pdf
- dcr binaries
- Library of Wavelength Calibrated lamps.

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
In case you want to get the full project you need to install ``git``, if you
don't, please go to `Get latest release <#get-latest-release>`__.

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

Get Latest Release
^^^^^^^^^^^^^^^^^^

Visit `this link <https://github.com/soar-telescope/goodman/tree/master/dist>`__
and download the newest file in the list.

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

Replace ``/home/user/`` for whatever the output is in the last command


4. Reload the environment variables. For this you can simply close and
   reopen the terminal or you can do:
   
.. code:: shell

    $ source ~/.bashrc


How to use it?
--------------

The pipeline is run by two separated scripts. *redccd* and *redspec*.
The ``--help`` argument will print the argument plus some short description.
Also check the manual included in the distribution package.

Acknowledge
-----------

`David Sanmartim <https://github.com/dsanmartim>`__ developed the first version
of the ccd reduction part before he moved to a new position. That first version
was included in the pipeline as a package and has evolved along all the code.