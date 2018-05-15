.. _dev-install-iraf:

Install Iraf on Centos 7
########################

Get the Files
*************

Get the files from the `iraf website <http://iraf.noao.edu/>`_ and get the
following files:

    ``iraf.lnux.x86_64.tar.gz``

and

    ``x11iraf-v2.0BETA-bin.redhat.tar.gz``

Install Iraf
************

Installing iraf is quite simple.

    ``sudo mkdir -p /iraf/iraf``

    ``cd /iraf/iraf/``

    ``sudo tar -xvf /full/path/to/iraf.lnux.x86_64.tar.gz``

    ``export iraf=/iraf/iraf/``

    ``sudo $iraf/install``

Install ``x11iraf`` Binaries
****************************

Go to the location of  ``x11iraf-v2.0BETA-bin.redhat.tar.gz`` and untar it

    ``tar -xvf x11iraf-v2.0BETA-bin.redhat.tar.gz``

    ``sudo ./install``

Since ``x11iraf is a 32bits distribution and most systems now are 64bits install
the following dependencies:

    ``sudo yum install glibc.i686``

    ``sudo yum install libXmu-1.1.2-2.el7.i686``


    ``sudo yum install ncurses-libs-5.9-13.20130511.el7.i686``

If you can't find the right package do one fo these commands:

    ``sudo yum provides '*libXmu.so.6'``

or

    ``sudo yum provides '*libncurses.so.5'``

And choose the 32bits version.

Final Adjustments
*****************

Go to your home directory and run:

    ``mkiraf``