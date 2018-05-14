Software Requirements
*********************
Using the pipeline remotely is the recommended method, in which case you don't need
to worry about software requirements.

However, for users who wish to go ahead with a local installation, we provide
simple instructions in the `Setup for local installation`_ section.

Setup for Remote Use
^^^^^^^^^^^^^^^^^^^^
The Goodman Spectroscopic Data Reduction Pipeline has been installed on a
dedicated computer at SOAR. The procedure requires to open a VNC session, for which
you need to be connected to the SOAR VPN. The credentials for the VPN are the
same you used for your observing run, provided by your *Support Scientist*, who
will also give you the information for the data reduction computer VNC
connection.

.. note:: IRAF is available in the data server at SOAR. Running ``iraf`` will
    open an *xgterm* and *ds9* windows. ``iraf-only`` will open *xgterm* but
    not *ds9*

Establish a VNC connection
~~~~~~~~~~~~~~~~~~~~~~~~~~
Separately, you should receive a server hostname, IP, display number and
VNC-password.

.. _`VNC Displays table`:
.. table:: VNC display number and working folder assigned to each partner.

   ========= ===================== ====================================
    Display    Partner/Institution     Folder
   ========= ===================== ====================================
       :1      NOAO                  ``/home/goodman/data/NOAO``
       :2      Brazil                ``/home/goodman/data/BRAZIL``
       :3      UNC                   ``/home/goodman/data/UNC``
       :4      MSU                   ``/home/goodman/data/MSU``
       :5      Chile                 ``/home/goodman/data/CHILE``
   ========= ===================== ====================================

For this tutorial we will call the vnc server host name as ``<vnc-server>``
the display number  is ``<display-number>`` and your password is ``<password>``.

The VNC connection should work with any VNC Client like TightVNC, TigerVNC,
RealVNC, etc. The first two run on Linux and can be used directly with the
``vncviewer`` command line.

VNC from the Terminal
~~~~~~~~~~~~~~~~~~~~~
Find the ``<display-number>`` that corresponds to you from the `VNC Displays table`_.
Open a terminal, and assuming you have installed ``vncviewer``.

    ``vncviewer <vnc-server>:<display-number>``

You will be asked to type in the ``<password>`` provided.

.. important::

    The real values for ``<vnc-server>`` and ``<password>``
    should be provided by your support scientist.

If the connection succeeds you will see a *Centos 7* Desktop using *Gnome*.

Setup for local installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The |pipeline name| is completely written in Python 3.x and relies on several
libraries like:

* NumPy
* SciPy
* MatPlotLib
* Pandas
* AstroPy
* AstroPy/ccdproc
* AstroPy/astroplan
* DCR

We **do not** recommend the installation of these libraries or the
|pipeline name| in your system since updates and upgrades may ruin it. We rather
recomment the use of Virtual Environments. If you are not familiar with this
term, please check the official documentation by visiting the link below:

    https://docs.python.org/3/tutorial/venv.html

    or

    http://docs.python-guide.org/en/latest/dev/virtualenvs/

Another option is to install **Conda**, a Virtual Environment Manager, or
**AstroConda**, the same but for astronomers. Everything you need to know
about installing both can be found in the link below:

    https://astroconda.readthedocs.io/

.. warning::

    You may find that `ccdproc` and `astroplan` do not come with Astroconda.
    They are not available on any Conda channel either. That means that you will
    have to install them separately. You can do so by downloading the source files
    and installing them by hand, or simply
    `activate your Virtual Environment <https://conda.io/docs/user-guide/tasks/manage-environments.html#activating-an-environment>`_ and
    then install these two packages using pip with

    ``pip install ccdproc astroplan``

System installation is not recommended because it can mess things up specially in
Linux and Mac OS. Before you proceed, make sure that your system has all
the required libraries, as described in `Setup for local installation`_.

Once you have Python running and all the libraries installed either using
Conda/AstroConda or not, you may download the last version available in the
following address:

    https://github.com/soar-telescope/goodman/releases/latest

Before continuing, make sure that your Virtual Environment is active if this is
the case. There are several ways of doing this but normally the command below
should work:

    ``$ source activate <my_environment_name>``

Where ``<my_environment_name>`` is the name of your Virtual Environment (e.g.
astroconda).

Now you can finally install the |pipeline name|. Download the file, decompress
it, and enter the directory created during the file decompression. Test the
installation by typing:

    ``$ python setup.py test``

If you have any errors, check the traceback. If you find difficulties carring
on at this poing, you may contact us by `opening a new issue <https://github.com/soar-telescope/goodman/issues>`_     or using the e-mail
`goodman-pipeline@ctio.noao.edu`.

If no error messages start popping up in your screen, you are good to carry
on with the installation.

    ``$ python setup.py install``

.. note::

    This will install the pipeline in the currently active Python version.
    If you have Virtual Environments, make sure that they are active. If not,
    you can add the ``--user`` option to install only for your user and avoid
    needing root access.

DCR (optional)
~~~~~~~~~~~~~~

.. admonition:: Acknowledgement Note

  Please cite: Pych, W., 2004, PASP, 116, 148

In terms of cosmic ray rejection we shifted to a non-python package because the
results were much better compared to LACosmic's implementation in Astropy.
LACosmic was not designed to work with spectroscopy.

The latest version of the Goodman Spectroscopic Pipeline uses a modified version
of ``dcr`` to help with the pipeline's workflow. It is included under

  ``<path_to_download_location>/goodman/pipeline/data/dcr-source/dcr/``

``goodman`` is the folder that will be created once you untar or unzip the latest
release of the |pipeline name|.

.. important::

    The changes we made to DCR include deletion of all ``HISTORY`` and ``COMMENT`` keywords,
    which we don't use in the pipeline. And addition of a couple of custom
    keywords, such as: ``GSP_FNAM``, which stores the name of the file being
    created. ``GSP_DCRR`` which stores the reference to the paper to cite.


You are still encouraged to visit the official
`Link <http://users.camk.edu.pl/pych/DCR/>`_. We remind again that users of the
Goodman Pipeline should cite the DCR paper with the reference indicated above.

Compiling DCR
-------------

Compiling ``dcr`` is actually very simple.

  ``cd <path_to_download_location>/goodman/pipeline/data/dcr-source/dcr/``

Then simply type:

  ``make``

This will compile `dcr` and also it will create other files. The executable
binary here is ``dcr``.

We have successfully compiled *dcr* right out the box in several platforms, such as:

- Ubuntu 16.04
- Centos 7.1, 7.4
- MacOS Sierra
- Solaris 11


Installing the DCR binary
-------------------------

This is a suggested method. If you are not so sure what you are doing, we
recommend you follow the steps shown below. If you are a more advanced user and
you want to do it your own way, all you have to achieve is to have the ``dcr``
executable binary in your ``$PATH`` variable.

1. Open a terminal
2. In your home directory create a hidden directory ``.bin`` (Home directory
   should be the default when you open a new terminal window)

   ``mkdir ~/.bin``

3. Move the binary of your choice and rename it ``dcr``. If you compiled it,
   most likely it's already called ``dcr`` so you can ignore the renaming part of
   this step.

   ``mv dcr.Ubuntu16.04 ~/.bin/dcr``

   Or

   ``mv dcr ~/.bin/dcr``

4. Add your ``$HOME/.bin`` directory to your ``$PATH`` variable. Open the file
   ``.bashrc`` and add the following line.

   ``export PATH=$PATH:/home/myusername/.bin``

   Where ``/home/myusername`` is of course your home directory.

5. Close and reopen the terminal or load the ``.bashrc`` file.

    ``source ~/.bashrc``
