.. _install:

Install
#######
Using the pipeline remotely is the recommended method, in which case you don't need
to worry about software requirements.

However, for users who wish to go ahead with a local installation, we provide
simple instructions in the current section.

Requirements
************

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

Using Conda
***********

We **do not** recommend the installation of these libraries or the
|pipeline name| in your system since updates and upgrades may ruin it. We rather
recommend the use of Virtual Environments. If you are not familiar with this
term, please check the official documentation by visiting the links below:

    https://docs.python.org/3/tutorial/venv.html

    or

    http://docs.python-guide.org/en/latest/dev/virtualenvs/

Another option is to install **Conda**, a Virtual Environment Manager, or
**AstroConda**, the same but for astronomers. Everything you need to know
about installing both can be found in the link below:

    https://astroconda.readthedocs.io/


.. include:: _working_with_virtualenv.rst

Using PIP
*********

.. warning::

    You may find that `ccdproc` and `astroplan` do not come with Astroconda.
    They are not available on any Conda channel either. That means that you will
    have to install them separately. You can do so by downloading the source files
    and installing them by hand, or simply
    `activate your Virtual Environment <https://conda.io/docs/user-guide/tasks/manage-environments.html#activating-an-environment>`_ and
    then install these two packages using pip with

    ``pip install ccdproc astroplan``



Setup for local installation
****************************
System installation is not recommended because it can mess things up specially in
Linux and Mac OS. Before you proceed, make sure that your system has all
the required libraries, as described in `Requirements`_.

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
    If you have Virtual Environments, make sure that it is active. If not,
    you can add the ``--user`` option to install only for your user and avoid
    needing root access.

.. include:: _install_dcr.rst

.. include:: _shortcuts.rst