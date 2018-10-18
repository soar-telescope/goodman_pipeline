.. _install:

Install
#######

We do not have the resources to provide installation support, thus we provide
a server with the latest and older versions installed that users with access rights
can use. However, the installation process is simple, and is described below (we hope that in enough detail)
so that users can follow the steps and end up with a running version of the pipeline.

Installation Overview
*********************

The required steps for a successful installation can vary depening on your OS and version, but here is a general
list.

Install *astroconda*.
  Visit `AstroConda site <https://astroconda.readthedocs.io/en/latest/>`_ and follow
  their instructions on how to install and basic setup. Make sure you add
  the `astroconda channel <https://astroconda.readthedocs.io/en/latest/installation.html#configure-conda-to-use-the-astroconda-channel>`_

Create a *virtual environment*.
  This will depend on the method of installation, but here is the official
  documentation on `virtual environment management <https://conda.io/docs/user-guide/tasks/manage-environments.html>`_.

Install all requirements.
  This can be done along with the virtual environment creation, but it might be
  necessary to manually install dependencies.

Run test code.
  This is optional and only possible if you are installing from source code or
  a released version.

Install the goodman pipeline.
  The most important step (of course). Go to installation_

Getting the code
****************

There are several ways to get the goodman pipeline code on your computer.

1. Clone the repository: Cloning the `repository <https://github.com/soar-telescope/goodman_pipeline>`_ will
   get you the latest version of the code, this could be a very unstable version,
   so is not recommended. Also it is very likely that you will download things you don't need.

2. Get a release: A release is like a snapshot of a particular stage of development.
   Though its not easy to guarantee that the code will not break, a released version should
   be stable in the sense that no experimental code should exist. Also, it is a
   special packaging of the necessary parts for operation. The releases are listed
   `here <https://github.com/soar-telescope/goodman_pipeline/releases>`_. make sure
   it has a green tag that says **Latest Release**.

3. Using pip: Since version :ref:`1.1.2 <v1.1.2>` the goodman pipeline can be installed using
   pip: ``pip install goodman-pipeline`` however this does not install astroconda
   neither creates the virtual environment.


Dependencies
************

We have made a significant effort to maintain dependencies under control, with
one exception. No worries! we have :ref:`detailed instruccions <dcr>` for
dealing with it.

For the ordinary dependencies we make use of the great tools provided by astroconda.
In one step we can create a virtual environment with all the dependencies on it.
This does not work for the pip installation, but either by cloning the repository or
by downloading a release you get a file called `environment.yml` that can be used by
conda in the following way:

  ``conda env create -f environment.yml``

This will create a virtual environment called ``goodman_pipeline`` that you can start
using by running.

  ``source activate goodman_pipeline``

.. include:: _install_dcr.rst

Running tests
*************

*Test code is code to test the code*, yes, that's right, and the |pipeline full name|
comes with a lot of it. Ideally, 100% of the code should be covered; we are not
there yet, but we are getting closer. Running the test code is not a requirement
but it helps a lot in identifying problems. Go to the extraction folder
(where the file ``setup.py`` exist) and execute the following instruction:

  ``python setup.py test``

This will detect the test code and execute it. Failure or success will be reported in a text message.

.. _installation:

Installation
************

If the previous steps ended successfully this next part should be easy. Go to the extraction
folder and execute the following instruction:

  ``python setup.py install``


Alternatively, after the creation of the virtual environment you can install the pipeline
using pip.

  ``pip install goodman-pipeline``
