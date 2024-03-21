.. _dcr:

Installing DCR
**************

.. admonition:: Acknowledgement Note

  Please cite: Pych, W., 2004, PASP, 116, 148

In terms of cosmic ray rejection we shifted to a non-python package because the
results were much better compared to LACosmic's implementation in Astropy.
LACosmic was not designed to work with spectroscopy. Though since version
:ref:`1.1.0 <v1.1.0>` we shifted from Astropy to Astroscrappy's implementation
of LACosmic.

The latest version of the Goodman Spectroscopic Pipeline uses a modified version
of ``dcr`` to help with the pipeline's workflow. It is included under

  ``<path_to_download_location>/goodman_pipeline/goodman_pipeline/data/dcr_source/dcr/``

``goodman_pipeline-<version>`` is the folder that will be created once you untar or unzip the latest
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
^^^^^^^^^^^^^

Compiling ``dcr`` is actually very simple.

  ``cd <path_to_download_location>/goodman_pipeline/goodman_pipeline/data/dcr_source/dcr/``

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
^^^^^^^^^^^^^^^^^^^^^^^^^

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

.. include:: _shortcuts.rst