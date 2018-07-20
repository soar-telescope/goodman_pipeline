.. _remote-access:

Setup for Remote Use
********************
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
^^^^^^^^^^^^^^^^^^^^^^^^^^
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

.. important::

    Please, help us to create an organized enviroment by creating a new folder
    using the format ``YYYY-MM-DD`` within your institution's directory and
    using it to process your data.

VNC from the Terminal
^^^^^^^^^^^^^^^^^^^^^
Find the ``<display-number>`` that corresponds to you from the `VNC Displays table`_.
Open a terminal, and assuming you have installed ``vncviewer``.

    ``vncviewer <vnc-server>:<display-number>``

You will be asked to type in the ``<password>`` provided.

.. important::

    The real values for ``<vnc-server>`` and ``<password>``
    should be provided by your support scientist.

If the connection succeeds you will see a *Centos 7* Desktop using *Gnome*.


.. include:: _shortcuts.rst