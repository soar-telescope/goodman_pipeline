.. _dev-vnc-setup:

Setting up VNC Servers
######################

The SOAR Telescope Data Reduction Servers run Centos 7, below are the
instructions to set it up. The original instructions where obtained from
`this website <https://www.howtoforge.com/vnc-server-installation-on-centos-7>`_

Login as root or administrative user
************************************

A non-root administrative user was provided for basic server management, I would
recommend not going further than installing necessary programs or system
configurations such as the one described in this document. For any other aspects
please contact CISS. For passwords please contact your supervisor to get the
appropriate documents.

    ``ssh sdevel@soardataX``

Where ``X`` can be ``1``, ``2`` or ``3``

Install GNOME
*************
For this application we have chosen GNOME which can be installed with the
following command.

    ``sudo yum groupinstall "GNOME Desktop"``

Answer ``yes`` or ``y`` to approve the installation

Install Tigervnc
****************

   ``sudo yum install tigervnc-server``

Configure the VNC user
**********************

For this example we will be using the username ``goodman``.

    ``sudo useradd goodman``

In order to add a password, do:

    ``sudo passwd goodman``

Type the password for the ``goodman`` user and then verify it by typing it again.


Copy the default vncserver configuration file with the following command.

    ``sudo cp /lib/systemd/system/vncserver@.service \``
        ``/etc/systemd/system/vncserver@:1.service``

And now edit it with your favorite text editor. My preferred text editor for a
terminal is ``vim``.

    ``sudo vim /etc/systemd/system/vncserver@:1.service``

replace all ``<USER>`` appearances by ``goodman``

.. literalinclude:: files/vncserver-service
    :linenos:
    :emphasize-lines: 42-43

So lines 42 and 43 will look like this now:


.. literalinclude:: files/vncserver-1.service
    :lines: 42-43
    :lineno-start: 42

Configuring the Firewall
************************

If you have a firewall you need to do the following commands:

    ``firewall-cmd --permanent --zone=public --add-service vnc-server``

    ``firewall-cmd --reload``

Setting the VNC password
************************
Login as ``goodman``

    ``su goodman``

And run vncserver, this will ask you the VNC password, don't confuse it with
your `user` password.

    ``vncserver``

Type your password and verify it by retyping it when requested.

Enable Service
**************

Exit the ``goodman`` user

    ``exit``

Now you should be logged in as the ``sdevel`` user.

    ``systemctl daemon-reload``

    ``systemctl enable vncserver@:1.service``

    ``reboot``

You will loose conection after this, wait for a reasonable time (1 minute or so)
and then log in again as ``sdevel``

    ``ssh sdevel@soardataX``

Remember that ``X`` can be ``1``, ``2`` or ``3`` and should be the same that you
used at the beginning.

    ``systemctl start vncserver@:1.service``

Log in using your new vnc account
*********************************

Using your vnc client of choice try to connect to the server. I will be using
``vncviewer`` for this example:

   ``vncviewer soardataX:1``

Remember that ``X`` has to be the same of the beginning (``1``, ``2`` or ``3``).