.. _dev-build-manuals:

User and Developer Manuals
##########################

The user and developer manuals are created using Sphinx. The source files are
written in *reStructuredText* a markup language with simple syntaxis that allows
for rich formatting. You can start `here <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`_
Although the language's syntaxis is simple you need to familiarize with it, or
you might get scared away, but don't worry, is really easy.

Preparation
***********

I recommend using `PyCharm <https://www.jetbrains.com/pycharm/>`_ I will not
provide instrucctions for installing it but is easy and you can find
instructions in their site. Also because the installation depends on the platform.

You will also need a `GitHub <https://github.com/>`_ account and very important
`Register your public ssh key <https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/>`_
on GitHub.

And finally define a directory structure, here I present my suggestion but feel
free to ignore this in case you already have yours. Keep in mind here that the
most important here is to have *something that works* for you.

    ``mkdir -p ~/documentation/soar``

    ``cd ~/documentation/soar``

Also you will need to install ``rst2pdf``


Get the source Files
********************

The source files are stored on GitHub.com, in order to obtain them simply do:

    ``git clone git@github.com:soar-telescope/goodman_user_manual.git``

All the repository will be downloaded to ``goodman_user_manual``

Open *PyCharm* and go to ``File > Open`` Navigate to
``~/documentation/soar/goodman_user_manual``
and click *Open*

.. image:: img/open_project.png

Choose *Open in a new window* when asked.

Editing the Manuals
*******************

This `rst cheatsheet <https://github.com/ralsina/rst-cheatsheet/blob/master/rst-cheatsheet.rst>`_
seems a very good resource to start. I will not go into further details here.

Building the PDF Files
**********************

In order to produce the pdfs, first you need to configure *PyCharm*. Go to
``Run > Edit Configurations...``

- Click the green plus sign and select ``Python Docs... > Sphinx task``.

- Edit the fields so they will look like this:

    - For the Developers Manual

        ``Name: create_pdf_dev_manual``

        ``Command: pdf``

        ``Input: ~/documentation/soar/goodman_user_manual/dev_manual/``

        ``Output: ~/documentation/soar/goodman_user_manual``

    - For the Users Manual:

        ``Name: create_pdf_user_manual``

        ``Command: pdf``

        ``Input: ~/documentation/soar/goodman_user_manual/user_manual/``

        ``Output: ~/documentation/soar/goodman_user_manual``


.. image:: img/configurations.png


Commit to GitHub
****************

Once you have completed a feature click this button |commit| to do a commit to
your repository.


.. |commit| image:: img/vcs_up.png
    :align: middle
    :width: 40

This will open a pop up window. In the upper left region you see the modified
files, below this you type a *Commit Message*, this has to be precise, not
too short not too long, just a reminder of what you did. At the bottom you
can see the changes you are introducing in the file selected in the top
left section of the window.

.. image:: img/commit_screen.png

If you select *Commit and Push* this will update the local repository as well as
the remote, i.e, the GitHub hosted version. If you choose *Commit*  only the
changes will not be uploaded. Later you can go to ``VCS > Git > Push..``.
