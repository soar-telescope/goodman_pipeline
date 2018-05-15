.. _dev-requirements:
How to create a ``requirements.txt`` file
#########################################

.. warning::

    DO THIS ALWAYS BEFORE A RELEASE

The file ``requirements.txt``  is necessary to facilitate other users the task
of setting up the system to get your *Shiny Python Code* working for them.
It is basically a list of libraries with their respective version numbers.

We came across a solution to do this automatically which is:

    ``pip freeze > requirements.txt``

This works fine but it saves ALL the libraries in the environment and not only
for the project (this might be useful to save the environment though).

I found a better way to do it, also automatically:

    ``sudo pip install pipreqs``

Once installed do:

    ``pipreqs /path/to/the/repository``

For the Goodman Pipeline the list also contained ``pygtk`` but it fails.
Simply delete it. PyGTK should come by default in most systems.
