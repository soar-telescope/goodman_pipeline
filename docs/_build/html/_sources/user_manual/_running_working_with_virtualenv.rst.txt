Working with Virtual Environments
*********************************

Virtual environments are a very useful tool, the main contribution of them being:

- Portability
- Protection to the host environment
- Flexibility

If you know nothing about them we recommend you to start in the `Conda site <https://conda.io/docs/index.html>`_.

For the purpose of this manual we will just say that a *Virtual Environment*
lets you have a custom set of libraries/tools in one place, and most importantly
is independent of your host system. Installation will not be discussed here but
you can visit `this link <https://conda.io/docs/user-guide/tasks/manage-environments.html>`_
for information.

Discover what environments exist in your system.
  ``conda env list``

  Will print a list where the first column is the name.

Activate (enter) the virtual Environment.
  ``source activate <venv-name>``

  Where ``<venv-name>`` is the name of your virtual environment. Your shell's
  prompt will change to:

  ``(<venv-name>) [user@hostname folder-name]$``


Deactivate (leave) the virtual environment.
  ``source deactivate``

  This time the prompt will change again to:

  ``[user@hostname folder-name]$``

