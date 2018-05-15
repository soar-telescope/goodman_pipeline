.. _dev-sphinx-docs:

Documentation using Sphinx
##########################

We have chosen to use Google's style for writing the Python Docstrings
because is easy to read in the code and also allows to generate
documentation in different formats. The guidelines for writing
Docstrings are stated in the `PEP257 <https://www.python.org/dev/peps/pep-0257/>`_
But no style is defined in there.

By default, Sphinx uses *reStructuredText* for Python's Docstrings. But it
is a bit hard to read, for instance:

.. code-block:: python

    def some_function(name, state):
        """This function does something.

        :param name: The name to use.
        :type name: str.
        :param state: Current state to be in.
        :type state: bool.
        :returns: int -- the return code.
        :raises: AttributeError, KeyError

        """
        return True

The same function with Google Style Docstrings would be:

.. code-block:: python

    def some_function(name, state):
        """This Function does something

        Args:
            name (str): The name to use.
            state (bol): Current State to be in.

        Returns:
            int: The return code.

        Raises:
            AttributeError: Some description
            KeyError: Some other description

        """
        return True

Everyone in our team agrees that Google style is easier to read. Sphinx's style
produces better formatted output and gives you more control on
formatting.

The good news is that it is possible to produce documentation using
Google Style and Sphinx. There are a few necessary steps to achieve
this:

Install Sphinx
**************

1. Install sphinx.

    ``pip install sphinx``


2. Install Napoleon Plugin. Allows to use Google style with Sphinx.

    ``pip install sphinxcontrib.napoleon``


3. Go to your source folder and run ``sphinx-quickstart`` and adapt the
 following answers to your needs. (some parts were omitted)

    .. code-block:: bash

        $ sphinx-quickstart
        Welcome to the Sphinx 1.4.5 quickstart utility.

        Please enter values for the following settings (just press Enter to
        accept a default value, if one is given in brackets).

        Enter the root path for documentation.
        > Root path for the documentation [.]: docs

        You have two options for placing the build directory for Sphinx output.
        Either, you use a directory "_build" within the root path, or you separate
        "source" and "build" directories within the root path.
        > Separate source and build directories (y/n) [n]: y

        Inside the root directory, two more directories will be created; "_templates"
        for custom HTML templates and "_static" for custom stylesheets and other static
        files. You can enter another prefix (such as ".") to replace the underscore.
        > Name prefix for templates and static dir [_]:

        The project name will occur in several places in the built documentation.
        > Project name: MyProject
        > Author name(s): The Author
        ...
        > Project version: 1.0
        > Project release [1.0]:
        ...
        > Create Makefile? (y/n) [y]:


4. A file ``conf.py`` has been created in the folder docs. In there you need to change some settings:

    - Add ``sphinxcontrib.napoleon``
    - Insert ``$PATH`` variable into the system
    - Add napoleon settings


    .. code-block:: python

        # conf.py

        # Add autodoc and napoleon to the extensions list
        extensions = ['sphinx.ext.autodoc', 'sphinxcontrib.napoleon']

        # Add source directory to PATH (where the modules are)
        import sys
        sys.path.append('/path/to/source')

        # Napoleon settings
        napoleon_google_docstring = True
        napoleon_numpy_docstring = True
        napoleon_include_init_with_doc = False
        napoleon_include_private_with_doc = False
        napoleon_include_special_with_doc = False
        napoleon_use_admonition_for_examples = False
        napoleon_use_admonition_for_notes = False
        napoleon_use_admonition_for_references = False
        napoleon_use_ivar = False
        napoleon_use_param = True
        napoleon_use_rtype = True
        napoleon_use_keyword = True

5. Use ``sphinx-apidoc`` to build your documentation.

    ``sphinx-apidoc -f -o docs/source projectdir``

    Or as I did from the home directory:

    ``sphinx-apidoc -f -o ~/development/soar/goodman/docs/source \``
        ``~/development/soar/goodman``


6. Make your documentation: Go to ``~/development/soar/goodman/docs`` and run the
command below. (There are more options that you can use).

    ``cd ~/development/soar/goodman/docs``

    ``make html``


Configure PyCharm to do it
**************************

7. In PyCharm go to ``Run > Edit Configurations...`` and then click
   the green **+** symbol in the upper left corner
   *(Pycharm Community Edition 2016.2)* and choose
   ``Python Docs > Sphinx Task``

8. A new dialog will open and it should look like this:

   ![docs-configuration](./docs/img/configurations.png)

   Apply or click OK. In the upper right corner of pycharm you should find
   the option for ``build-documentation`` if you used the same name for this
   particular configuration.

Local Management - How to use github pages
******************************************

.. warning::

    This section is incomplete if you need assistance in this regard
    please contact me

In oder to minimize the chance of making a mistake I suggest the following
configuration.

Just to be sure we are talking the same, I want to remind you that the symbol
``~`` represent your home directory, for instance, if your user name is ``user``
its home directory will be ``/home/user``. In GNU/Linux, this is by default but
beware that it might change, a simple *trick* to find out your home directory, in
case you are in doubt, is to do the following:

Open a terminal and type these two commands

    ``cd``

By default this will take you to your home directory

    ``pwd``

This command does *Print Working Directory*


Create Directory Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^
I suggest using separate directories for *code development* and for *building
documentation*. The command sequence is as follows.

For the **code**, create the folder  ``~/development/soar/``

    ``mkdir -p ~/development/soar/``

For the **documentation** is similar:

    ``mkdir -p ~/documentation/soar``

Clone the repository
^^^^^^^^^^^^^^^^^^^^
At this point we need to run ``git clone`` in both location.

Go to the ``code`` folder:

    ``cd ~/development/soar/``

Clone the repository, for this you need to setup your github account first.

    ``git clone git@github.com:soar-telescope/goodman.git``

You will have now a folder named ``goodman``

Now go to the ``documentation`` folder and repeat the ``clone`` command:

    ``cd ~/documentation/soar/``

    ``git clone git@github.com:soar-telescope/goodman.git``

Again you will have the ``goodman`` directory but there is an extra step to
do here.

Go into the ``goodman`` directory and checkout the branch ``gh-pages``

    ``cd goodman``

    ``git checkout gh-pages``

This will bring all the files related to `GitHub Pages <https://pages.github.com/>`_

Inside your code development directory is a folder called ``docs``. Here
is where you run ``sphinx-quickstart``. After this it will contain a file
named ``Makefile`` and two folders ``source`` and ``build``


Modifying the Makefile
^^^^^^^^^^^^^^^^^^^^^^

Modify ``Makefile`` to build the documentation in a different location. Open
the file using your favourite editor and change the ``BUILDDIR`` variable to:

    ``BUILDDIR      = /HOME/USER/documentation/soar/goodman``

It should look like this: (Note that I changed the name of the home directory
for security reasons, change ``/HOME/USER`` by whatever matches your system)

::

    # Minimal makefile for Sphinx documentation
    #

    # You can set these variables from the command line.
    SPHINXOPTS    =
    SPHINXBUILD   = sphinx-build
    SPHINXPROJ    = GoodmanPipelines
    SOURCEDIR     = source
    BUILDDIR      = /HOME/USER/documentation/soar/goodman

    # Put it first so that "make" without argument is like "make help".
    help:
            @$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

    .PHONY: help Makefile

    # Catch-all target: route all unknown targets to Sphinx using the new
    # "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
    %: Makefile
            @$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)


How to include Markdown README to index
***************************************

Although I will recommend to write an .rst file instead, I will include the
necessary steps to include a Markdown README. I didn't get the results I wanted
with Markdown.

1. Install ``recommonmark``

    ``sudo pip2.7 install recommonmark``

2. Edit Sphinx's ``conf.py`` and add the following lines.

    .. code-block:: python

        from recommonmark.parser import CommonMarkParser

        source_suffix = ['.rst', '.md']

        source_parsers = {
            '.md' : CommonMarkParser,
        }

4. go to ``docs/sources`` and add the following line to ``index.rst``:

    ::

        .. include:: README.md


    So it will look like this:

    ::

        Welcome to Goodman Spectroscopic Tools's documentation!
        ===========================================================

        .. include:: README.md

        Contents:

        .. toctree::
           :maxdepth: 2
        
        Indices and tables
        ======================

        * :ref:`genindex`
        * :ref:`modindex`
        * :ref:`search`


Convert from Markdown to reStructuredText
*****************************************

But you can also convert your Markdown README to .rst using ``pandoc``.

    ``sudo yum install pandoc``

And you can use it with:

    ``pandoc --from=markdown --to=rst --output=README.rst README.md``

Converting Markdown to .rst worked better for me. Also I would recommend some
editing since you don't want all the README information to go to the web page,
such as the link to the page itself for instance. Also some checks are needed to
make sure the translation from *Markdown* to *reStructuredText* was correct.


Building the ``html`` Documentation
***********************************

Use ``make html`` in your documentation's root folder (where the Makefile is located)

    ``cd ~/development/soar/goodman/docs``

    ``make html``

The new html documentation will be located at
``~/documentation/soar/goodman/html``. You need to move all its content to the
parent directory i.e. ``~/documentation/soar/goodman``

    ``cd ~/documentation/soar/goodman``

    ``cp -vr html/* .``

6. Finally in add commit and push

Add all new files to the repository

    ``git add *``

Commit the new changes. Try to use a meaningful message.

    ``git commit -m "updating html documentation"``

Push to the repository

    ``git push``

Then you can see the page in
`https://soar-telescope.github.io/goodman/ <https://soar-telescope.github.io/goodman/>`_


Final Thoughts
**************

These tools are highly customizable so expect some troubles setting this
up. Keep your eyes open and read very well the debug or feedback on the
problem.