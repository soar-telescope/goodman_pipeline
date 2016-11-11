# Development Notes

## Notes on how to build the Documentation using Sphinx

We have chosen to use Google's style for writing the Python Docstrings
because is easy to read in the code and also allows to generate
documentation in different formats. The guidelines for writing
Docstrings are stated in the [PEP257](https://www.python.org/dev/peps/pep-0257/)
But no style is defined in there.

By default, Sphinx uses its own format for Python's Docstrings. But it
is a bit hard to read, for instance:

```python
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

```

The same function with Google Style Docstrings would be:

```python
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
```

Everyone agrees that Google style is easier to read. Sphinx's style
produces better formatted output and gives you more control on
formatting.

The good news is that it is possible to produce documentation using 
Google Style and Sphinx. There are a few necessary steps to achieve
this:

### Install Sphinx

1. Install sphinx.

   ```shell
   $ pip install sphinx
   ```

2. Install Napoleon Plugin. Allows to use Google style with Sphinx.

   ```shell
   $ pip install sphinxcontrib.napoleon
   ```

3. Go to your source folder and run `sphinx-quickstart` and adapt the
following answers to your needs. (some parts were omitted)

   ```shell
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
   ```

4. A file `conf.py` has been created in the folder docs. In there
   you need to change some settings:
      1. Add sphinxcontrib.napoleon
      2. Insert $PATH variable into the system
      3. Add napoleon settings

   ```python
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
   ```

5. Use `sphinx-apidoc` to build your documentation.
   ```shell
   $ sphinx-apidoc -f -o docs/source projectdir
   ```

6. Make your documentation. There are more options for the output.
   ```shell
   $ make html
   ```

### Configure PyCharm to do it.

7. In PyCharm go to `Run > Edit Configurations...` and then click
   the green **+** symbol in the upper left corner 
   _(Pycharm Community Edition 2016.2)_ and choose
   `Python Docs > Sphinx Task`
   
8. A new dialog will open and it should look like this:
   
   ![docs-configuration](./docs/img/configurations.png)
   
   Apply or click OK. In the upper right corner of pycharm you should find
   the option for `build-documentation` if you used the same name for this
   particular configuration.

## Local Management - How to use github pages

**This section is incomplete if you need assistance in this regard
please contact me**

In oder to minimize the chance of making a mistake I set up the following
configuration.

The code development directory is in `~/development/soar/goodman` and
the documentation is in `~/documentation/soar/goodman`. The git
repository is basically the same with the only difference being that in 
the documentation folder I checkout `gh-pages` and in development I 
work with the `master` branch or any of the code branches.
 
### Procedure
1. Create directories and clone from github
    ```shell
    $ mkdir ~/development/soar
    $ cd !*
    $ git clone git@github.com:simontorres/goodman.git
    $ 
    $ cd
    $ mkdir ~/documentation/soar/
    $ cd !*
    $ git clone git@github.com:simontorres/goodman.git
    $ cd goodman
    $ git checkout gh-pages
    ```

2. Inside your code development directory is a folder called `docs`. Here
is where you run `sphinx-quickstart`. After this it will contain a file
named `Makefile` and two folders `source` and `build`

### How to include Markdown README to index

Although I will recomend to write an .rst file instead, I will include the
necessary steps to achieve that. I didn't get the results I wanted with
Markdown.

1. Install recommonmark
    ```shell
    $ sudo pip2.7 install recommonmark
    `
2. Edit Sphinx's `conf.py` and add the following lines.
    ```python
    from recommonmark.parser import CommonMarkParser
    
    source_suffix = ['.rst', '.md']
    
    source_parsers = {
        '.md' : CommonMarkParser,
    }
    ```
 
4. go to `docs/sources` and add the following line to `index.rst`:
    ```shell
    .. include:: README.md
    ```

So it will look like this:

    ```text
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
    ```

5. Use `make html` in your documentation's root folder (where the Makefile is located)

But you can also convert your Markdown README to .rst using `pandoc`.
    
    ```shell
    $ sudo yum install pandoc
    ```
    
And you can use it with:
    
    ```shell
    $ pandoc --from=markdown --to=rst --output=README.rst README.md 
    ```

### Final Thoughts

These tools are highly customizable so expect some troubles setting this 
up. Keep your eyes open and read very well the debug or feedback on the
problem.