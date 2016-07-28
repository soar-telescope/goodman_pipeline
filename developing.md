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
following answers to your needs. (some parts are omitted)

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

### Final Thoughts

These tools are highly customizable so expect some troubles setting this 
up. Keep your eyes open and read very well the debug or feedback on the
problem.