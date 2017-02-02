# Goodman High Throughput Spectrograph
**WARNING** This is the first release of this pipeline, although we have 
thoroughly tested it there might still be bugs. Please let me know by an
e-mail to storres [at] ctio noao edu.

The Goodman High Throughput Spectrograph is an imaging spectrograph,
 if you wish to know more about the instrument please check the 
 [SOAR website](http://www.ctio.noao.edu/soar/content/goodman-high-throughput-spectrograph)
 
To see full documentation please go to the GitHub hosted site for [Goodman](https://simontorres.github.io/goodman/)

## What is contained in this package?

This repository contains tools for processing Goodman's spectrograph data, from the ccd image on. 
It is separated into two main components. The **CCD** part is done by _redccd_ originally developed by 
[David Sanmartim's](https://github.com/dsanmartim/goodman_ccdreduction) and currently maintained by our team. It does
the standard _ccd image reduction_, i.e. trim, bias and flat correction. Currently the ccd reduction pipeline has been
integrated into this package so there is no need to look for it separtely. The **Spectroscopic** processing is done by
_redspec_ and includes the following features:


- [x] Identify targets in images
- [x] Trace the target
- [x] Extract the target with background subtraction
- [x] Find the wavelength Solution (Requires User Input)
- [x] Linearize data
- [x] Write wavelength solution to a FITS header
- [x] Create a new file for the wavelength calibrated 1D spectrum

There is also a library of calibrated spectrum in FITS format. Different configurations are provided but is very easy
to add your own.

## How to install it?

Either clone or download this code. If you decide to clone it, just do:

```shell
git clone https://github.com/simontorres/goodman.git
```

Or you can simply go and [click here](https://github.com/simontorres/goodman/archive/master.zip)
for download a zip file with all the script files.

### Requirements

This software was developed on Python 2.7, use the `requirements.txt` file to install all dependencies.

```shell
sudo -H pip2.7 install -r requirements.txt
```

Remember that if you do not have super user privileges, you can still install all the requirements by adding 
the `--users` flag. There is also an option of [installing it using different virtual environments](http://docs.python-guide.org/en/latest/dev/virtualenvs/) 
or even using [Anaconda](https://www.continuum.io/downloads).


You may have some issues when using MatPlotLib on `virtualenvs` and on Mac OSX. If so, you can try to follow 
the instructions on [this site](http://matplotlib.org/faq/osx_framework.html#osxframework-faq) and, then add the 
following line on your `.bashrc` or `.bash_profile` file.

```bash
alias python='pythonw' 
```

## How to use it?
 
To get a list of the possible arguments do:

```shell
/path/to/this/repo/redspec.py --help
```

The simplest way of running this pipeline is to go to your data folder,
already processed with [goodman_ccdreduction](https://github.com/dsanmartim/goodman_ccdreduction)
and execute `redspec.py`

```shell
/path/to/this/repo/redspec.py
```

Will run the following defaults:
- [ ] Observing Mode **0**: One solution per configuration applied to all night
- [ ] Interactive Mode **True**
- [ ] Data Path **./**
- [ ] Destination folder for processed data **./**
- [ ] Search Pattern **fzh_**
- [ ] Output prefix **g**
- [ ] Reference Files Path [/path/to/this/repo/]**refdata/**
- [ ] Plots Enabled **False**
- [ ] Reference Lamp **None**

