# Goodman High Throughput Spectrograph
**WARNING** This code is being developed and is not ready for
scientific use. Always check the branch **other than master**.
 
 If you are interested in this software 

The Goodman High Throughput Spectrograph is an imaging spectrograph
 if you wish to know more about the instrument please check the 
 [SOAR website](http://www.ctio.noao.edu/soar/content/goodman-high-throughput-spectrograph)
 
To see full documentation please go to the GitHub hosted site for [Goodman](https://simontorres.github.io/goodman/)

## What is contained in this package?

This repository contains tools for spectroscopy, but after the data is 
reduced, i.e. bias and flat corrected, for that part please use 
[David Sanmartim's](https://github.com/dsanmartim/goodman_ccdreduction) github repository.

#### redspec.py  
 Spectra extraction and wavelength calibration

## How to install it?

Either clone or download this code. If you decide to clone it just do ...

Or you can simply go and click _here_ for download a zip file with all the script files.

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
- [ ] Observing Mode **0**: One solution applied to all night
- [ ] Interactive Mode **True**
- [ ] Data Path **./**
- [ ] Destination folder for processed data **./**
- [ ] Search Pattern **fzh.**
- [ ] Output prefix **g**
- [ ] Reference Files Path [/path/to/this/repo/]**refdata/**
- [ ] Plots Enabled **False**
- [ ] Reference Lamp **None**

## Important Notes

Needs python2.7 and a newer version of numpy1.12.0 otherwise there will
be problems with numpy.linspace
