# Goodman High Throughput Spectrograph
**WARNING** This code is being developed and is not ready for
 scientific use. Always check the branch other than master.

The Goodman High Throughput Spectrograph is an imaging spectrograph
 if you wish to know more about the instrument please check the 
 [SOAR website](http://www.ctio.noao.edu/soar/content/goodman-high-throughput-spectrograph)

## What is contained in this package?

This repository contains tools for spectroscopy, but after the data is 
reduced, i.e. bias and flat corrected, for that part please use 
[David Sanmartim's](link) github repository.

#### redspec.py  
 Spectra extraction and wavelength calibration

## How to use it?

Place the files in your system's **$PATH** variable and then you can call it from anywhere in your system

## Important Notes

Needs python2.7 and a newer version of numpy1.12.0 otherwise there will
be problems with numpy.linspace