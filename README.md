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

## How to use it?

Place the files in your system's **$PATH** variable and then you can call it from anywhere in your system

## Important Notes

Needs python2.7 and a newer version of numpy1.12.0 otherwise there will
be problems with numpy.linspace

# Development Status

The pipeline is organized in three files (this might increase):

## redspec.py
This module basically organizes the data and sort of makes a plan for processing the night.

- [x] Read and Parse Arguments (uses argparse). Also checks for consistency.
- [x] Define Night Class. Reads header information from and creates a class that contain important _night_ information.
- [ ] Set mode to be used at the telescope (while observing). Is not defined yet if this is really necessary.
- [ ] Organize Full Night. Edits the Night Class by adding SienceObject wich will be defined depending on one of the following modes.
  * [x] Mode 0: One lamp for all targets in the night.
  * [x] Mode 1: One or more lamps per science target.
  * [ ] Mode 2: A text file defines what lamps will be used in what targets.
  * [ ] Mode 3: No lamps, solution will be calculated from skylines.

For every ScienceObject the module Process is called wich will do the actual work.

## process.py
In order to work this needs the source's path and a ScienceObject (class).

- [x] Read target data and header.
- [ ] For any of the modes that use comparison lamps:
  1. [x] Read lamps data and header
  2. [x] Identify single targets (one target in the slit)
  3. [ ] Identify multiple targets (most of the work done but not tested recently)
  4. [x] Trace spectrum.
  5. [x] Extract data.
- [ ] Nothing done for the case when skylines will be used as wavelength calibrators.

After the extraction of the spectrum it will be packed and then parsed to the next module.

## wavelengthCalibration.py
**TODO**: find a better name

- [x] Interpret science pack (parsed from previous module)
- [ ] Automatic Wavelength Calibration
  - [x] Identify lines
  - [x] Estimate wavelength solution from information in the header
  - [ ] Find wavelength solution (can be tricky and will take some time)
- [ ] Automatic Wavelength Calibration by template (can also be tricky since dispersion is not linear)
- [ ] Find wavelength solution interactively.
- [ ] Put wavelength solution in the header and write down files.


# Documentation Status Estimate
