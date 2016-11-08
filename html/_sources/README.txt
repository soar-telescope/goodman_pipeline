Goodman High Throughput Spectrograph
====================================

**WARNING** This code is being developed and is not ready for scientific
use. Always check the branch **other than master**.

If you are interested in this software

The Goodman High Throughput Spectrograph is an imaging spectrograph if
you wish to know more about the instrument please check the `SOAR
website <http://www.ctio.noao.edu/soar/content/goodman-high-throughput-spectrograph>`__

To see full documentation please go to the GitHub hosted site for
`Goodman <https://simontorres.github.io/goodman/>`__

What is contained in this package?
----------------------------------

This repository contains tools for spectroscopy, but after the data is
reduced, i.e. bias and flat corrected, for that part please use `David
Sanmartim's <https://github.com/dsanmartim/goodman_ccdreduction>`__
github repository.

redspec.py
^^^^^^^^^^

Spectra extraction and wavelength calibration

How to use it?
--------------

Place the files in your system's **$PATH** variable and then you can
call it from anywhere in your system

Important Notes
---------------

Needs python2.7 and a newer version of numpy1.12.0 otherwise there will
be problems with numpy.linspace

Development Status
==================

The pipeline is organized in three files (this might increase):

redspec.py
----------

This module basically organizes the data and sort of makes a plan for
processing the night.

-  [x] Ability to operate Blue and Red camera. (transparent for the
   user)
-  [x] Read and Parse Arguments (uses argparse). Also checks for
   consistency.
-  [x] Define Night Class. Reads header information from and creates a
   class that contain important *night* information.
-  [ ] Set mode to be used at the telescope (while observing). Is not
   defined yet if this is really necessary.
-  [ ] Organize Full Night. Edits the Night Class by adding SienceObject
   wich will be defined depending on one of the following modes.
-  [x] Mode 0: One lamp for all targets in the night.
-  [x] Mode 1: One or more lamps per science target.
-  [ ] Mode 2: A text file defines what lamps will be used in what
   targets.
-  [ ] Mode 3: No lamps, solution will be calculated from skylines.

For every ScienceObject the module Process is called wich will do the
actual work.

process.py
----------

In order to work this needs the source's path and a ScienceObject
(class).

-  [x] Read target data and header.
-  [ ] For any of the modes that use comparison lamps:
-  [x] Read lamps data and header
-  [x] Identify single targets (one target in the slit)
-  [ ] Identify multiple targets (most of the work done but not tested
   recently)
-  [x] Trace spectrum.
-  [x] Extract data.

   -  [x] Normal Extraction
   -  [ ] Optimal extraction

-  [ ] Nothing done for the case when skylines will be used as
   wavelength calibrators.

After the extraction of the spectrum it will be packed and then parsed
to the next module.

wavelengthCalibration.py
------------------------

**TODO**: find a better name

-  [x] Interpret science pack (parsed from previous module)
-  [ ] Automatic Wavelength Calibration
-  [x] Identify lines
-  [x] Estimate wavelength solution from information in the header
-  [ ] Find wavelength solution (can be tricky and will take some time)
-  [ ] Automatic Wavelength Calibration by template (can also be tricky
   since dispersion is not linear)
-  [x] Find wavelength solution interactively.
-  [x] Linearize spectrum
-  [x] Put wavelength solution in the header and write down files.

Documentation Status Estimate
=============================

In **Bold** are the classes and in *italics* are their respective
methods. The values presented here are based only on my personal
perception and satisfaction about the status of progress The class
percentage doesn't take in account the completeness of

redspec.py
----------

-  **MainApp** 100%

-  \_ **init** \_ 0%
-  *get*\ args\_ 100%
-  *set*\ night\_ 100%
-  *organize*\ full\_night\_ 100%
-  *convert*\ time\_ 100%
-  *ra*\ dec\_to\_deg\_ 100%
-  *print*\ spacers\_ 100% - Most likely will be removed
-  *print*\ progress\_ 100% - Most likely will be removed

-  **Night** 70%
-  \_ **init** \_ 0%
-  *add*\ sci\_ 100%
-  *add*\ lamp\_ 100%
-  *add*\ sci\_object\_ 100%
-  *is*\ telescope\_ 100%

-  **ScienceObject** 100%
-  \_ **init** \_ 0%
-  *add*\ lamp\_ 100%
-  *print*\ all\_ 100%

process.py
----------

-  **Process** 100%

-  \_ **init** \_ 0%
-  *identify*\ spectra\_ 100%
-  *trace* 100%
-  *mask* 100%
-  *extract* 100%

-  **IdentifiedTarget** 50%
-  \_ **init** \_ 0%

wavelength.py
-------------

-  **WavelengthCalibration** 0%

-  \_ **init** \_ 0%
-  \_ **call** \_ 0%
-  *get*\ wsolution\_ 0%
-  *get*\ calibration\_lamp\_ 0%
-  *get*\ wavelength\_solution\_ 0%
-  *get*\ line\_limits\_ 100%
-  *get*\ line\_centers\_ 100%
-  *get*\ spectral\_characteristics\_ 0%
-  *interpolate* 100%
-  *recenter*\ line\_by\_model\_ 0%
-  *recenter*\ line\_by\_data\_ 0%
-  *predicted*\ wavelength\_ 0%
-  *automatic*\ wavelength\_solution\_ 0%
-  *pixel*\ axis\_cross\_correlate\_ 0%
-  *interactive*\ wavelength\_solution\_ 0%
-  *on*\ click\_ 0%
-  *find*\ more\_lines\_ 0%
-  *update*\ clicks\_plot\_ 0%
-  *plot*\ raw\_over\_reference\_ 0%
-  *evaluate*\ solution\_ 0%
-  *fit*\ pixel\_to\_wavelength\_ 0%
-  *linearize*\ spectrum\_ 0%
-  *add*\ wavelength\_solution\_ 0%


