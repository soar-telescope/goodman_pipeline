|Build Status| |Coverage Status| |Documentation Status| |astropy|

The Goodman High Throughput Spectrograph (Goodman HTS) Data-Reduction
Pipeline is the SOAR Telescope's official data reduction pipeline for
*Goodman HTS*.

It has been fully developed in Python 3.5 and uses mostly astropy
affiliated packages with the exception of
`dcr <http://users.camk.edu.pl/pych/DCR/>`__ which is an external tool
that does cosmic ray identification and correction. The reason for using
it instead of LACosmic is that it works very well for spectroscopic data
and the results are evidently superior. Some of the negative aspects of
using this external (meaning outside of Python domains) software were:
The integration into the pipeline's workflow and the use of an external
``dcr.par`` parameter file. Such parameters have to be changed by hand
and can't be integrated into the pipeline's workflow itself. In
particular for binning 2x2 and custom ROI those parameters contained in
*dcr.par* has to be specifically tuned.

The pipeline is currently in its first beta release version. We have
thoroughly tested it but we know there are still bugs to be found. If
you found one or have any comment please contact Simón Torres at storres
[at] ctio noao edu.

The Goodman High Throughput Spectrograph is an imaging spectrograph, if
you wish to know more about the instrument please check the `SOAR
website <http://www.ctio.noao.edu/soar/content/goodman-high-throughput-spectrograph>`__

What is contained in the project?
---------------------------------

The Goodman HTS Pipeline is distributed as a single package with two
main scripts, *redccd* dedicated to the basic ccd reduction process and
*redspec* for advanced processing of spectroscopic data such as target
identification, extraction and wavelength calibration.

They are supposed to work in sequence and they are called from two
indepentent scripts, *redccd* and *redspec*.

Other files included are:
- Readme.md
- requirements.txt
- user\_manual.pdf
- dcr binaries
- Library of Wavelength Calibrated lamps.

How to install it?
------------------

In principle you don't need to install it since SOAR Telescope is
providing dedicated servers where you can run the pipeline.

If you are a *user* you should follow the instructions in the
``user_manual.pdf`` file, where you will find instruction for accessing
the server and using the pipeline, even for installing it.

The user manual file is included in the distribution.

If you want to try it out as a developer please request the developer's
manual.

Common Failures
===============

-  The pipeline fails to detect a target: Usually is the keyword OBSTYPE
   is wrong and is trying to detect targets where there is any.

-  The wavelength solution RMS Error is too high (larger than 3): The
   lamp was mislabeled. For instance it says HgArNe when is HgAr.

Acknowledge
===========

Contributors:
-------------

`David Sanmartim <https://github.com/dsanmartim>`__ developed the first
version of the ccd reduction part before he moved to a new position.
That first version was included in the pipeline as a package and has
evolved along all the code.

`Tina Armond <https://github.com/tarmond>`__ Helped with the
identification of comparison lamps for the *reference lamps library*.

Citations:
----------

This pipeline makes extensive use of Astropy therefore you should cite
as suggested on `Astropy Citation
Page <https://github.com/astropy/astropy/blob/master/CITATION>`__ as
follows:

::

    This research made use of Astropy, a community-developed core Python package
    for Astronomy (Astropy Collaboration, 2013, 2018).

It also uses `DCR <http://users.camk.edu.pl/pych/DCR/>`__ for cosmic
rays identification and removal. You should cite `this
paper <http://adsabs.harvard.edu/abs/2004PASP..116..148P>`__

::

     Pych, W., 2004, PASP, 116, 148

.. |Build Status| image:: https://travis-ci.org/soar-telescope/goodman.svg?branch=master
   :target: https://travis-ci.org/soar-telescope/goodman
.. |Coverage Status| image:: https://coveralls.io/repos/github/soar-telescope/goodman/badge.svg?branch=master
   :target: https://coveralls.io/github/soar-telescope/goodman?branch=master
.. |Documentation Status| image:: https://readthedocs.org/projects/goodman/badge/?version=latest
   :target: http://goodman.readthedocs.io/en/latest/?badge=latest
.. |astropy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
   :target: http://www.astropy.org/
