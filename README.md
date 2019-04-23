# Goodman High Throughput Spectrograph Data Reduction Pipeline

[![Build Status](https://travis-ci.org/soar-telescope/goodman_pipeline.svg?branch=master)](https://travis-ci.org/soar-telescope/goodman_pipeline)
[![Coverage Status](https://coveralls.io/repos/github/soar-telescope/goodman_pipeline/badge.svg?branch=master)](https://coveralls.io/github/soar-telescope/goodman_pipeline?branch=master)
[![Documentation Status](https://readthedocs.org/projects/goodman/badge/?version=latest)](http://goodman.readthedocs.io/en/latest/?badge=latest)
[![pypi](https://img.shields.io/pypi/v/goodman_pipeline.svg?style=flat)](https://pypi.org/project/goodman-pipeline/)
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

## Overview
The Goodman High Throughput Spectrograph (Goodman HTS) Data-Reduction Pipeline
is the SOAR Telescope's official data reduction pipeline for *Goodman HTS*.

It has been fully developed in Python 3.5 and uses mostly astropy affiliated packages
with the exception of [dcr](http://users.camk.edu.pl/pych/DCR/) which is an external tool
that does cosmic ray identification and correction. The reason for using it
instead of LACosmic is that it works very well for spectroscopic data and the
results are evidently superior. Some of the negative aspects of using this
external (meaning outside of Python domains) software were: The integration into
the pipeline's workflow and the use of an external `dcr.par` parameter file.
 Such parameters have to be changed by hand and can't be integrated into the
 pipeline's workflow itself. In particular for binning 2x2 and custom ROI those
 parameters contained in _dcr.par_ has to be specifically tuned.

## Documentation

You will find a user manual on [goodman.readthedocs.org](http://goodman.readthedocs.io/en/latest/)

If you wish to know more about the instrument please check the 
[SOAR website](http://www.ctio.noao.edu/soar/content/goodman-high-throughput-spectrograph)

## Having trouble?

If you are having trouble operating the Goodman Pipeline we suggest the following
procedure.

* Check [existing issues](https://github.com/soar-telescope/goodman_pipeline/issues) or 
open a [new Issue](https://github.com/soar-telescope/goodman_pipeline/issues/new) on GitHub.

## Development Team

- [Simón Torres](https://github.com/simontorres) (SOAR Telescope Data Analyst - main code developer)
- [César Briceño](https://github.com/cbaorion) (SOAR Telescope Scientist - team lead)
- [Bruno Quint](https://github.com/b1quint) (Brazil Support Astronomer - code development adviser)


## Acknowledgements

We acknowledge the important contribution of  [David Sanmartim](https://github.com/dsanmartim), who developed
the initial incarnation of the redccd module. We thank [Tina Armond](https://github.com/tarmond) for her
invaluable help in adding calibrated comparison lamps to the library of
reference comparison lamps for wavelngth solution.

Our work would not be possible without the friendly work atmosphere at CTIO
headquarters in La Serena, were we can interact with our SOAR and CTIO
colleagues in lively and useful discussions that have been important in making
the Goodman pipeline possible.  We also acknowledge fruitful discussions and
suggestions from our colleagues Bart Dunlop, Chris Clemens, and Erik Dennihy,
at University of North Carolina at Chapel Hill.
  
## Citations:
  This pipeline makes extensive use of Astropy therefore you should cite as suggested
  on [Astropy Citation Page](https://github.com/astropy/astropy/blob/master/CITATION) as follows:
  
    This research made use of Astropy, a community-developed core Python package
    for Astronomy (Astropy Collaboration, 2013, 2018).
    
  It also uses [DCR](http://users.camk.edu.pl/pych/DCR/) for cosmic rays identification
  and removal. You should cite [this paper](http://adsabs.harvard.edu/abs/2004PASP..116..148P)
  
     Pych, W., 2004, PASP, 116, 148

