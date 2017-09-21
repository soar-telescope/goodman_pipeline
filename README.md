# Goodman High Throughput Spectrograph Data Reduction Pipeline

The Goodman High Throughput Spectrograph (Goodman HTS) Data-Reduction Pipeline
is the SOAR Telescope's official data reduction pipeline for _Goodman HTS_.

It has been fully developed in Python 2.7 and mostly astropy affiliated packages
with the exception of [dcr](http://users.camk.edu.pl/pych/DCR/) an external tool
that does cosmic ray identification and correction. The reason for using it
instead of LACosmic was that it works very well for spectroscopic data and the
results are evidently superior. Some of the negative aspects of using this
external (meaning outside of Python domains) software were: The integration into
the pipeline's workflow and the use of an external `dcr.par` parameter file.
 Such parameters have to be changed by hand and can't be integrated into the
 pipeline's workflow itself. In particular for binning 2x2 and custom ROI those
 parameters contained in _dcr.par_ has to be specifically tuned.


The pipeline is currently in its first beta release version. We have thoroughly
tested it but we know there are still bugs to be found. If you found one or have
any comment please contact Simón Torres at storres [at] ctio noao edu.


The Goodman High Throughput Spectrograph is an imaging spectrograph, if you wish
to know more about the instrument please check the 
[SOAR website](http://www.ctio.noao.edu/soar/content/goodman-high-throughput-spectrograph)
 
To see full documentation please go to the GitHub hosted site for
[Goodman](https://soar-telescope.github.io/goodman/)

## What is contained in the project?

The Goodman HTS Pipeline has two packages, _goodman_ccd_ dedicated to the basic
ccd reduction process and _goodman_spec_ for advanced processing for 
spectroscopic data such as target identification and extraction for instance.

They are supposed to work in sequence and they are called from two indepentent 
scripts, _redccd_ and _redspec_. 

Other files included are:
- Readme.md (this)
- requirements.txt
- user_manual.pdf
- dcr binaries
- Library of Wavelength Calibrated lamps.

## How to install it?

In principle you don't need to install it since SOAR Telescope is providing
dedicated servers where you can run the pipeline.

If you are a _user_ you should follow the instructions in the `user_manual.pdf`
file, where you will find instruction for accessing the server and using the
pipeline, even for installing it.
  
The user manual file is included in the distribution.

If you want to try it out as a developer please request the developer's manual.


# Acknowledge

[David Sanmartim](https://github.com/dsanmartim) developed the first version of
the ccd reduction part before he moved to a new position. That first version was
included in the pipeline as a package and has evolved along all the code.