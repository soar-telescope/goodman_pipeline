---
layout: archive
title: "About"
permalink: /about/
author_profile: false
---

This is a Beta Release Candidate version of this pipeline,
although we have thoroughly tested it there might (and will) still be bugs.
Please let us know if you find any.
My e-mail is storres [at] ctio noao edu.

The Goodman High Throughput Spectrograph is an imaging spectrograph,
if you wish to know more about the instrument please check the
[SOAR website](http://www.ctio.noao.edu/soar/content/goodman-high-throughput-spectrograph)

To see full documentation please go to the GitHub hosted site for
[Goodman](https://soar-telescope.github.io/goodman/)

## Package contents?

This repository contains tools for processing Goodman's spectrograph data.
It is separated into two main components. The **CCD** part is done by _redccd_
originally developed by [David Sanmartim](https://github.com/dsanmartim)
and currently maintained by our team. It does
the standard _ccd image reduction_, i.e. trim, bias and flat correction.
Currently the ccd reduction pipeline has been integrated into this package so
there is no need to look for it separately. The **Spectroscopic** processing is
done by _redspec_ and includes the following features:


- [x] Identify targets in images
- [x] Trace the target
- [x] Extract the target with background subtraction
- [x] Find the wavelength Solution, interactive and automatically
- [x] Linearize data (resample)
- [x] Write wavelength solution to a FITS header
- [x] Create a new file for the wavelength calibrated 1D spectrum

There is also a library of calibrated spectrum in FITS format. Different
configurations are provided but is very easy
to add your own.

## How to get it?

  Check our [download page](/download/) to get more information about how to
  get the source files.

## Requirements

  Check our [requirement page](/requirements/) in order to learn how to make
  sure that you have all the required packages in your computer.

## How to install?

  The installation guidelines can be found in the [installation pages](/install/).
  Please, follow them and you will be soon ready to run your pipeline!

## Running the pipeline

  The [run page](/run/) gives you a brief description about how to run the pipeline.
  We will add more information and examples in a near future here.
