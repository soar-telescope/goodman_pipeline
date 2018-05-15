.. _dev-versioning:

Version Numbering: Proposal
###########################

Introduction
************

Version number can be a confusing matter, specially because there are many ways
to do it and at the end it is a matter of personal (or team choice).
After doing some digging in python PEPs and internet in general I found this
`Wikipedia site on Software Versioning <https://en.wikipedia.org/wiki/Software_versioning#Designating_development_stage>`_
that to me makes a lot of sense and what I propose here is the following:

Pipeline Versioning
*******************

Let's define our goal to be the version 1.0 as the one that will be stable,
portable and also must contain an automatic wavelength calibration module as
well as a flux calibration module. Besides, the list of issues tagged as bug,
or enhancements should be at least 70 ~ 90% closed. That is a large gap but some
enhancements might be of less importance, bugs will be naturally prioritized.

According to the `Wikipedia site <https://en.wikipedia.org/wiki/Software_versioning#Designating_development_stage>`_
and the current status of our development the number would be something like
this:

``1.0a1`` which can be read as **"we are developing towards version 1.0 and
currently we are working on the alpha version of our pipeline"**.

For the first beta release the version number will be ``1.0b1`` and as we add
improvements and fix bugs we can do for instance ``1.0b2`` which would be
considered a new beta release, or in other words *beta version 2*

Right now numbering is a mere formality and I have been experimenting in
introducing them since some time ago, but they will become very important
when we go public, for the first beta release.

Summary and Schedule
********************

- ``1.0a1`` until the end of May 2017.
- ``1.0b1`` From beginning of June until end of August (last digit may vary)
- ``1.0`` Sometime between September and before the end of the year (ideally October).