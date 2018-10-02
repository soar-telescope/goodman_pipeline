.. _questions-and-answers:

Questions & Answers
###################

1. What is the Goodman High Throughput Spectrograph?.

  This is thoroughly documented in `SOAR web site <http://www.ctio.noao.edu/soar/content/goodman-high-throughput-spectrograph>`_
  and links within.


2. How does the pipeline select the reference lamp?.

   The lamps are selected comparing two keywords. ``OBJECT`` and ``WAVMODE``

3. How should I organize the data?.

  More than organized your data should be *cleaned* of undesired files. There
  are some general assumptions in the implementation of the pipeline's data
  organization system that might get confused by files that are not supposed to
  be there.

4. What is *slit trim*?.
  Is a process to trim the 2D spectroscopic images to the
  *slit illuminated area* only. It works by fitting a box function to the
  dispersion-axis-collapsed spatial profile.

  The box function is :class:`~astropy.modeling.functional_models.Box1D` .
  The reason for doing it is because the non-illuminated area causes all sorts of
  problems in later steps, such as: existence of ``nan`` in master flats.