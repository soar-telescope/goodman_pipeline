.. _questions-and-answers:

Questions & Answers
###################

1. What is the Goodman High Throughput Spectrograph?.

  answer.

2. How does the pipeline select the reference lamp?.

  answer.

3. How should I organize the data?.

  answer.

4. What is *slit trim*?.
  Is a process to trim the 2D spectroscopic images to the
  *slit illuminated area* only. It works by fitting a box function to the
  dispersion-axis-collapsed spatial profile.
  The box function is ``astropy.modeling.models.Box1D``.

  The reason for doing it is because the non-illuminated area causes all sorts of
  problems in later steps, such as: existance of ``nan`` in master flats.