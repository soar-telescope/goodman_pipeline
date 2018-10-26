.. _extraction-methods:

Extraction Methods
******************

The argument ``--extraction <method>`` has two options but only ``fractional``
is implemented.

``fractional``:
  *Fractional pixel extraction* differs from a simple and rough extraction
  in how it deals with the edges of the region.
  :meth:`goodman_pipeline.core.core.extract_fractional_pixel`

``optimal``:
  Unfortunately this method has not been implemented yet.