.. _observing:

Observing Guidelines
********************

In order to be able to process your data with the |pipeline full name| you need
to follow some guidelines, we do not intend to tell you how to do your science.
Here are some basic hints.

- Make sure you have a good observing plan as well as a good backup plan too.
- Put special attention to the calibration files that are needed for the data
  that you are planning to obtain, for instance, you can process your
  spectroscopic data without bias because using overscan will give you a good
  enough approximation, but Imaging does not have overscan therefore you MUST
  obtain bias frames.
- Keep a detailed log of things that happened while you were observing,
  mistakes that you made, exposures repeated, etc. An observing log is not an
  extraction of header information. Well, it can be, but it will be useless.
- If you are unsure about the required steps to achieve your science goals ask
  your PI, not the support scientist, Her/His job is to assist you on how to get
  good quality data not what data you need in order to achieve your scientific
  goals.

For using the pipeline you don't need to use any special file naming convention,
in fact all the information is obtained from the headers. As of version
``1.2.0`` you need to use a reference lamp naming convention though. Not the
file but the field that goes into OBJECT. It is actually very simple:

.. _table-lamp-names:

.. table:: Convention names for comparison lamps

    ======================= ===========================
     Lamp name               Convention
    ======================= ===========================
     Argon                   Ar
     Neon                    Ne
     Copper                  CuHeAr
     Iron                    FeHeAr
     Mercury Argon           HgAr
     Mercury Argon Neon      HgArNe
    ======================= ===========================

This is to ensure the pipeline is able to recognize them. This will no be the
case in future versions but for now this is how it works.

Observing for Radial Velocity
*****************************

Radial velocity measurements are possible with the |goodman HTS| but you have
to be careful. A very detailed description of the procedures and what you can
expect was prepared and is available  `here <https://noirlab.edu/science/sites/default/files/media/archives/documents/scidoc0489.pdf>`_
and `here <https://noirlab.edu/science/sites/default/files/media/archives/documents/scidoc0490.pdf>`_ .

Please read it carefully so you don't find any surprises when trying to reduce
your data.

