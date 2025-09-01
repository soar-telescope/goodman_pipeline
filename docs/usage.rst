.. _`usage`:

General Guidelines
##################

The |pipeline full name| is designed to be simple to use, however simple does
not always is the best case for everyone, thus |pipeline name| is also
flexible.

Getting Help
************

This manual is intended to be the preferred method to get help. However the quickest option is using ``-h`` or ``--help``

``redccd --help``

Will print the list of arguments along with a quick explanation and default values.

It is the same for ``redspec``

``redspec --help``

For astrometric processing:

``redastrometry --help``

And for photometric processing:

``redphotometry --help``

All commands will display their respective parameters, options, and usage examples with detailed explanations and default values.

.. include:: _observing.rst

.. include:: _running_prepare_data_for_reduction.rst


.. _`processing_data`:

Processing Data
###############

.. include:: _running_redccd.rst

.. include:: _running_redspec.rst

.. include:: _running_redastrometry.rst

.. include:: _running_redphotometry.rst

.. include:: _running_new_keywords.rst

.. include:: _cosmic_ray_removal.rst

.. include:: _flat_normalization.rst

.. include:: _extraction_methods.rst

.. include:: _file_prefixes.rst

.. include:: _file_suffixes.rst

.. include:: _shortcuts.rst

.. include:: _common_issues.rst
