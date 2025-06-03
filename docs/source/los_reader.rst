.. _`using the Line of Sight Reader`:

Line of Sight Reader
====================

This page explains how to use the Line of Sight Reader to visualize Line of sight IDS data structures.


Supported IDSs
--------------

Currently, the following IDS and structures are supported in the Line of Sight Reader:

.. list-table::
   :widths: auto
   :header-rows: 1

   * - IDS
     - Structure
   * - ``bolometer``
     - `Channels <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/bolometer.html#bolometer-channel-line_of_sight>`__
   * - ``bremsstrahlung_visible``
     - `Channels <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/bremsstrahlung_visible.html#bremsstrahlung_visible-channel-line_of_sight>`__
   * - ``ece``
     - `Channels <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/ece.html#ece-line_of_sight>`__,
       `All Channels <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/ece.html#ece-line_of_sight>`__
   * - ``hard_x_rays``
     - `Channels <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/hard_x_rays.html#hard_x_rays-channel-line_of_sight>`__
   * - ``interferometer``
     - `Channels <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/interferometer.html#interferometer-channel-line_of_sight>`__
   * - ``mse``
     - `Channels <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/mse.html#mse-channel-line_of_sight>`__
   * - ``polarimeter``
     - `Channels <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/polarimeter.html#polarimeter-channel-line_of_sight>`__
   * - ``refractometer``
     - `Channels <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/refractometer.html#refractometer-channel-line_of_sight>`__
   * - ``soft_x_rays``
     - `Channels <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/soft_x_rays.html#soft_x_rays-channel-line_of_sight>`__
   * - ``spectrometer_uv``
     - `Channels <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/spectrometer_uv.html#spectrometer_uv-channel-line_of_sight>`__
   * - ``spectrometer_visible``
     - `Channels <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/spectrometer_visible.html#spectrometer_visible-channel-line_of_sight>`__

Using the Line of Sight Reader
------------------------------

The Line of Sight Reader functions similarly to the GGD Reader, with the same interface and data loading workflow. 
This means that the steps for loading an URI, an IDS, and selecting attributes are identical. 
Refer to the :ref:`using the GGD Reader` for detailed instructions on:

- :ref:`Loading an URI <loading-an-uri>`: How to provide the file path or select a dataset.
- :ref:`Loading an IDS <loading-an-ids>`: How to load a dataset and display the grid.
- :ref:`Selecting attribute arrays <selecting-ggd-arrays>`: How to choose and visualize attributes.


Setting the Scaling Factor
--------------------------

The Line of Sight Reader allows you the change the length the of the line-of-sight, using the 
``Scaling Factor``. This factor allows you to scale the line-of-sight to be either longer or shorter 
than its original length. Additionally, the Scaling Factor can be set to a negative value, 
which will invert the direction of the line-of-sight.


.. list-table::
   :widths: 50 49
   :header-rows: 0

   * - .. figure:: images/los_1.png
         :alt: Line of Sight at Scaling Factor 1
     - .. figure:: images/los_1_5.png
         :alt: Line of Sight at Scaling Factor 1.5
   * - Line of sights of vacuum vessel cameras in bolometer IDS, with ``Scaling Factor`` at 1 (default).
     - Line of sights of vacuum vessel cameras in bolometer IDS, with ``Scaling Factor`` at 1.5.
