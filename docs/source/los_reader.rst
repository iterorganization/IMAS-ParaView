.. _`using the LoS Reader`:

LoS (Line of Sight) Reader
==========================

This page explains how to use the LoS (Line-of-Sight) Reader to visualize Line of sight IDS data structures.


Supported IDSs
--------------

Currently, the following IDSs are supported in the LoS Reader.

- ``bolometer``
- ``bremsstrahlung_visible``
- ``ece``
- ``hard_x_rays``
- ``interferometer``
- ``mse``
- ``polarimeter``
- ``refractometer``
- ``soft_x_rays``
- ``spectrometer_uv``
- ``spectrometer_visible``

Using the LoS Reader
--------------------

The LoS Reader functions similarly to the GGD Reader, with the same interface and data loading workflow. 
This means that the steps for loading an URI, an IDS, and selecting attributes are identical. 
Refer to the :ref:`using the GGD Reader` for detailed instructions on:

- :ref:`Loading an URI <loading-an-uri>`: How to provide the file path or select a dataset.
- :ref:`Loading an IDS <loading-an-ids>`: How to load a dataset and display the grid.
- :ref:`Selecting attribute arrays <selecting-ggd-arrays>`: How to choose and visualize attributes.


Setting the Scaling Factor
--------------------------

The LoS allows you the change the length the of the line-of-sight, using the 
``Scaling Factor``. This factor allows you to scale the line-of-sight to be either longer or shorter 
than its original length. Additionally, the Scaling Factor can be set to a negative value, 
which will invert the direction of the line-of-sight.


.. list-table::
   :widths: 50 49
   :header-rows: 0

   * - .. figure:: images/los_1.png
         :alt: LoS at Scaling Factor 1
     - .. figure:: images/los_1_5.png
         :alt: LoS at Scaling Factor 1.5
   * - Line of sights of vacuum vessel cameras in bolometer IDS, with ``Scaling Factor`` at 1 (default).
     - Line of sights of vacuum vessel cameras in bolometer IDS, with ``Scaling Factor`` at 1.5.
