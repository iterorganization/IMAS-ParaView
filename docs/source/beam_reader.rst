.. _`using the Beam Reader`:

Beam Reader
===========

This page explains how to use the Beam Reader to visualize Beam IDS data structures. 

.. note:: The EC beam directions are visualized as straight lines from their initial 
   launching position and angle. This is not the actual path of the EC waves through the plasma.

Supported IDSs
--------------

Currently, the following IDS and structures are supported in the Beam Reader:

     
.. list-table::
   :widths: auto
   :header-rows: 1

   * - IDS
     - Structure
   * - ``ec_launchers``
     - `Beam <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/ec_launchers.html#ec_launchers-beam>`__

Using the Beam Reader
---------------------

The Beam Reader functions similarly to the GGD Reader, with the same interface and data loading workflow. 
This means that the steps for loading an URI, an IDS, and selecting attributes are identical. 
Refer to the :ref:`using the GGD Reader` for detailed instructions on:

- :ref:`Loading an URI <loading-an-uri>`: How to provide the file path or select a dataset.
- :ref:`Loading an IDS <loading-an-ids>`: How to load a dataset and display the grid.
- :ref:`Selecting attribute arrays <selecting-ggd-arrays>`: How to choose and visualize attributes.



Setting the Beam Distance
-------------------------

The Beam Reader allows you to specify the length of the beams (in meters). This is the 
distance from the launching point of the beam, in the beam direction.

.. list-table::
   :widths: 50 50
   :header-rows: 0

   * - .. figure:: images/beam_2m.png

         beam structures in ec_launchers with ``Beam Distance`` set to 2 meters.
     - .. figure:: images/beam_10m.png

         beam structures in ec_launchers with ``Beam Distance`` set to 10 meters.
