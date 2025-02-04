.. _`using the Position Reader`:

Position Reader
===============

This page explains how to use the Position Reader to visualize position IDS data structures.


Supported IDSs
--------------

Currently, the following IDS and structures are supported in the 1D Profiles Reader:

.. list-table::
   :widths: auto
   :header-rows: 1

   * - IDS
     - Structure
   * - ``barometry``
     - `Gauges <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/barometry.html#barometry-gauge-position>`__
   * - ``langmuir_probes``
     - `Embedded probes <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/langmuir_probes.html#langmuir_probes-embedded-position>`__
   * - ``magnetics``
     - `Flux loops <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/magnetics.html#magnetics-flux_loop-position>`__,
       `Rogowski coils <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/magnetics.html#magnetics-rogowski_coil-position>`__,
       `Poloidal field probes <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/magnetics.html#magnetics-b_field_pol_probe-position>`__,
       `Toroidal field probes <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/magnetics.html#magnetics-b_field_phi_probe-position>`__

Using the Position Reader
-------------------------

The Position Reader functions similarly to the GGD Reader, with the same interface and data loading workflow. 
This means that the steps for loading an URI, an IDS, and selecting attributes are identical. 
Refer to the :ref:`using the GGD Reader` for detailed instructions on:

- :ref:`Loading an URI <loading-an-uri>`: How to provide the file path or select a dataset.
- :ref:`Loading an IDS <loading-an-ids>`: How to load a dataset and display the grid.
- :ref:`Selecting attribute arrays <selecting-ggd-arrays>`: How to choose and visualize attributes.


Visualize the positions
-----------------------

To visualize the selected positions, under ``Representation`` select ``Point Gaussian``. 
The ``Gaussian Radius`` can be changed to a suitable size such that the points are 
clearly visible.

.. figure:: images/magnetics_positions.png

   Positions of sensors in a magnetics IDS.
