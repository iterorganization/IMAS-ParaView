.. _`using the PF Reader`:

PF Reader
=========

This page explains how to use the PF Reader to visualize the axisymmetric active poloidal 
field coils, as well as axisymmetric passive conductor structures.

Supported IDSs
--------------

The following IDS and structures are supported in the PF Reader:

.. list-table::
   :widths: auto
   :header-rows: 1

   * - IDS
     - Structure
   * - ``pf_active``
     - `coil geometries <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/pf_active.html#pf_active-coil-element-geometry>`__
   * - ``pf_passive``
     - `loop geometries <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/pf_passive.html#pf_passive-loop-element-geometry>`__

Using the PF Reader
-------------------

The PF Reader functions similarly to the GGD Reader, with the same interface and data loading workflow. 
This means that the steps for loading an URI, an IDS, and selecting attributes are identical. 
Refer to the :ref:`using the GGD Reader` for detailed instructions on:

- :ref:`Loading an URI <loading-an-uri>`: How to provide the file path or select a dataset.
- :ref:`Loading an IDS <loading-an-ids>`: How to load a dataset and display the grid.
- :ref:`Selecting attribute arrays <selecting-ggd-arrays>`: How to choose and visualize attributes.
