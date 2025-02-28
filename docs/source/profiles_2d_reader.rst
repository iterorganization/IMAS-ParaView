.. _`using the 2d Profiles Reader`:

2D Profiles Reader
==================

This page explains how to use the 2D Profiles Reader to visualize 2D profile data.


Supported IDSs
--------------

Currently, the following IDS and structures are supported in the 1D Profiles Reader:

.. list-table::
   :widths: auto
   :header-rows: 1

   * - IDS
     - Structure
   * - ``core_profiles``
     - `Core plasma poloidal cross sections <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/core_profiles.html#core_profiles-profiles_2d>`__
   * - ``equilibrium``
     - `Equilibrium 2D profiles <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/equilibrium.html#equilibrium-time_slice-profiles_2d>`__
   * - ``plasma_initiation``
     - `2D poloidal profiles <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/plasma_initiation.html#plasma_initiation-profiles_2d>`__
   * - ``plasma_profiles``
     - `2D poloidal profiles <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/plasma_profiles.html#plasma_profiles-profiles_2d>`__

Using the 2D Profiles Reader
----------------------------

The 2D Profiles Reader functions similarly to the GGD Reader, with the same interface and data loading workflow. 
This means that the steps for loading an URI, an IDS, and selecting attributes are identical. 
Refer to the :ref:`using the GGD Reader` for detailed instructions on:

- :ref:`Loading an URI <loading-an-uri>`: How to provide the file path or select a dataset.
- :ref:`Loading an IDS <loading-an-ids>`: How to load a dataset and display the grid.
- :ref:`Selecting attribute arrays <selecting-ggd-arrays>`: How to choose and visualize attributes.

Plotting 2D profiles
--------------------

The 2D profile reader will try to load the first valid ``profiles_2d`` node that is encountered. 
A ``profiles_2d`` node is valid if both the grid coordinates (``r``- and ``z``-coordinates), and at least one profile is found. 
If the `radial <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/equilibrium.html#equilibrium-time_slice-profiles_2d-r>`_ and `height <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/equilibrium.html#equilibrium-time_slice-profiles_2d-z>`_ coordinates are filled, these will be loaded.
Alternatively, if these coordinates are empty, it will try to load the radial and height coordinates from the `dim1 <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/equilibrium.html#equilibrium-time_slice-profiles_2d-grid-dim1>`_ and `dim2 <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/equilibrium.html#equilibrium-time_slice-profiles_2d-grid-dim2>`_ coordinates, respectively. Note, this is only supported for rectangular grids.

The VTK plugin will output a ``vtkMultiBlockDataSet`` where each block contains a single profile stored as a ``vtkUnstructuredGrid``. The following image shows a 2D profile from an equilibrium IDS.

.. figure:: images/profiles_2d.png

   Contour plot of a 2D profile of poloidal magnetic flux, with the divertor and first wall as reference.
