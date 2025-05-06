.. _`using the 1d Profiles Reader`:

1D Profiles Reader
==================

This page explains how to use the 1D Profiles Reader to visualize 1D profile data.


Supported IDSs
--------------

Currently, the following IDS and structures are supported in the 1D Profiles Reader:

.. list-table::
   :widths: auto
   :header-rows: 1

   * - IDS
     - Structure
   * - ``core_profiles``
     - `Core plasma radial profiles <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/core_profiles.html#core_profiles-profiles_1d>`__
   * - ``core_sources``
     - `Source profiles <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/core_sources.html#core_sources-source-profiles_1d>`__

Using the 1D Profiles Reader
----------------------------

The 1D Profiles Reader functions similarly to the GGD Reader, with the same interface and data loading workflow. 
This means that the steps for loading an URI, an IDS, and selecting attributes are identical. 
Refer to the :ref:`using the GGD Reader` for detailed instructions on:

- :ref:`Loading an URI <loading-an-uri>`: How to provide the file path or select a dataset.
- :ref:`Loading an IDS <loading-an-ids>`: How to load a dataset and display the grid.
- :ref:`Selecting attribute arrays <selecting-ggd-arrays>`: How to choose and visualize attributes.

Mapping 1D profiles to 2D
-------------------------

For information on how to map 1D profiles to a 2D grid, refer to the section :ref:`using the 1d Profiles Mapper`.

Plotting 1D profiles
--------------------

The 1D Profiles Reader outputs vtkTable data, this can be plotted in a 1D plot.
After loading the attribute arrays using the steps above, apply the ``Plot Data`` Paraview filter 
(For more information on applying filters in Paraview, have a look at the `Paraview documentation <https://docs.paraview.org/en/latest/UsersGuide/filteringData.html>`_).
Within the ``Plot Data`` filter options, disable the ``Use Index For X Axis`` option, 
and in the drop-down menu, select the array to use as X-axis (e.g. ``rho_tor_norm``).
The profiles that should be plotted can now be selected.


.. list-table::
   :widths: 50 50
   :header-rows: 0

   * - .. figure:: images/profiles_plot_options.png

         ``rho_tor_norm`` is chosen for the X-axis and the ``Ion Density`` is chosen as array to plot
     - .. figure:: images/profiles_plot.png

         The plotted ``Ion Density`` as a function of ``rho_tor_norm``
