.. _`using the Magnetics Reader`:

Magnetics Reader
================

.. versionadded:: 2.4.0

This page explains how to use the Magnetics Reader to visualize magnetic diagnostic
geometries from the ``magnetics`` IDS.

Supported IDSs
--------------

The following structures are supported in the Magnetics Reader:

.. list-table::
   :widths: auto
   :header-rows: 1

   * - IDS
     - Structure
   * - ``magnetics``
     - `flux_loop <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/magnetics.html#magnetics-flux_loop>`__,
       `b_field_pol_probe <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/magnetics.html#magnetics-b_field_pol_probe>`__,
       `b_field_phi_probe <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/magnetics.html#magnetics-b_field_phi_probe>`__,
       `rogowski_coil <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/magnetics.html#magnetics-rogowski_coil>`__

Using the Magnetics Reader
--------------------------

The Magnetics Reader functions similarly to the GGD Reader, with the same interface and data loading workflow.
This means that the steps for loading an URI, an IDS, and selecting attributes are identical.
Refer to the :ref:`using the GGD Reader` for detailed instructions on:

- :ref:`Loading an URI <loading-an-uri>`: How to provide the file path or select a dataset.
- :ref:`Loading an IDS <loading-an-ids>`: How to load a dataset and display the grid.
- :ref:`Selecting attribute arrays <selecting-ggd-arrays>`: How to choose and visualize attributes.

Visualized Elements
-------------------

Flux Loops and Rogowski Coils
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Flux loops and Rogowski coils are visualized as closed polylines representing the wire geometry.
Each element is displayed using the position coordinates defined in the IDS.

B-Field Probes
^^^^^^^^^^^^^^

Poloidal and toroidal B-field probes are visualized as arrows indicating the sensor orientation.
The arrow direction is determined by the poloidal and toroidal angles defined for each probe.

The Probe Arrow Length parameter controls the length of the arrows drawn for each B-field probe
(in meters). The arrow tip is placed at the specified distance from the probe center along the
sensor normal axis.

Example Case
------------

The following figure uses the Magnetics Reader to visualize the geometries of the 
flux loops, rogowski coils, and B-field probes of ``magnetics`` from the following 
machine description NetCDF file: `iter_md_magnetics_150100_5.nc <https://zenodo.org/records/15525525>`_.

.. figure:: ../images/magnetics_reader.png

   Visualization of flux loops (blue), Rogowski coils (green), poloidal B-field probes (red),
   and toroidal B-field probes (white) from a machine description ``magnetics`` IDS.

