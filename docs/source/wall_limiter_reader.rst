.. _`using the Wall Limiter Reader`:

Wall Limiter Reader
===================

This page explains how to use the Wall Limiter Reader to visualize 
`limiter structures of description 2D structures <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/wall.html#wall-description_2d-limiter>`_ in a Wall IDS.
The `outline <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/wall.html#wall-description_2d-limiter-unit-outline>`_ 
of the limiter units is loaded into Paraview. 

.. note:: The ``phi_extensions`` property of the limiters is ignored, therefore this plugin is only applicable for axisymmetric limiter units.


Supported IDSs
--------------

Currently, the following IDS and structures are supported in the Wall Limiter Reader:

.. list-table::
   :widths: auto
   :header-rows: 1

   * - IDS
     - Structure
   * - ``wall``
     - `Limiter units <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/wall.html#wall-description_2d-limiter-unit>`__

Using the Wall Limiter Reader
-----------------------------

The Wall Limiter Reader functions similarly to the GGD Reader, with the same interface and data loading workflow. 
This means that the steps for loading an URI, an IDS, and selecting attributes are identical. 
Refer to the :ref:`using the GGD Reader` for detailed instructions on:

- :ref:`Loading an URI <loading-an-uri>`: How to provide the file path or select a dataset.
- :ref:`Loading an IDS <loading-an-ids>`: How to load a dataset and display the grid.
- :ref:`Selecting attribute arrays <selecting-ggd-arrays>`: How to choose and visualize attributes.


The following image shows the limiter structures of a description 2D structure of a wall IDS.

.. figure:: images/wall_limiter.png

   The first wall and divertor of a wall IDS.
