.. _`using the Non-Axisymmetric Coils Reader`:

Non-Axisymmetric Coils Reader
=============================

This page explains how to use the Non-Axisymmetric Coils Reader to visualize the 
Non-axisymmetric active coil systems.

Supported IDSs
--------------

The following IDS and structures are supported in the Non-Axisymmetric Coils Reader:

.. list-table::
   :widths: auto
   :header-rows: 1

   * - IDS
     - Structure
   * - ``coils_non_axisymmetric``
     - `conductor <https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/coils_non_axisymmetric.html#coils_non_axisymmetric-coil-conductor>`__

Using the Non-Axisymmetric Coils Reader
---------------------------------------

The Non-Axisymmetric Coils Reader functions similarly to the GGD Reader, with the same interface and data loading workflow. 
This means that the steps for loading an URI, an IDS, and selecting attributes are identical. 
Refer to the :ref:`using the GGD Reader` for detailed instructions on:

- :ref:`Loading an URI <loading-an-uri>`: How to provide the file path or select a dataset.
- :ref:`Loading an IDS <loading-an-ids>`: How to load a dataset and display the grid.
- :ref:`Selecting attribute arrays <selecting-ggd-arrays>`: How to choose and visualize attributes.

Setting the Geometry Resolution
-------------------------------

The Non-Axisymmetric Coils Reader allows you to change two resolution parameters that 
control how coil geometry is discretized for visualization.

Resolution
^^^^^^^^^^

Defines the number of interpolation points used to represent curved conductor elements (arc_of_circle and full circle). 
Increasing this value produces smoother coil elements at the cost of higher geometric complexity.


.. figure:: images/coils_non_axi_res.png

   Arc of circle element (top) and a full circle element (bottom) with a resolution of 5 (left) and 50 (right)

Cross-sectional Resolution
^^^^^^^^^^^^^^^^^^^^^^^^^^

Defines the number of interpolation points used to approximate the conductor cross-section when a 
circular cross-sectional geometry is present. Increasing this value produces a smoother, 
more circular cross-section.

Example Case
------------

The following figure uses the Non-Axisymmetric Coils Reader to visualize the 
geometries of the ``coils_non_axisymmetric`` IDSs of the following machine description URIs:

- ``imas:hdf5?path=/work/imas/shared/imasdb/ITER_MD/3/111003/2``
- ``imas:hdf5?path=/work/imas/shared/imasdb/ITER_MD/3/115001/2``
- ``imas:hdf5?path=/work/imas/shared/imasdb/ITER_MD/3/115002/2``
- ``imas:hdf5?path=/work/imas/shared/imasdb/ITER_MD/3/115003/2``

.. figure:: images/non_axisymmetric.png
   :align: center

