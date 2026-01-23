.. _`using the 0D Time-Dependent Data Reader`:

0D Time-Dependent Data Reader
=============================

This page explains how to use the 0D Time-Dependent Data Reader to visualize 0D data changing over time.


Supported IDSs
--------------

Currently, all IDSs are supported in the 1D Profiles Reader. The reader will automatically
scan the provided IDS for quantities which are the following: 

- Time-dependent 
- Either floating point or integer
- Either 0D quantity inside a ``time_slice`` Array of Structure, or a 1D quantity with time as its coordinate.

Using the Position Reader
-------------------------

The Position Reader functions similarly to the GGD Reader, with the same interface and data loading workflow. 
This means that the steps for loading an URI, an IDS, and selecting attributes are identical. 
Refer to the :ref:`using the GGD Reader` for detailed instructions on:

- :ref:`Loading an URI <loading-an-uri>`: How to provide the file path or select a dataset.
- :ref:`Loading an IDS <loading-an-ids>`: How to load a dataset and display the grid.
- :ref:`Selecting attribute arrays <selecting-ggd-arrays>`: How to choose and visualize attributes.


Visualize the 0D time-dependent data
------------------------------------

The 0d Time-Dependent Data Reader outputs vtkTable data, which can be plotted in a 1D plot.
After loading the attribute arrays using the steps above, you can visualize the data using
the Line Chart View in ParaView. See `Paraview's documentation <https://docs.paraview.org/en/latest/UsersGuide/displayingData.html#line-chart-view>`_ 
to learn more about how to apply this and other views.

Within the Line Chart View options, disable the ``Use Index For X Axis`` option, 
and in the drop-down menu, select the ``Time [s]`` array to use as X-axis.
The data that should be plotted can now be selected.

.. tip:: The Line Chart View can be opened side-by-side with other ParaView views. 
   So this reader will allow you to visualize how 0D quantities (e.g. the plasma current) 
   evolve over time, while visualizing other time-dependent geometries (see for example
   the :ref:`visualizations of JOREK data <training_jorek>` in the training material.)
