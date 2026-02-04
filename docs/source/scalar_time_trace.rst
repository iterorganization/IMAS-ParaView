.. _`using the Scalar Time Trace Reader`:

Scalar Time Trace Reader
========================

This page explains how to use the Scalar Time Trace Reader to visualize 0D data changing over time.


Supported IDSs
--------------

All IDSs are supported in the Scalar Time Trace Reader. Upon loading an IDS, 
the reader will automatically scan the IDS for quantities which are: 

- Time-dependent 
- Either floating point or integer
- Either a 0D quantity inside a time dependent Array of Structure, or a 1D quantity with time as its coordinate.

Using the Scalar Time Trace Reader
----------------------------------

The Scalar Time Trace Reader functions similarly to the GGD Reader, 
with the same interface and data loading workflow. 
This means that the steps for loading an URI, an IDS, and selecting attributes are identical. 
Refer to the :ref:`using the GGD Reader` for detailed instructions on:

- :ref:`Loading an URI <loading-an-uri>`: How to provide the file path or select a dataset.
- :ref:`Loading an IDS <loading-an-ids>`: How to load a dataset and display the grid.
- :ref:`Selecting attribute arrays <selecting-ggd-arrays>`: How to choose and visualize attributes.


Visualizing a time trace
------------------------

The Scalar Time Trace Reader outputs vtkTable data, containing the time traces of the
selected quantities. These quantities can be visualized in a 1d line chart plot.
After loading the the selected quantities using the steps above, you can visualize the data using
the 'Line Chart View' in ParaView. See `Paraview's documentation <https://docs.paraview.org/en/latest/UsersGuide/displayingData.html#line-chart-view>`_ 
to learn more about how to apply this view.

Within the Line Chart View options, disable the ``Use Index For X Axis`` option, 
and in the drop-down menu, select the ``Time [s]`` array to use it as X-axis in the line chart.

.. tip:: The Line Chart View can be opened side-by-side with other ParaView views. 
   So this reader will allow you to visualize how 0D quantities (e.g. the plasma current) 
   evolve over time, while visualizing other time-dependent 3D geometries. 

This reader supports two ways of visualizing data:

#. Partial Time Trace Up to Selected Time (Default)

   When the ``Show Full Time Trace`` option is disabled in the plugin properties, 
   the reader outputs only the time steps up to the currently selected time. 
   This provides the possibility to inspect the evolution of the selected quantities.
#. Full Time Trace with Highlighted Current Time

   When the ``Show Full Time Trace`` option is enabled, the reader outputs the complete 
   time series for all selected quantities. In addition, a ``Time Marker`` column is added, 
   which marks the currently selected time step. This allows the user to see the entire 
   time series from start to finish.

.. list-table::
   :widths: 50 49
   :header-rows: 0

   * - .. figure:: images/full_time_trace_disabled.png
     - .. figure:: images/full_time_trace_enabled.png
   * - Plasma current over time with ``Show Full Time Trace`` disabled.
     - Plasma current over time with ``Show Full Time Trace`` enabled.
