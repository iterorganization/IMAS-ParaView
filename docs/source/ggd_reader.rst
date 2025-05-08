.. _`using the GGD Reader`:

GGD Reader
==========
This page will go over how to use the GGD Paraview plugin.

.. tip:: More information about the usage of specific plugin UI elements can be obtained by hovering
   the mouse over the element.

Supported IDSs
--------------

Currently, the following IDSs are supported in the GGD Reader.

- ``edge_profiles``
- ``edge_sources``
- ``edge_transport``
- ``mhd``
- ``radiation``
- ``runaway_electrons``
- ``wall``
- ``equilibrium``
- ``distribution_sources``
- ``distributions``
- ``tf``
- ``transport_solver_numerics``
- ``waves``
- ``plasma_profiles``,
- ``plasma_sources``,
- ``plasma_transport``,

.. _loading-an-uri:

Loading an URI
--------------
If Paraview is correctly installed according to the :ref:`installation docs <installing>`, you should
see the GGD reader plugin under `Sources > IMAS Tools > GGD Reader`. Selecting this source
will create a GGD Reader source in the Paraview pipeline.

The first step is to select which URI you want to load. The GGD Reader provides three different
methods of supplying an URI.

1. A string containing the URI.
2. A file browser with which you can select the file to load.
3. Manual selection of the backend, database, pulse, run, user, and version.

The different input options are shown in the figure below:

.. list-table::
   :widths: 33 33 33
   :header-rows: 0

   * - .. figure:: images/input_uri.png

         1\. String input.
     - .. figure:: images/input_file.png

         2\. File browser.
     - .. figure:: images/input_legacy.png

         3\. Manual input.

.. _loading-an-ids:

Loading an IDS
--------------
When the URI input is selected, and filled accordingly, you can press `Apply` to load the URI.
If the URI loads successfully, an IDS can be selected from the drop-down list. The reader automatically
detects which IDSs are available for the given URI, and only shows applicable IDSs for this plugin.

.. figure:: images/ids_dropdown.png

   The drop-down list to select an IDS.

.. |ico1| image:: images/rotate_axis.png

When the desired IDS is selected, click on `Apply` to load it. After the IDS is loaded, 
the grid should appear in the viewport. If the grid does not show up, it is
possible that the axes are not aligned properly. You can use the buttons above to align your viewpoint
with the data. For example, to align the viewpoint in the positive Y direction, press: |ico1|.

You should now be able to see the mesh of the GGD.

.. figure:: images/mesh.png

   The mesh of the GGD.

.. |ico2| image:: images/filter_icon.png

.. _selecting-ggd-arrays:

Selecting GGD arrays
--------------------
Attributes which are defined on the GGD grid can be loaded using the Attribute Selection window. 
The selection window can be found on the left under `Select attribute Arrays`.
Any number of attributes can be selected, and the list
can be filtered using the |ico2|-icon in the top right of the selection window. 
When the attributes are selected, click on `Apply` to load the selected attributes from the backend.
The selected quantities can now be visualized using Paraview's selection drop-down menu.

.. figure:: images/electron_pressure.png

   Visualize the `Electron Pressure` of the GGD.

.. note:: By default, the data is lazy loaded. This means that the data will only be fetched from
   the backend as soon it is required by the plugin. This helps loading times, as we don't need to
   fetch the entire dataset before handling any data. If you plan on loading the entire dataset, it
   may be beneficial to load all data up front. This can be done by enabling the `Preload Data`
   checkbox in the plugin settings.

