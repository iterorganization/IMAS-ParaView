.. _`cli`:

Using the CLI
=============
This page go over how to use the CLI interface to convert a GGD to VTK.

Convert GGD to VTK using the CLI
--------------------------------
IMAS-ParaView supports converting the GGD of an IDS to VTK format from a command-line interface (CLI).
This can be done by using the `ggd2vtk` tool. 

.. tip:: Detailed usage of the CLI tool can be found by running ``ggd2vtk --help``

Example usage:

.. code-block:: bash

    ggd2vtk imas:hdf5?path=testdb#edge_profiles output_dir

This will create the following tree file structure:

.. code-block:: bash 

    output_dir/
    ├── edge_profiles_0/
    │   ├── edge_profiles_0_0_0.vtu
    │   ├── edge_profiles_0_1_0.vtu
    │   ├── edge_profiles_0_2_0.vtu
    │   …
    └── edge_profiles_0.vtpc

Here, the 0 in ``edge_profiles_0`` denotes that is has converted the first time step. 
The vtpc file contains a vtkPartitionedCollection, which can directly be loaded into Paraview
for visualization. The vtu files contain a vtkUnstructuredGrid, and can also be loaded into
Paraview directly.

Converting specific time steps
------------------------------

The CLI tools allows for conversion of only a single time step, a range of time steps or
all timesteps. The tool allows for index and time-based time slicing.

No time slicing
^^^^^^^^^^^^^^^

If no time slicing option are provided, only the first time step will be converted. To 
convert all time steps instead, use the ``-a`` or ``--all-times`` flag.

.. code-block:: bash

    $ ggd2vtk imas:hdf5?path=testdb#edge_profiles test_dir -a


.. note:: If you intend to convert a large chunk of the entire GGD dataset, it may be 
   beneficial to disable lazy loading. For this, the ``--no-lazy`` flag can be set. 
   By default, lazy loading is enabled, except when all time steps are converted, in 
   which case it is disabled. To enforce the use of lazy loading, use the ``--lazy`` flag.

Index-based slicing
^^^^^^^^^^^^^^^^^^^

To convert specific time step indices or a range of time step indices, the ``-i`` or 
``--index`` flag can be used.

To convert index 5:

.. code-block:: bash

    $ ggd2vtk imas:hdf5?path=testdb#edge_profiles test_dir -i 5

To convert indices 5, 8, and 9:

.. code-block:: bash

    $ ggd2vtk imas:hdf5?path=testdb#edge_profiles test_dir -i 5,8,9

To convert a range of indices, such as 2 upto 6 (including 2 and 6):

.. code-block:: bash

    $ ggd2vtk imas:hdf5?path=testdb#edge_profiles test_dir -i 2:6

Time-based slicing
^^^^^^^^^^^^^^^^^^

To convert specific time step or a range of time steps, the ``-t`` or 
``--time`` flag can be used. The time values will be interpreted in seconds.

To convert time step at 5.5s:

.. code-block:: bash

    $ ggd2vtk imas:hdf5?path=testdb#edge_profiles test_dir -t 5.5

To convert time steps 5.5s, 8s, and 9.1s:

.. code-block:: bash

    $ ggd2vtk imas:hdf5?path=testdb#edge_profiles test_dir -t 5.5,8,9.1

To convert all time steps that fall between 2.2s and 6.6s:

.. code-block:: bash

    $ ggd2vtk imas:hdf5?path=testdb#edge_profiles test_dir -t 2.2:6.6

.. note:: If the specified time step is not found in the IDS, the time step before the
   specified time step will be used instead.

