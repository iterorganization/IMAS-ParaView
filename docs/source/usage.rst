.. _`usage`:

Using GGD-VTK
=============
TODO


Visualize GGD in Paraview
-------------------------
TODO

Convert GGD to VTK using the CLI
--------------------------------
GGD-VTK supports converting the GGD of an IDS to VTK format from a command-line interface (CLI).
This can be by using the `ggd2vtk` tool. Detailed usage, as well as examples can be found by running

  .. code-block:: bash

    vtkggdtools ggd2vtk --help

Example usage:

  .. code-block:: bash

    vtkggdtools ggd2vtk imas:hdf5?path=testdb#edge_profiles output_dir

