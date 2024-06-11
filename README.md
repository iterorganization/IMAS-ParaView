# VTKGGDTools

This project contains a ParaView plugin that exposes GGD readers and
GGD Writers for use in the ParaView pipeline.

## How to use

After installation, export the installation path of
vtkggdtools/plugins to PV_PLUGIN_PATH and launch ParaView. You can
also follow the "Developer instructions" below to have this
automatically done for you.

### View GGD (Grid and Plasma State) in ParaView

In ParaView, go to Sources->VTKGGDTools and choose the input plug-in
you want to use. Then fill the fields in the "Properties" panel with
the correct details of the IDS and press "Apply" (if the panel is not
visible, select it in the "View" menu). After a few seconds you should
see the grid-ggd colored by type of grid. You can now press the
"vtkBlockColors" button under the "Coloring" section and select the
plasma state quantity you want to view.

#### Example of IDSs to try:

You might find the list below useful to test. The reference order is
Pulse/Run/Tokamak/User.

- **134174/117/ITER/public** 
  
  Use TimeIdx 4
  - edge_profiles OK
  - edge_sources OK
  - edge_transport OK

- **135011/7/ggdtest/panchuj**

  TimeIdx 399, use for equilibrium
  - equilibrium OK

- **123217/1/ITER/public**
  - edge_profiles OK
  - edge_sources OK
  - edge_transport OK
  - radiation OK

- **122408/3/ITER/public**
  - edge_profiles OK
  - edge_sources OK
  - edge_transport OK
  - radiation OK

- **116100/1001/smiter/kosl**

  `~simicg/MESH_DIRECTORY/FULL_BLANKET_MESH/FullTokamak.med` (1.7 M triangles)
  - wall (geometry only) MDSPLUS - first wall - low resolution

- **116100/2001/smiter/kosl**

  `~simicg/MESH_DIRECTORY/DIVERTOR/Divertor.med` (0.6 M triangles) 
  - wall (geometry only) MDSPLUS - divertor

- **116100/3001/smiter/kosl**

  `~kosl/public/mesh/FirstWall.med` (9.1 M triangles)
  - wall (geometry only) MDSPLUS (0.3 GB) - first wall - high resolution
  It takes about 20 min for Paraview to open this IDS.

- **135913/4/smiter/kosl**

  Clipped right side of target from `~kosl/public/mesh/Powcal_powx_Temperatures.vtk`
  - wall target (geometry, temperature and power density) MDSPLUS

- **135913/5/smiter/kosl**

  `~kosl/public/mesh/Powcal_powx.vtk`
  - wall target (geometry and power density) MDSPLUS

- **135913/6/smiter/kosl**

  `~kosl/public/mesh/fingerleft8.vtk.vtk`
  - wall finger from target (geometry , temperature and power density) MDSPLUS

- **1/10/jorek/kosl**

  Bezier FEM ggd from jorek.
  - MHD

- **111114/1/test/artolaj**

  Bezier FEM ggd from jorek (K-Star).
  - MHD
  - radiation


### GGD (Grid and Plasma State) to VTK

Incomplete for transport_solver_numerics IDS.

Wall IDS grid (nodes, edges, faces, volumes) is read into one
ParaView's "partition dataset". This enables us to show values that are written
on cells (one value per cell) or show directly interpolated values from points
(to cells). If nodes and cells are not in same partition, just values specified
on cells can be shown on cells. Other subsets have their own "partition dataset". 
This "feature" is not implemented for other IDSs.

For MHD and radiation IDS reading 2D and 3D grids is enabled. For reading 3D
cases chosen number of planes must be more than zero.

### VTK to GGD (Grid and Plasma State)

Works for wall and (synthetic) equilibrium IDS.

WallGgdWriter writes temperature and power density from 'Temperature' and 'Q'
fields only (on nodes or cells). Names must be the same as mentioned in other
case just grid is written to IDS.

When writing grid for bigger meshes not including nodes, edges, faces and
volumes to subsets should be considered. This reduces IDS writting time and
it's size plus there is 'no need to replicate the grid elements in the
grid_subset structure' as said in data dictionary.

To be fully implemented.

## Tasks

- [x] Write plugin code for paraview. This will allow to read GGD from ParaView (VTKGGDTools.py)
- [x] Convert the geometric grid from `grids_ggd/*` to VTK
  following [ReadUALEdge/readGeometry.cxx](https://git.iter.org/projects/BND/repos/solps-gui/browse/src/plugins/paraview/readGmtryEdge.cxx)
  with some deviation. (read_geom.py)
- [x] Convert the geometric grid from VTK to `grids_ggd/*` (write_geom.py)
- [ ] Convert the plasma state description in `ggd/*` to VTK
  following [ReadUALEdge/readPsEdge.cxx](https://git.iter.org/projects/BND/repos/solps-gui/browse/src/plugins/paraview/readPsEdge.cxx) (
  read_ps.py)
- [ ] Convert the plasma state description from VTK to `ggd/*` (write_ps.py)

|IDS|  grid_ggd->VTK| VTK->grid_ggd| ggd->VTK| VTK->ggd|
|---|---|---|---|---|
|distribution_sources|done|done|done|n/a
|distributions | done |done|done|n/a
|equilibrium| done |done|done|done
|edge_profiles| done  |done|done|n/a
|edge_sources| done  |done|done|n/a
|edge_transport| done  |done|done|n/a
|mhd| done  |done|done|n/a
|radiation| done  |done|done|n/a
|tf| done  |done|done|n/a
|transport_solver_numerics| done  |done|incomplete|n/a
|wall| done  |done|done|done
|waves| done  |done|done|n/a

*Note*: The ReadUALEdge plugin uses `vtkMultiBlockDataSet`: a
deprecated VTK data structure and is limited to the `edge_profiles`
and a few edge related IDSs. This project uses
`vtkPartitionedDataSetCollection` and hopes to cover base all IDSs.

## Developer instructions

Load the pre-requisite modules, create and activate a virtual environment and install
the project files.

```bash
  # Load compatible IMAS and ParaView modules, like:
  # AL5 and ParaView 5.12 (recommended on RHEL9):
  $ module load IMAS-AL-Python/5.2.1-intel-2023b-DD-3.41.0 ParaView/5.12.0-foss-2023b
  # or AL4 and ParaView 5.10 (compatibility):
  # $ module load IMAS/3.41.0-4.11.9-foss-2020b ParaView/5.10.0-foss-2020b-mpi
  $ python -m venv --system-site-packages --clear --prompt vtkggddev .venv
  $ source .venv/bin/activate
  $ source install.sh .venv
  # Either launch paraview and test the plugins
  $ paraview
  # Or open up your IDE/code editor and begin development.
```
Upon modifying the source, run `source install.sh .venv` again.
