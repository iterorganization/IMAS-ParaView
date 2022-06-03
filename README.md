# VTKGGDTools

This project contains a ParaView plugin that exposes GGD readers and GGD Writers for use in the ParaView pipeline.

## How to use

After installation, export the installation path of vtkggdtools/plugins to PV_PLUGIN_PATH and launch ParaView. You can also follow the "Developer instructions" below to have this automatically done for you.

### View GGD (Grid and Plasma State) in ParaView

In ParaView, go to Sources->VTKGGDTools and choose the input plug-in you want to use. Then fill the fields in the "Properties" panel with the correct details of the IDS and press "Apply" (if the panel is not visible, select it in the "View" menu). After a few seconds you should see the grid-ggd colored by type of grid. You can now press the "vtkBlockColors" button under the "Coloring" section and select the plasma state quantity you want to view.

#### Example of IDSs to try:
You might find the list below useful to test. The reference order is Pulse/Run/Tokamak/User.
- **134174/117/ITER/public**<br>
  Use TimeIdx 4
  - edge_profiles OK
  - edge_sources OK
  - edge_transport OK
- **135011/7/ggdtest/panchuj**<br>
  TimeIdx 399, use for equilibrium
  - equilibrium OK
- **123217/1/ITER/public**<br>
  - edge_profiles OK
  - edge_sources OK
  - edge_transport OK
  - radiation OK
- **122408/3/ITER/public**<br>
  - edge_profiles OK
  - edge_sources OK
  - edge_transport OK
  - radiation OK

### Export GGD (Grid and Plasma State) to VTK

To be written.

### Export VTK to GGD (Grid and Plasma State)

To be written.

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
|wall| done  |done|incomplete|n/a
|waves| done  |done|done|n/a

*Note*: The ReadUALEdge plugin uses `vtkMultiBlockDataSet`: a deprecated VTK data structure and is limited to
the `edge_profiles` and a few edge related IDSs. This project uses `vtkPartitionedDataSetCollection` and hopes to cover
base all IDSs.

## Developer instructions

Load the pre-requisite modules, create and activate a virtual environment and install
the project files.
  ```bash
  $ module load IMAS ParaView
  $ python -m venv --system-site-packages --clear --prompt vtkggddev .venv
  $ source .venv/bin/activate
  $ source install.sh .venv
  # Either launch paraview and test the plugins
  $ paraview
  # Or open up your IDE/code editor and begin development.
  ```
Upon modifying the source, run `source install.sh .venv` again.
