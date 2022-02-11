VTKGGDTools
===========

This project contains a ParaView plugin that exposes GGD readers and GGD Writers for use in the ParaView pipeline.

After installation, export the installation path of vtkggdtools/plugins to PV_PLUGIN_PATH.

Tasks
-----

- [x] Write plugin code for paraview. This will allow to read GGD from ParaView (VTKGGDTools.py)
- [x] Translate `grids_ggd/*` -> VTK
  following [ReadUALEdge/readGeometry.cxx](https://git.iter.org/projects/BND/repos/solps-gui/browse/src/plugins/paraview/readGmtryEdge.cxx)
  with some deviation. (read_geom.py)
- [x] Translate VTK -> `grids_ggd/*` (write_geom.py)
- [ ] Translate `ggd/*` -> VTK
  following [ReadUALEdge/readPsEdge.cxx](https://git.iter.org/projects/BND/repos/solps-gui/browse/src/plugins/paraview/readPsEdge.cxx) (
  read_ps.py)
- [ ] Translate VTK -> `ggd/*` (write_ps.py)

|IDS|  read grid_ggd| write grid_ggd| read plasma state| write plasma state|
|---|---|---|---|---|
|distribution_sources|done|done|done|n/a
|distributions | done |done|done|n/a
|equilibrium| done |done|done|done
|edge_profiles| done  |done|incomplete|n/a
|edge_sources| done  |done|incomplete|n/a
|edge_transport| done  |done|n/a|n/a
|mhd| done  |done|n/a|n/a
|radiation| done  |done|n/a|n/a
|tf| done  |done|done|n/a
|transport_solver_numerics| done  |done|incomplete|n/a
|wall| done  |done|done|n/a
|waves| done  |done|done|n/a

*Note*: The ReadUALEdge plugin uses `vtkMultiBlockDataSet`: a deprecated VTK data structure and is limited to
the `edge_profiles` and a few edge related IDSs. This project uses `vtkPartitionedDataSetCollection` and hopes to cover
base all IDSs.
