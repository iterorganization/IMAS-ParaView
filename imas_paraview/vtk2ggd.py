import logging

import imas
import numpy as np
import vtk
from imas import identifiers
from imas.ids_structure import IDSStructure
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonDataModel import (
    vtkPartitionedDataSetCollection,
    vtkUnstructuredGrid,
)

from imas_paraview.util import create_first_grid, int32array, vtk_cells_to_nodes

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class VTK2GGDConverter:
    """Handles the conversion of multiple vtkPartitionedDataSetCollections into a single
    IMAS GGD structure with multiple time steps.
    """

    # TODO: docstrings, extra unit tests, cleanup of grid subset logic

    def __init__(
        self,
        vtk_objects: list[vtkPartitionedDataSetCollection],
        ids_name,
        dd_version=None,
    ):
        self.ids_name = ids_name
        self.factory = imas.IDSFactory(version=dd_version)
        self.dd_version = self.factory.version
        self.vtk_objects = vtk_objects

    def convert(self):
        ids, grid_ggd = self._setup_ids()

        # Each partitioned dataset collection in the list contains the grid for a
        # single time step
        for time_idx, vtk_pdsc in enumerate(self.vtk_objects):
            # TODO: Only support linear grids
            grid_ggd[time_idx].identifier = identifiers.ggd_identifier.linear

            self._set_time_step(ids, time_idx, vtk_pdsc)
            space = self._setup_space(grid_ggd[time_idx], vtk_pdsc)

            self._fill_space(space, vtk_pdsc)
            self._build_grid_subsets(grid_ggd[time_idx], vtk_pdsc)

        return ids

    def _setup_ids(self):
        ids = self.factory.new(self.ids_name)
        first_grid = create_first_grid(ids)
        if not first_grid:
            raise RuntimeError(f"IDS '{self.ids_name}' does not have a GGD grid.")
        grid_ggd = imas.util.get_parent(first_grid)

        num_steps = len(self.vtk_objects)
        ids.time = np.zeros(num_steps)
        grid_ggd.resize(num_steps)

        # TODO: only support HOMOGENEOUS TIME for now
        ids.ids_properties.homogeneous_time = imas.ids_defs.IDS_TIME_MODE_HOMOGENEOUS

        # Set version_put properties (version_put was added in DD 3.22)
        if hasattr(ids.ids_properties, "version_put"):
            version_put = ids.ids_properties.version_put
            version_put.data_dictionary = self.dd_version
            version_put.access_layer_language = f"IMAS-Python {imas.__version__}"
        return ids, grid_ggd

    def _set_time_step(self, ids, time_idx, vtk_pdsc):
        info = vtk_pdsc.GetInformation()
        if info.Has(vtk.vtkDataObject.DATA_TIME_STEP()):
            time = info.Get(vtk.vtkDataObject.DATA_TIME_STEP())
            ids.time[time_idx] = time
            logger.info("Converting time step %d: t = %f", time_idx, time)
        else:
            logger.warning(
                "The vtkPartitionedDataSetCollection at time index %d does not have a "
                "time value. Using index as the time value instead.",
                time_idx,
            )
            ids.time[time_idx] = time_idx

    def read_vtk(self, vtk_pdsc):
        points = None
        objects_by_dimension = {}

        for i in range(vtk_pdsc.GetNumberOfPartitionedDataSets()):
            partitioned_dataset = vtk_pdsc.GetPartitionedDataSet(i)
            # GGD2VTK only fills first partition of a partitioned dataset
            ugrid = partitioned_dataset.GetPartition(0)
            if ugrid is None or not isinstance(ugrid, vtkUnstructuredGrid):
                logger.warning("partition %d does not contain an unstructured grid", i)
                continue
            if points is None and ugrid.GetNumberOfPoints() > 0:
                points = vtk_to_numpy(ugrid.GetPoints().GetData()).astype(np.float64)
            if ugrid.GetNumberOfCells() > 0:
                dim = ugrid.GetCell(0).GetCellDimension()
                if dim > 0:
                    objects_by_dimension.setdefault(dim, []).extend(
                        vtk_cells_to_nodes(ugrid.GetCells())
                    )
        return points, objects_by_dimension

    def _fill_space(self, space, vtk_pdsc):
        points, objects_by_dimension = self.read_vtk(vtk_pdsc)
        # Fill nodes
        if points is not None:
            space.objects_per_dimension[0].object.resize(len(points))
            for i, p in enumerate(points):
                obj = space.objects_per_dimension[0].object[i]
                obj.geometry.resize(3)
                # Wall converts coordinates differently from other IDSs
                if self.ids_name == "wall":
                    obj.geometry[:] = [p[0], p[1], p[2]]
                else:
                    obj.geometry[:] = [p[0], p[2], p[1]]

        # Fill cells
        for dim, cells in objects_by_dimension.items():
            space.objects_per_dimension[dim].object.resize(len(cells))
            for i, nodes in enumerate(cells):
                obj = space.objects_per_dimension[dim].object[i]
                obj.nodes = int32array(nodes)

    def _setup_space(self, grid_ggd, vtk_pdsc):
        # TODO: Only support creating a single space, I think multiple spaces are
        # lost upon conversion
        grid_ggd.space.resize(1)
        space = grid_ggd.space[0]
        # NOTE: Only support primary standard spaces
        space.identifier = identifiers.ggd_space_identifier.primary_standard
        # NOTE: Only support non-fourier grids
        space.geometry_type.index = 0
        # coordinates_type changed from INT_1D to an AoS of identifiers in DD4.0.0
        space.coordinates_type.resize(3)
        # NOTE: Only support X,Y,Z coordinates
        coord = identifiers.coordinate_identifier
        for i, c in enumerate([coord.x, coord.y, coord.z]):
            if isinstance(space.coordinates_type, IDSStructure):
                space.coordinates_type[i] = c
            else:
                space.coordinates_type[i] = c.index
        max_dim = self._get_max_dimension(vtk_pdsc)
        space.objects_per_dimension.resize(max_dim + 1)
        return space

    def _get_max_dimension(self, vtk_pdsc):
        max_dimension = 0
        for i in range(vtk_pdsc.GetNumberOfPartitionedDataSets()):
            pds = vtk_pdsc.GetPartitionedDataSet(i)
            ug = pds.GetPartition(0) if pds.GetNumberOfPartitions() > 0 else None
            if ug and ug.GetNumberOfCells() > 0:
                max_dimension = max(max_dimension, ug.GetCell(0).GetCellDimension())
        return max_dimension

    def _build_grid_subsets(self, grid_ggd, vtk_pdsc):
        num_partitions = vtk_pdsc.GetNumberOfPartitionedDataSets()
        grid_ggd.grid_subset.resize(num_partitions)

        # TODO add pytests for grid_subsets

        # Use current_grid's space
        used_offsets = dict.fromkeys(
            range(len(grid_ggd.space[0].objects_per_dimension)), 0
        )

        for p in range(num_partitions):
            subset = grid_ggd.grid_subset[p]
            name = vtk_pdsc.GetMetaData(p).Get(vtk_pdsc.NAME())
            pds = vtk_pdsc.GetPartitionedDataSet(p)
            ug = pds.GetPartition(0) if pds.GetNumberOfPartitions() > 0 else None

            if ug is None:
                subset.element.resize(0)
                subset.dimension = 1
                continue

            n_cells = ug.GetNumberOfCells()
            if n_cells == 0 and ug.GetNumberOfPoints() > 0:
                cell_dim = 0
                n_elem = ug.GetNumberOfPoints()
            else:
                cell_dim = ug.GetCell(0).GetCellDimension() if n_cells > 0 else 0
                n_elem = n_cells

            subset.dimension = cell_dim + 1
            identifier_name = name.lower().replace("-", "_").replace(" ", "_")
            if identifier_name in [m.name for m in identifiers.ggd_subset_identifier]:
                subset.identifier = identifiers.ggd_subset_identifier[identifier_name]
            else:
                subset.identifier.name = identifier_name
                subset.identifier.index = 0
                subset.identifier.description = "Unknown subset identifier"

            subset.element.resize(n_elem)
            for i in range(n_elem):
                el = subset.element[i]
                el.object.resize(1)
                obj = el.object[0]
                obj.space = 1
                obj.dimension = subset.dimension
                if cell_dim == 0:
                    obj.index = i + 1
                else:
                    obj.index = used_offsets[cell_dim] + i + 1

            used_offsets[cell_dim] += n_elem
