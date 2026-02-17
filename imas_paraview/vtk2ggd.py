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
    """Handles the conversion of multiple vtkPartitionedDataSetCollections, where
    each vtkPartitionedDataSetCollection contains a grid for a single time step."""

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
        """Converts a list of vtkPartitionedDataSetCollections into a single IDS. Each
        vtkPartitionedDataSetCollection must contain a single partition containing a
        vtkUnstructuredGrid.

        Returns:
            An IDS containing the converted GGD grid for all time steps.
        """
        ids, grid_ggd = self._setup_ids()

        # Each partitioned dataset collection in the list contains the grid for a
        # single time step
        for time_idx, vtk_pdsc in enumerate(self.vtk_objects):
            # NOTE: We only support linear grids
            grid_ggd[time_idx].identifier = identifiers.ggd_identifier.linear

            time = self._get_time_step(vtk_pdsc)
            if time is not None:
                logger.info("Converting time step %d: t = %f", time_idx, time)
            else:
                logger.warning(
                    "The vtkPartitionedDataSetCollection at time index %d does not "
                    "have a time value. Using index as the time value instead.",
                    time_idx,
                )
                time = time_idx
            ids.time[time_idx] = time

            space = self._setup_space(grid_ggd[time_idx], vtk_pdsc)

            self._fill_space(space, vtk_pdsc)
            self._build_grid_subsets(grid_ggd[time_idx], vtk_pdsc)

        return ids

    def _setup_ids(self):
        """Create and initialize the IDS and its GGD grid structure.

        Returns:
            ids: The initialized IDS.
            grid_ggd: The GGD grid AoS for the created IDS.
        """
        ids = self.factory.new(self.ids_name)
        first_grid = create_first_grid(ids)
        if not first_grid:
            raise RuntimeError(f"IDS '{self.ids_name}' does not have a GGD grid.")
        grid_ggd = imas.util.get_parent(first_grid)

        num_steps = len(self.vtk_objects)
        ids.time = np.zeros(num_steps)
        grid_ggd.resize(num_steps)

        # NOTE: We only support IDSs with a homogeneous time mode
        ids.ids_properties.homogeneous_time = imas.ids_defs.IDS_TIME_MODE_HOMOGENEOUS

        return ids, grid_ggd

    def _get_time_step(self, vtk_pdsc):
        """Return the time value of the vtkPartitionedDataSetCollection, or None if it
        does not contain any time information."""
        info = vtk_pdsc.GetInformation()
        if info.Has(vtk.vtkDataObject.DATA_TIME_STEP()):
            time = info.Get(vtk.vtkDataObject.DATA_TIME_STEP())
            return time
        return None

    def _setup_space(self, grid_ggd, vtk_pdsc):
        """Initialize and configure the GGD space definition for a grid.

        Args:
            grid_ggd: GGD grid structure for the current time step.
            vtk_pdsc: VTK partitioned dataset collection used for dimension inference.

        Returns:
            Initialized GGD space structure.
        """
        # NOTE: We only support creating a single GGD space
        grid_ggd.space.resize(1)
        space = grid_ggd.space[0]
        # NOTE: We only support primary standard spaces
        space.identifier = identifiers.ggd_space_identifier.primary_standard
        # NOTE: We only support non-fourier grids
        space.geometry_type.index = 0

        # NOTE: We only support X,Y,Z coordinates
        space.coordinates_type.resize(3)
        coord = identifiers.coordinate_identifier

        # GGD2VTK converts coordinates differently for wall than from other IDSs
        if self.ids_name == "wall":
            coord_identifiers = [coord.x, coord.z, coord.y]
        else:
            coord_identifiers = [coord.x, coord.y, coord.z]

        # coordinates_type changed from INT_1D to an AoS of identifiers in DD4.0.0
        for i, coord_identifier in enumerate(coord_identifiers):
            if isinstance(space.coordinates_type, IDSStructure):
                space.coordinates_type[i] = coord_identifier
            else:
                space.coordinates_type[i] = coord_identifier.index
        max_dim = self._get_max_dimension(vtk_pdsc)
        space.objects_per_dimension.resize(max_dim + 1)
        return space

    def _get_max_dimension(self, vtk_pdsc):
        """Returns the maximum cell dimension in the vtkPartitionedDataSetCollection."""
        max_dimension = 0
        for i in range(vtk_pdsc.GetNumberOfPartitionedDataSets()):
            pds = vtk_pdsc.GetPartitionedDataSet(i)
            ugrid = pds.GetPartition(0) if pds.GetNumberOfPartitions() > 0 else None
            if ugrid and ugrid.GetNumberOfCells() > 0:
                max_dimension = max(max_dimension, ugrid.GetCell(0).GetCellDimension())
        return max_dimension

    def _fill_space(self, space, vtk_pdsc):
        """Populate a GGD space with node coordinates and cell connectivities.

        Args:
            space: GGD space structure to populate.
            vtk_pdsc: vtkPartitionedDataSetCollection to read data from.
        """
        points, objects_per_dimension = self.extract_data_from_pdsc(vtk_pdsc)

        # Fill nodes
        if points is not None:
            space.objects_per_dimension[0].object.resize(len(points))
            for point, obj in zip(points, space.objects_per_dimension[0].object):
                # GGD2VTK converts coordinates differently for wall than from other IDSs
                if self.ids_name != "wall":
                    point = [point[0], point[2], point[1]]
                obj.geometry = point

        # Fill cells
        for dim, cells in objects_per_dimension.items():
            space.objects_per_dimension[dim].object.resize(len(cells))
            for i, nodes in enumerate(cells):
                obj = space.objects_per_dimension[dim].object[i]
                obj.nodes = int32array(nodes)

    def extract_data_from_pdsc(self, vtk_pdsc):
        """Extract point coordinates and cell connectivities grouped by dimension from
        a vtkPartitionedDataSetCollection.

        Args:
            vtk_pdsc: vtkPartitionedDataSetCollection to extract data from.

        Returns:
            points: Array containing the x,y,z-points
            objects_per_dimension: Dictionary mapping dimension to node lists.
        """
        points = None
        objects_per_dimension = {}

        for i in range(vtk_pdsc.GetNumberOfPartitionedDataSets()):
            partitioned_dataset = vtk_pdsc.GetPartitionedDataSet(i)
            # NOTE: We only read the first partition here, as GGD2VTK only fills first
            # partition of a partitioned dataset
            ugrid = partitioned_dataset.GetPartition(0)
            if ugrid is None or not isinstance(ugrid, vtkUnstructuredGrid):
                logger.warning("partition %d does not contain an unstructured grid", i)
                continue
            if points is None and ugrid.GetNumberOfPoints() > 0:
                points = vtk_to_numpy(ugrid.GetPoints().GetData()).astype(np.float64)
            if ugrid.GetNumberOfCells() > 0:
                dim = ugrid.GetCell(0).GetCellDimension()
                if dim > 0:
                    objects_per_dimension.setdefault(dim, []).extend(
                        vtk_cells_to_nodes(ugrid.GetCells())
                    )
        return points, objects_per_dimension

    def _build_grid_subsets(self, grid_ggd, vtk_pdsc):
        """Create grid subsets mapping partitions of a vtkPartitionedDataSetCollection
        to GGD subset definitions.

        Args:
            grid_ggd: GGD grid structure for the current time step, in which the
                grid subsets will be stored.
            vtk_pdsc: the vtkPartitionedDataSetCollection to load partitions from.
        """
        num_partitions = vtk_pdsc.GetNumberOfPartitionedDataSets()
        grid_ggd.grid_subset.resize(num_partitions)

        used_offsets = dict.fromkeys(
            range(len(grid_ggd.space[0].objects_per_dimension)), 0
        )

        for p in range(num_partitions):
            subset = grid_ggd.grid_subset[p]
            name = vtk_pdsc.GetMetaData(p).Get(vtk_pdsc.NAME())
            pds = vtk_pdsc.GetPartitionedDataSet(p)
            ugrid = pds.GetPartition(0) if pds.GetNumberOfPartitions() > 0 else None

            if ugrid is None:
                continue

            # Get number of elements for this subset
            n_cells = ugrid.GetNumberOfCells()
            if n_cells == 0 and ugrid.GetNumberOfPoints() > 0:
                cell_dim = 0
                n_elem = ugrid.GetNumberOfPoints()
            else:
                cell_dim = ugrid.GetCell(0).GetCellDimension() if n_cells > 0 else 0
                n_elem = n_cells
            subset.dimension = cell_dim + 1

            identifier_name = name.lower().replace("-", "_").replace(" ", "_")
            try:
                subset.identifier = identifiers.ggd_subset_identifier[identifier_name]
            except KeyError:
                subset.identifier.name = identifier_name
                subset.identifier.index = 0
                subset.identifier.description = "Unknown subset identifier"

            subset.element.resize(n_elem)
            for i, element in enumerate(subset.element):
                element.object.resize(1)
                obj = element.object[0]
                obj.space = 1
                obj.dimension = subset.dimension
                if cell_dim == 0:
                    obj.index = i + 1
                else:
                    obj.index = used_offsets[cell_dim] + i + 1

            used_offsets[cell_dim] += n_elem
