import logging

import numpy as np
from imas.ids_data_type import IDSDataType
from imas.ids_struct_array import IDSStructArray
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonCore import vtkDoubleArray
from vtkmodules.vtkCommonDataModel import vtkCellData, vtkPointData, vtkUnstructuredGrid

from imas_paraview.ids_util import (
    create_name_recursive,
    get_arrays_from_ids,
    recursive_ggd_path_search,
)

# We'll need these below when we create some units manually:
from imas_paraview.util import format_units

logger = logging.getLogger("imas_paraview")
u_pre = "["
u_post = "]"


class PlasmaStateReader:
    def __init__(self, ids):
        """Initializes plasma state reader and retrieves all filled GGD scalar and
        vector arrays from the IDS.

        Args:
            ids: The IDS to load GGD arrays from
        """

        # _cache stores names for each node to avoid recomputing them. It checks if
        # a node's name is already cached before generating it, speeding up the process
        # and ensuring names are computed only once.
        self._cache = {}
        self._ids = ids
        self.scalar_array_list = []
        self.vector_array_list = []

    def load_paths_from_ids(self, ggd_idx=0, return_empty=False):
        """Retrieves scalar and vector array paths from the IDS metadata by performing a
        recursive search through GGD paths.

        Args:
            ggd_idx: Retrieve filled paths of a specific time step.
            return_empty: Whether to return the paths of empty GGD arrays. Defaults to
                False.

        Returns:
            scalar_array_paths: A list of paths of GGD scalar arrays in the IDS
            vector_array_paths: A list of paths of GGD vector arrays in the IDS
        """

        logger.debug("Retrieving GGD paths from IDS metadata")
        all_scalar_paths = []
        all_vector_paths = []
        filled_scalar_paths = []
        filled_vector_paths = []
        recursive_ggd_path_search(
            self._ids.metadata, all_scalar_paths, all_vector_paths
        )

        if return_empty:
            logger.debug("Retrieving all GGD arrays, including empty ones.")
        else:
            logger.debug("Only retrieving filled GGD arrays.")
            filled_scalar_paths = self._get_filled_arrays(
                all_scalar_paths, ggd_idx=ggd_idx
            )
            filled_vector_paths = self._get_filled_arrays(
                all_vector_paths, ggd_idx=ggd_idx
            )

        return (
            all_scalar_paths,
            all_vector_paths,
            filled_scalar_paths,
            filled_vector_paths,
        )

    def load_arrays_from_path(self, ggd_idx, scalar_array_paths, vector_array_paths):
        """Retrieves scalar and vector arrays from the IDS based on the provided time
        index and array paths, and stores them in the respective lists.

        Args:
            ggd_idx: The time index for which to load GGD arrays
            scalar_array_paths: A list of paths of GGD scalar arrays in the IDS
            vector_array_paths: A list of paths of GGD vector arrays in the IDS
        """
        logger.debug("Retrieving GGD arrays from IDS")
        self.scalar_array_list, self.vector_array_list = get_arrays_from_ids(
            self._ids,
            ggd_idx=ggd_idx,
            scalar_array_paths=scalar_array_paths,
            vector_array_paths=vector_array_paths,
        )
        logger.debug(
            f"Found {len(self.scalar_array_list)} scalar arrays and "
            f"{len(self.vector_array_list)} vector arrays in the IDS."
        )

    def read_plasma_state(self, subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:
        """Reads plasma state data arrays from the ggd node in the IDS. These arrays are
        added as point data or cell data to the unstructured grid.

        Args:
            subset_idx: an index into grid_ggd/grid_subset AoS
            ugrid: the unstructured grid instance
        """
        # TODO: GGD-fast

        # Read scalar arrays
        logger.debug("Converting scalar GGD arrays to VTK field data")
        for scalar_array in self.scalar_array_list:
            name = self._create_name_with_units(scalar_array)
            self._add_aos_scalar_array_to_vtk_field_data(
                scalar_array, subset_idx, name, ugrid
            )

        # Read vector arrays
        logger.debug("Converting Vector GGD arrays to VTK field data")
        for vector_array in self.vector_array_list:
            name = self._create_name_with_units(vector_array)
            self._add_aos_vector_array_to_vtk_field_data(
                vector_array, subset_idx, name, ugrid
            )

    def _get_filled_arrays(self, path_list, ggd_idx):
        """Look through a list of IDSPaths containing GGD arrays and return a list
        with only the IDSPaths of GGD arrays which are filled.

        Args:
            path_list: A list of IDSPaths containing GGD arrays.
            ggd_idx: Index for selecting specific GGD array.

        Returns:
            A list containing only the paths of filled GGD arrays.
        """
        filled_arrays = []
        for path in path_list:
            if self._is_path_filled(path, ggd_idx):
                filled_arrays.append(path)
        return filled_arrays

    def _is_path_filled(self, path, ggd_idx):
        """Checks if the provided GGD path for a given time index is filled.

        Args:
            path: IDSPath containing a GGD array.
            ggd_idx: Index for selecting specific GGD array.
        """

        def _check_subnode(node, parts, idx):
            """
            Recursive helper function
            """

            if idx == len(parts):
                return node.size > 0
            part = parts[idx]
            if isinstance(node, IDSStructArray):
                name = getattr(node.metadata, "name", None)
                if name == "ggd":
                    if len(node) <= ggd_idx:
                        return False
                    node = node[ggd_idx]
                    return _check_subnode(node, parts, idx)

                elif len(node) >= 1:
                    for subnode in node:
                        if _check_subnode(subnode, parts, idx):
                            return True
            else:
                node = node[part]
                return _check_subnode(node, parts, idx + 1)

        node = self._ids[path.parts[0]]
        return _check_subnode(node, path.parts, 1)

    def _create_name_with_units(self, array):
        """Creates a name for the GGD array based on its path and units.

        Args:
            array: The ggd scalar or vector array to create the name for

        Returns:
            name_with_units: The name and units of the provided GGD array
        """

        # Format the path in a neat format
        name = self._create_ggd_name(array)

        # Get units for this quantity
        units = format_units(array)

        # Combine name and units
        name_with_units = f"{name} {units}"
        return name_with_units

    def _create_ggd_name(self, array):
        """Generates a name for the GGD array. The name of the metadata
        of each parent node is stored as well as the identifier, name or labels of the
        node, which are added in brackets, if applicable.

        Args:
            array: The GGD array

        Returns:
            Name of the GGD scalar or vector
        """
        node_id = id(array)
        if node_id in self._cache:
            return self._cache[node_id]
        name = create_name_recursive(array)

        self._cache[node_id] = name
        return name

    def _add_scalar_array_to_vtk_field_data(
        self, array: np.ndarray, name: str, ugrid: vtkUnstructuredGrid
    ) -> None:
        """Add a named array as scalars to a vtkUnstructuredGrid instance

        Args:
            array: the numpy array.
            name: the name string.
            ugrid: an instance of vtkUnstructuredGrid
        """

        logger.debug(f"           {name}...")
        point_data: vtkPointData = ugrid.GetPointData()
        num_points = ugrid.GetNumberOfPoints()
        cell_data: vtkCellData = ugrid.GetCellData()
        num_cells = ugrid.GetNumberOfCells()

        # Split complex arrays into real and imaginary parts
        if array.metadata.data_type is IDSDataType.CPX:
            array_real = np.real(array)
            array_imag = np.imag(array)
            vtk_arr_real = dsa.numpyTovtkDataArray(array_real, name=f"{name}_real")
            vtk_arr_imag = dsa.numpyTovtkDataArray(array_imag, name=f"{name}_imag")
        else:
            vtk_arr = dsa.numpyTovtkDataArray(array, name)

        # To see interpolation of point data on cells, just point data is necessary.
        if len(array) == num_points:
            if array.metadata.data_type is IDSDataType.CPX:
                point_data.AddArray(vtk_arr_real)
                point_data.AddArray(vtk_arr_imag)
            else:
                point_data.AddArray(vtk_arr)
        if len(array) == num_cells:
            if array.metadata.data_type is IDSDataType.CPX:
                cell_data.AddArray(vtk_arr_real)
                cell_data.AddArray(vtk_arr_imag)
            else:
                cell_data.AddArray(vtk_arr)

    def _add_aos_scalar_array_to_vtk_field_data(
        self, aos_scalar_node, subset_idx: int, name: str, ugrid: vtkUnstructuredGrid
    ) -> None:
        """Add the array under the aos_scalar_node to the unstructured grid.

        Args:
            aos_scalar_node: A node with scalar values for each grid subset.
            subset_idx: an index into aos_scalar_node
            name: this becomes the array name in VTK
            ugrid: an unstructured grid instance
        """
        logger.debug(f"           {name}...")
        if subset_idx >= len(aos_scalar_node):
            return

        # For wall IDS nodes, edges, cells, volumes in one partition.
        if subset_idx == -1:
            for i in range(4):
                try:
                    if hasattr(aos_scalar_node[i], "values") and len(
                        aos_scalar_node[i].values
                    ):
                        self._add_scalar_array_to_vtk_field_data(
                            aos_scalar_node[i].values, name, ugrid
                        )
                except IndexError:
                    logger.warning(
                        f"           no index {i} for subset {subset_idx}..."
                    )
                except AttributeError:
                    logger.warning(
                        f"           no index {i} for subset {subset_idx}..."
                    )
        else:
            if hasattr(aos_scalar_node[subset_idx], "values") and len(
                aos_scalar_node[subset_idx].values
            ):
                self._add_scalar_array_to_vtk_field_data(
                    aos_scalar_node[subset_idx].values, name, ugrid
                )

    def _add_aos_vector_array_to_vtk_field_data(
        self, aos_vector_node, subset_idx: int, name: str, ugrid: vtkUnstructuredGrid
    ) -> None:
        """Add the array under the aos_vector_node to the unstructured grid.

        Args:
            aos_vector_node: A node with component vectors for each grid subset.
            subset_idx: an index into aos_vector_node
            name: this becomes the array name in VTK
            ugrid: an unstructured grid instance
        """
        logger.debug(f"           {name}...")
        if subset_idx >= len(aos_vector_node):
            return

        point_data: vtkPointData = ugrid.GetPointData()
        num_points = ugrid.GetNumberOfPoints()
        cell_data: vtkCellData = ugrid.GetCellData()
        num_cells = ugrid.GetNumberOfCells()

        # Only add the components that have data:
        components = dict()  # name and values
        for component_name in [
            "radial",
            "diamagnetic",
            "parallel",
            "poloidal",
            "toroidal",
            "r",
            "z",
        ]:
            try:
                values = getattr(aos_vector_node[subset_idx], component_name)
            except (IndexError, AttributeError):
                continue
            if len(values):
                components[component_name] = values

        vtk_arr = vtkDoubleArray()
        vtk_arr.SetName(name)
        vtk_arr.SetNumberOfComponents(len(components))
        num_tuples = 0

        for i, component_name in enumerate(components):
            vtk_arr.SetComponentName(i, component_name.capitalize())
            scalar_arr = dsa.numpyTovtkDataArray(
                components[component_name], name + "-" + component_name.capitalize()
            )
            if num_tuples == 0:
                num_tuples = scalar_arr.GetNumberOfTuples()
                vtk_arr.SetNumberOfTuples(num_tuples)

            vtk_arr.CopyComponent(i, scalar_arr, 0)

        if num_tuples == num_points:
            point_data.AddArray(vtk_arr)
        if num_tuples == num_cells:
            cell_data.AddArray(vtk_arr)
