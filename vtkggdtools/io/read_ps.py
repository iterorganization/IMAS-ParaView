import logging

import numpy as np
from imaspy.ids_data_type import IDSDataType
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonCore import vtkDoubleArray
from vtkmodules.vtkCommonDataModel import vtkCellData, vtkPointData, vtkUnstructuredGrid

from vtkggdtools.ids_util import get_arrays_from_ids

# We'll need these below when we create some units manually:
from vtkggdtools.util import format_units

logger = logging.getLogger(__name__)
u_pre = "["
u_post = "]"


class PlasmaStateReader:
    def __init__(self, ids):
        """Initializes plasma state reader and retrieves all filled GGD scalar and
        vector arrays from the IDS.

        Args:
            ids: The IDS to load GGD arrays from
        """
        self.ids = ids

        # Retrieve all GGD scalar and vector arrays from IDS
        logger.debug("Retrieving GGD arrays from IDS")
        self.scalar_array_list, self.vector_array_list = get_arrays_from_ids(ids)

    def read_plasma_state(self, subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:
        """Reads plasma state data arrays from the ggd node in the IDS. These arrays are
        added as point data or cell data to the unstructured grid.

        Args:
            subset_idx: an index into grid_ggd/grid_subset AoS
            ugrid: the unstructured grid instance
        """
        # TODO: properly handle time indexing
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

    def _create_name_with_units(self, array):
        """Creates a name for the GGD array based on its path and units.

        Args:
            ids: The IDS that this array belongs to
            array: The ggd scalar or vector array to create the name for

        Returns:
            name_with_units: The name and units of the provided GGD array
        """

        # Format the path in a neat format
        name = self._create_name(array)

        # Get units for this quantity
        units = format_units(array)

        # Combine name and units
        name_with_units = f"{name} {units}"
        return name_with_units

    def _create_name(self, array):
        """Generates a name for the GGD array. The parent nodes of the array are
        traversed until the IDS toplevel is reached. The name of the metadata of each
        parent node is stored as well as the identifier, name or labels of the node,
        which are added in brackets, if applicable.

        Args:
            array: The ggd scalar or vector array to create the name for

        Returns:
            Name of the ggd scalar or vector
        """
        current_node = array
        name_current_node = ""
        name = []
        # Traverse through the parent nodes until the IDS top level is reached
        while not hasattr(current_node, "ids_properties"):
            previous_name = name_current_node
            name_current_node = current_node.metadata.name

            # Do not nodes contain 'GGD' or 'time_slice' to the final name.
            # Also, avoid adding a name if it's the same as the previous name.
            # This occurs for example when searching through edge_sources. First it
            # The loop may first encounter 'source[0]' and then 'source', both of which
            # have the same metadata name. We only include 'source[0]' because this
            # entry contains the identifier/label name.
            if (
                "ggd" != name_current_node
                and "time_slice" not in name_current_node
                and previous_name != name_current_node
            ):
                name_appendix = ""

                # Check if node has an identifier.name
                if hasattr(current_node, "identifier") and hasattr(
                    current_node.identifier, "name"
                ):
                    name_appendix = str(current_node.identifier.name).strip()

                # Check if node has a name
                elif hasattr(current_node, "name"):
                    name_appendix = str(current_node.name).strip()

                # Check if node has a label
                elif hasattr(current_node, "label"):
                    name_appendix = str(current_node.label.value).strip()

                # Add identifier/name/label in between brackets to the full name
                if name_appendix:
                    full_name = f"{name_current_node} ({name_appendix})"
                else:
                    full_name = name_current_node

                name.append(full_name)

            # Set current node to the parent node
            current_node = current_node._parent
        return " ".join(reversed(name))

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
                    logger.warn(f"           no index {i} for subset {subset_idx}...")
                except AttributeError:
                    logger.warn(f"           no index {i} for subset {subset_idx}...")
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
