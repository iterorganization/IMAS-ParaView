import logging
import re

import numpy as np
from imaspy.ids_data_type import IDSDataType
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonCore import vtkDoubleArray
from vtkmodules.vtkCommonDataModel import vtkCellData, vtkPointData, vtkUnstructuredGrid

# We'll need these below when we create some units manually:
from vtkggdtools.util import format_units

logger = logging.getLogger(__name__)
u_pre = "["
u_post = "]"


def read_plasma_state(
    ids,
    subset_idx: int,
    ugrid: vtkUnstructuredGrid,
) -> None:
    """Reads plasma state data arrays from the ggd node. These arrays are added as point
    data or cell data to the unstructured grid.
    This is just a dispatching routine, the real work is done in the called
    funtions read_<IDS_name>(ids_obj, aos_index_values, subset_idx, ugrid)

    Args:
        ids: an ids object
        subset_idx: an index into grid_ggd/grid_subset AoS
        ugrid: the unstructured grid instance

    Returns:
        None
    """
    # TODO: read_plasma_state is getting called many times for every subset
    # in vtkggdtools.py, which means some functionality (such as getting the scalar and
    # vector arrays from the IDS) is done multiple times unnecessarily. In the future,
    # change this to an OOP approach so this is only done once to improve loading times.
    read_ids(ids, subset_idx, ugrid)


def read_ids(ids, subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:
    """Reads IDS and converts GGD to VTK data

    Args:
        ids: The IDS to be converted
        subset_idx: an index into grid_ggd/grid_subset AoS
        ugrid: the unstructured grid instance
    """

    # TODO: properly handle time indexing
    # TODO: GGD-fast

    # Retrieve all GGD scalar and vector arrays from IDS
    logger.debug("Retrieving GGD arrays from IDS")
    scalar_array_list, vector_array_list = get_arrays_from_ids(
        ids, get_empty_arrays=False
    )

    # Read scalar arrays
    logger.debug("Converting scalar GGD arrays to VTK field data")
    for scalar_array in scalar_array_list:
        name = _get_name(ids, scalar_array)
        _add_aos_scalar_array_to_vtk_field_data(scalar_array, subset_idx, name, ugrid)

    # Read vector arrays
    logger.debug("Converting Vector GGD arrays to VTK field data")
    for vector_array in vector_array_list:
        name = _get_name(ids, vector_array)
        _add_aos_vector_array_to_vtk_field_data(vector_array, subset_idx, name, ugrid)


def get_arrays_from_ids(ids, get_empty_arrays):
    """Fetches all GGD scalar and vector arrays that reside in the IDS.

    Args:
        ids: The IDS to retrieve the arrays from

    Returns:
        scalar_array_list: The GGD scalar arrays (real & complex)
        vector_array_list: The GGD vector arrays
    """
    scalar_array_list = []
    vector_array_list = []
    _recursive_array_search(ids, scalar_array_list, vector_array_list, get_empty_arrays)
    return scalar_array_list, vector_array_list


def _recursive_array_search(
    quantity, scalar_array_list, vector_array_list, get_empty_arrays
):
    """Recursively searches through the IDS node for scalar (real & complex) and vector
    arrays, and appends these to the scalar_array_list and vector_array_list
    respectively.

    Args:
        quantity: The IDS node to search from
        scalar_array_list: The GGD scalar arrays (real & complex)
        vector_array_list: The GGD vector arrays
        get_empty_arrays (bool): Whether to return empty arrays
    """
    for subquantity in quantity:
        # Only checkout subquantity if it is non-empty
        if not subquantity.has_value and not get_empty_arrays:
            continue

        metadata = subquantity.metadata
        # If subquantity is a struct array
        if metadata.data_type == IDSDataType.STRUCT_ARRAY:
            # Get scalar array quantities
            if metadata.structure_reference == "generic_grid_scalar":
                scalar_array_list.append(subquantity)
            # TODO: waves IDS can have complex arrays that are being read in by
            # read_waves. This gives a warning for the pytests, due to casting complex
            # to real. Need to investigate this further
            # Get complex scalar array quantity
            elif metadata.structure_reference == "generic_grid_scalar_complex":
                scalar_array_list.append(subquantity)
            # Get vector array quantities
            elif metadata.structure_reference == "generic_grid_vector_components":
                vector_array_list.append(subquantity)
            # Get rzphi-vector array quantities
            elif metadata.structure_reference == "generic_grid_vector_components_rzphi":
                vector_array_list.append(subquantity)
            # Recursively search
            else:
                _recursive_array_search(
                    subquantity, scalar_array_list, vector_array_list, get_empty_arrays
                )

        # If subquantity is a structure
        elif metadata.data_type == IDSDataType.STRUCTURE:
            # Skip "grid" quantity, this can occur if the grid is stored within the GGD
            # e.g. in distribution_sources distributions IDSs
            if metadata.name != "grid":
                _recursive_array_search(
                    subquantity, scalar_array_list, vector_array_list, get_empty_arrays
                )


def _get_name(ids, array):
    """Creates a name for the GGD array based on its path and units.

    Args:
        ids: The IDS that this array belongs to
        array: The specific scalar or vector array to create the name for

    Returns:
        name_with_units: The name and units of the provided GGD array
    """
    # TODO: properly format the names
    # Get full path of array vector
    path = array._path

    # Format the path in a neat format
    name = _create_name_from_path(ids, path)

    # Get units for this quantity
    units = format_units(array)

    # Combine name and units
    name_with_units = f"{name} {units}"

    return name_with_units


def _create_name_from_path(ids, array_path):
    """Rewrites the path of the array to a neat format.

    Example:
        Given path="source[4]/ggd[0]/ion[2]/state[13]/energy"

        This gets split into its constituents:
            ["source[4]", "ggd[0]", "ion[2]", "state[13]", "energy"]

        The part containing ggd gets removed:
            ["source[4]", "ion[2]", "state[13]", "energy"]

        For each part that contains brackets, the identifier name or label is
        substituted:
            ["Atomic ionization", "Ar", "Ar+14", "energy"]

        These get appended to result in the final returned name:
            name = "Atomic ionization Ar Ar+14 energy"

    Args:
        ids: The IDS that this array belongs to
        array_path: The full path of the GGD array

    Returns:
        name: The name of the GGD array
    """

    # Split path into parts
    split_path, accum_split_name = _split_and_accumulate_path(array_path)

    # Remove the path part containing ggd, as this is not needed for the name
    split_path, accum_split_name = _remove_ggd_from_split_path(
        split_path, accum_split_name
    )

    name_segments = []
    for i in range(len(split_path)):
        split_part = split_path[i]
        split_accum_part = accum_split_name[i]

        # Check if name can be rewritten using identifier.name or label values
        if "[" in split_part and "]" in split_part:

            # Check if node has an identifier.name
            if hasattr(ids[split_accum_part], "identifier") and hasattr(
                ids[split_accum_part].identifier, "name"
            ):
                name_segments.append(
                    f"{str(ids[split_accum_part].metadata.name).strip()} "
                    f"({str(ids[split_accum_part].identifier.name).strip()})"
                )

            # Check if node has a name
            elif hasattr(ids[split_accum_part], "name"):
                name_segments.append(
                    f"{str(ids[split_accum_part].metadata.name).strip()} "
                    f"({str(ids[split_accum_part].name).strip()})"
                )

            # Check if node has a label
            elif hasattr(ids[split_accum_part], "label"):
                name_segments.append(
                    f"{str(ids[split_accum_part].metadata.name).strip()} "
                    f"({str(ids[split_accum_part].label.value).strip()})"
                )

            # Otherwise just use the split part as is
            else:
                name_segments.append(split_part)
        else:
            name_segments.append(split_part)

    name = " ".join(name_segments)
    return name


def _remove_ggd_from_split_path(split_path, accum_split_path):
    """Removes the element from split_path that matches "ggd[idx]" and also removes
    the element at the same index from accum_split_path.

    Example:
        Given the following split_path and accum_split_path:

        split_path=["source[11]", "ggd[0]", "neutral[3]", "momentum"]

        accum_split_path=[
            "source[11]",
            "source[11]/ggd[0]",
            "source[11]/ggd[0]/neutral[3]",
            "source[11]/ggd[0]/neutral[3]/momentum"]

        We obtain:

        split_path_no_ggd = ["source[11]", "neutral[3]", "momentum"]

        accum_split_path_no_ggd = [
            "source[11]",
            "source[11]/ggd[0]/neutral[3]",
            "source[11]/ggd[0]/neutral[3]/momentum"]

    Args:
        split_path: The split path of the GGD array
        accum_split_path: The splitted and accumulated path of the GGD array

    Returns:
        split_path_no_ggd: The split path of the GGD array without the GGD part
        accum_split_path_no_ggd: The splitted and accumulated path of the GGD array,
        without GGD part
    """

    # Define the pattern to match "ggd[idx]"
    ggd_pattern = re.compile(r"ggd\[\d+\]")

    # Required for wall IDS
    description_ggd_pattern = re.compile(r"description_ggd\[\d+\]")

    # Required for equilbrium IDS
    time_slice_pattern = re.compile(r"time_slice\[\d+\]")

    # Create new lists for the output
    split_path_no_ggd = []
    accum_split_path_no_ggd = []

    # Iterate through both lists and only add elements that do not match the pattern
    for i in range(len(split_path)):
        if (
            not ggd_pattern.match(split_path[i])
            and not description_ggd_pattern.match(split_path[i])
            and not time_slice_pattern.match(split_path[i])
        ):
            split_path_no_ggd.append(split_path[i])
            accum_split_path_no_ggd.append(accum_split_path[i])

    return split_path_no_ggd, accum_split_path_no_ggd


def _split_and_accumulate_path(path):
    """Splits the path and returns each path segment separately and accumulated.

    Example:
        Given input path:

        "ggd[0]/neutral[0]/state[0]/density"

        We get the following output:

        split_path = ["ggd[0]", "neutral[0]", "state[0]", "density_fast"]

        accumulated_split_path = [
        "ggd[0]",
        "ggd[0]/neutral[0]",
        "ggd[0]/neutral[0]/state[0]",
        "ggd[0]/neutral[0]/state[0]/density_fast"]

    Args:
        path: The path of a GGD array

    Returns:
        split_path: A list of strings containing the path segments
        accumulated_split_path: A list of strings containing the accumulated path
        segments
    """
    split_path = path.split("/")
    accumulated_split_path = []
    for i in range(1, len(split_path) + 1):
        accumulated_split_path.append("/".join(split_path[:i]))
    return split_path, accumulated_split_path


def _add_scalar_array_to_vtk_field_data(
    array: np.ndarray, name: str, ugrid: vtkUnstructuredGrid
) -> None:
    """
    Add a named array as scalars to a vtkUnstructuredGrid instance
    :param array: the numpy array.
    :param name: the name string.
    :param ugrid: an instance of vtkUnstructuredGrid
    :return: None
    """
    logger.debug(f"           {name}...")
    point_data: vtkPointData = ugrid.GetPointData()
    num_points = ugrid.GetNumberOfPoints()
    cell_data: vtkCellData = ugrid.GetCellData()
    num_cells = ugrid.GetNumberOfCells()

    # To see interpolation of point data on cells, just point data is necessary.
    vtk_arr = dsa.numpyTovtkDataArray(array, name)
    if len(array) == num_points:
        point_data.AddArray(vtk_arr)
    if len(array) == num_cells:
        cell_data.AddArray(vtk_arr)


def _add_aos_scalar_array_to_vtk_field_data(
    aos_scalar_node, subset_idx: int, name: str, ugrid: vtkUnstructuredGrid
):
    """
    Add the array under the aos_scalar_node to the unstructured grid.
    :param aos_scalar_node: A node with scalar values for each grid subset.
    :param subset_idx: an index into aos_scalar_node
    :param name: this becomes the array name in VTK
    :param ugrid: an unstructured grid instance
    :return: None
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
                    _add_scalar_array_to_vtk_field_data(
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
            _add_scalar_array_to_vtk_field_data(
                aos_scalar_node[subset_idx].values, name, ugrid
            )


def _add_aos_vector_array_to_vtk_field_data(
    aos_vector_node, subset_idx: int, name: str, ugrid: vtkUnstructuredGrid
):
    """
    Add the array under the aos_vector_node to the unstructured grid.
    :param aos_vector_node: A node with component vectors for each grid subset.
    :param subset_idx: an index into aos_vector_node
    :param name: this becomes the array name in VTK
    :param ugrid: an unstructured grid instance
    :return: None
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
