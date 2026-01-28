import imas
from imas.ids_data_type import IDSDataType
from imas.ids_struct_array import IDSStructArray
from imas.ids_toplevel import IDSToplevel

from imas_paraview._ids_util import _get_nodes_from_path


def get_arrays_from_ids(
    ids,
    ggd_idx=None,
    get_empty_arrays=False,
    scalar_array_paths=None,
    vector_array_paths=None,
    create_empty_structs=False,
):
    """Fetches GGD scalar and vector arrays that reside in the IDS. If lists of
    IDSPaths of GGD arrays are provided through scalar_array_paths or
    vector_array_paths, only these will be loaded. If either of them is not provided,
    instead all GGD arrays in the IDS will be returned.

    Args:
        ids: The IDS from which to fetch GGD arrays
        ggd_idx: The GGD time step to load. Defaults to None, which corresponds with
            loading all timesteps.
        get_empty_arrays (bool): Whether to return empty GGD arrays
        scalar_array_paths: A list of IDSPaths of GGD scalar arrays to search through.
        vector_array_paths: A list of IDSPaths of GGD vector arrays to search through.
        create_empty_structs: If this flag is enabled and an empty structure is
            encountered through which has to be traversed to reach a GGD array, said
            structure is resized to have length 1.
    Returns:
        A tuple with two elements.

        - scalar_array_list: The GGD scalar arrays (real & complex)
        - vector_array_list: The GGD vector arrays (normal & rphiz)
    """
    if scalar_array_paths is None or vector_array_paths is None:
        # Recursively search the IDS for GGD paths
        scalar_array_paths = []
        vector_array_paths = []
        recursive_ggd_path_search(
            ids.metadata,
            scalar_array_paths,
            vector_array_paths,
        )

    # Find scalar and vector GGD arrays in the IDS from the paths
    scalar_array_list = []
    vector_array_list = []
    for scalar_path in scalar_array_paths:
        scalar_array_list.extend(
            _get_nodes_from_path(
                ids, scalar_path, get_empty_arrays, ggd_idx, create_empty_structs
            )
        )

    for vector_path in vector_array_paths:
        vector_array_list.extend(
            _get_nodes_from_path(
                ids, vector_path, get_empty_arrays, ggd_idx, create_empty_structs
            )
        )

    return scalar_array_list, vector_array_list


def recursive_ggd_path_search(
    quantity_metadata, scalar_array_paths, vector_array_paths
):
    """Recursively searches through the metadata of an IDS node for scalar GGD arrays
    (real & complex) and vector GGD arrays (regular and rphiz), and appends the paths of
    these to the scalar_array_paths and vector_array_paths respectively.

    Args:
        quantity_metadata: The metadata of an IDS node
        scalar_array_paths: The IDSPaths of GGD scalar arrays (real & complex)
        vector_array_paths: The IDSPaths of GGD vector arrays (regular and rphiz)
    """
    for subquantity_metadata in quantity_metadata:
        if subquantity_metadata.data_type == IDSDataType.STRUCT_ARRAY:
            # Get scalar and complex scalar array quantities
            if subquantity_metadata.structure_reference in [
                "generic_grid_scalar",
                "generic_grid_scalar_complex",
            ]:
                scalar_array_paths.append(subquantity_metadata.path)

            # Get vector and rzphi-vector array quantities
            # From DDv4 onward `generic_grid_vector_components_rzphi` will be
            # replaced by `generic_grid_vector_components_rphiz`
            elif subquantity_metadata.structure_reference in [
                "generic_grid_vector_components",
                "generic_grid_vector_components_rzphi",
                "generic_grid_vector_components_rphiz",
            ]:
                vector_array_paths.append(subquantity_metadata.path)

        recursive_ggd_path_search(
            subquantity_metadata,
            scalar_array_paths,
            vector_array_paths,
        )


def is_child_of_time_dependent_aos(node):
    """Returns True if the provided IDS node is a child of a time-dependent AoS,
    and False otherwise.

    Examples of children of a time-dependent AoSs:
    - equilibrium/time_slice[0]
    - core_sources/source[0]/global_quantities[1]
    - core_profiles/profiles_1d[2]
    """
    parent = imas.util.get_parent(node)
    return (
        hasattr(node.metadata, "coordinate1")
        and node.metadata.coordinate1.is_time_coordinate
        and hasattr(parent.metadata, "coordinate1")
        and parent.metadata.coordinate1.is_time_coordinate
    )


def create_name_recursive(node):
    """Generates a name for an IDS node. The parents of the node are
    searched recursively until the IDS toplevel is reached. The name of the metadata
    of each parent node is stored as well as the identifier, name or labels of the
    node, which are added in brackets, if applicable.

    Args:
        node: The IDS node

    Returns:
        Name of the IDS node
    """
    name_current_node = node.metadata.name
    name = ""
    parent = imas.util.get_parent(node)
    # skip time slice quantity
    if not is_child_of_time_dependent_aos(node):
        name_appendix = ""

        # Check if node has an identifier.name
        if hasattr(node, "identifier") and hasattr(node.identifier, "name"):
            name_appendix = str(node.identifier.name).strip()

        # Check if node has a name
        elif hasattr(node, "name"):
            name_appendix = str(node.name).strip()

        # Check if node has a label
        elif hasattr(node, "label"):
            name_appendix = str(node.label.value).strip()

        # Add identifier if AoS has no name to ensure uniqueness
        elif isinstance(parent, IDSStructArray):
            # NOTE: This is O(N) in the size of the AoS
            for index, child in enumerate(parent):
                if child is node:
                    name_appendix = f"#{index + 1}"
                    break

        # Add identifier/name/label in between brackets to the full name
        if name_appendix != "":
            name = f"{name_current_node.capitalize()} ({name_appendix.capitalize()})"
        else:
            name = name_current_node.capitalize()

    if parent.metadata is node.metadata:
        parent = imas.util.get_parent(parent)
    if not isinstance(parent, IDSToplevel):
        name = f"{create_name_recursive(parent)} {name}"

    name = name.strip()
    return name
