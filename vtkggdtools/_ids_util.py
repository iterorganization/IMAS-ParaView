from imaspy.ids_data_type import IDSDataType
from imaspy.ids_struct_array import IDSStructArray


def _recursive_ggd_path_search(
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

        _recursive_ggd_path_search(
            subquantity_metadata,
            scalar_array_paths,
            vector_array_paths,
        )


def _get_nodes_from_path(node, path, get_empty_arrays, ggd_idx=None):
    """Retrieve a list of nodes from a given IDSPath.

    Args:
        node: The starting node to navigate.
        path: An IDSPath to traverse.
        get_empty_arrays (bool): Whether to return empty GGD arrays
        ggd_idx: The GGD time step to load. Defaults to None, which corresponds with
        loading all timesteps.

    Returns:
        A list of nodes obtained from the specified path.
    """
    return list(_iter_nodes_from_path(node, path.parts, get_empty_arrays, ggd_idx))


def _iter_nodes_from_path(node, path_parts, get_empty_arrays, ggd_idx):
    """Recursively iterate through nodes of an IDS node based on path parts.

    Args:
        node: The current node being traversed.
        path_parts: A list of IDSPath segments (parts).
        get_empty_arrays (bool): Whether to return empty GGD arrays
        ggd_idx: The GGD time step to load.

    Yields:
        The next node in the structure corresponding to the current path part.
    """
    child_node = node[path_parts[0]]
    if len(path_parts) == 1:
        # The path_parts refer to nodes that have a defined length, such as struct
        # arrays
        if len(child_node) > 1 or get_empty_arrays:
            yield child_node
    elif isinstance(child_node, IDSStructArray):
        # Only load specific timeidx from ggd node
        if ggd_idx is not None and path_parts[0] == "ggd":
            if len(child_node) > ggd_idx:
                structure = child_node[ggd_idx]
                yield from _iter_nodes_from_path(
                    structure, path_parts[1:], get_empty_arrays, ggd_idx
                )
        else:
            for structure in child_node:
                yield from _iter_nodes_from_path(
                    structure, path_parts[1:], get_empty_arrays, ggd_idx
                )
    else:
        yield from _iter_nodes_from_path(
            child_node, path_parts[1:], get_empty_arrays, ggd_idx
        )
