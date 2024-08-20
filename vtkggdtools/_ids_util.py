from imaspy.ids_data_type import IDSDataType
from imaspy.ids_struct_array import IDSStructArray


def _recursive_ggd_path_search(
    quantity_metadata, scalar_array_paths, vector_array_paths
):
    """Recursively searches through the metadata of an IDS  node for scalar GGD arrays
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
            # TODO: From DDv4 onward `generic_grid_vector_components_rzphi` will be
            # replaced by `generic_grid_vector_components_rphiz`
            # Get vector and rzphi-vector array quantities
            elif subquantity_metadata.structure_reference in [
                "generic_grid_vector_components",
                "generic_grid_vector_components_rzphi",
            ]:
                vector_array_paths.append(subquantity_metadata.path)
            else:
                _recursive_ggd_path_search(
                    subquantity_metadata,
                    scalar_array_paths,
                    vector_array_paths,
                )

        elif subquantity_metadata.data_type == IDSDataType.STRUCTURE:
            # Skip "grid" quantity, this can occur if the grid is stored within the
            # GGD. e.g. in distribution_sources distributions IDSs
            if subquantity_metadata.name != "grid":
                _recursive_ggd_path_search(
                    subquantity_metadata,
                    scalar_array_paths,
                    vector_array_paths,
                )


def _get_nodes_from_path(node, path):
    """Retrieve a list of nodes from a given IDSPath.

    Args:
        node: The starting node to navigate.
        path: An IDSPath to traverse.

    Returns:
        A list of nodes obtained from the specified path.
    """
    return list(_iter_nodes_from_path(node, path.parts))


def _iter_nodes_from_path(node, path_parts):
    """Recursively iterate through nodes in an IDS mpde based on path parts.

    Args:
        node: The current node being traversed.
        path_parts: A list of IDSPath segments (parts).

    Yields:
        The next node in the structure corresponding to the current path part.
    """
    child_node = node[path_parts[0]]
    if len(path_parts) == 1:
        yield child_node
    elif isinstance(child_node, IDSStructArray):
        for structure in child_node:
            yield from _iter_nodes_from_path(structure, path_parts[1:])
    else:
        yield from _iter_nodes_from_path(child_node, path_parts[1:])
