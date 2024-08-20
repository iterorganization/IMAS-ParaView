from imaspy.ids_data_type import IDSDataType
from imaspy.ids_struct_array import IDSStructArray


def _recursive_array_search(
    quantity, scalar_array_list, vector_array_list, get_empty_arrays
):
    """Recursively searches through the IDS node for scalar (real & complex) and
    vector arrays, and appends these to the scalar_array_list and vector_array_list
    respectively.

    Args:
        quantity: The IDS node to search from
        scalar_array_list: The GGD scalar arrays (real & complex)
        vector_array_list: The GGD vector arrays
        get_empty_arrays (bool): Whether to return empty arrays
    """
    for subquantity in quantity:
        # # Only checkout subquantity if it is non-empty
        # if not subquantity.has_value and not get_empty_arrays:
        #     continue

        metadata = subquantity.metadata
        # If subquantity is a struct array
        if metadata.data_type == IDSDataType.STRUCT_ARRAY:
            # Get scalar and complex scalar array quantities
            if metadata.structure_reference in [
                "generic_grid_scalar",
                "generic_grid_scalar_complex",
            ]:
                scalar_array_list.append(subquantity)
            # TODO: From DDv4 onward `generic_grid_vector_components_rzphi` will be
            # replaced by `generic_grid_vector_components_rphiz`
            # Get vector and rzphi-vector array quantities
            elif metadata.structure_reference in [
                "generic_grid_vector_components",
                "generic_grid_vector_components_rzphi",
            ]:
                vector_array_list.append(subquantity)
            # Recursively search
            else:
                _recursive_array_search(
                    subquantity,
                    scalar_array_list,
                    vector_array_list,
                    get_empty_arrays,
                )

        # If subquantity is a structure
        elif metadata.data_type == IDSDataType.STRUCTURE:
            # Skip "grid" quantity, this can occur if the grid is stored within the
            # GGD. e.g. in distribution_sources distributions IDSs
            if metadata.name != "grid":
                _recursive_array_search(
                    subquantity,
                    scalar_array_list,
                    vector_array_list,
                    get_empty_arrays,
                )


def _recursive_path_search(quantity, scalar_array_paths, vector_array_paths):
    """Recursively searches through the IDS node for scalar (real & complex) and
    vector arrays, and appends these to the scalar_array_list and vector_array_list
    respectively.

    Args:
        quantity: The IDS node to search from
        scalar_array_list: The GGD scalar arrays (real & complex)
        vector_array_list: The GGD vector arrays
        get_empty_arrays (bool): Whether to return empty arrays
    """
    for _, subquantity in quantity._children.items():
        # If subquantity is a struct array
        if subquantity.data_type == IDSDataType.STRUCT_ARRAY:
            # Get scalar and complex scalar array quantities
            if subquantity.structure_reference in [
                "generic_grid_scalar",
                "generic_grid_scalar_complex",
            ]:
                scalar_array_paths.append(subquantity.path)
            # TODO: From DDv4 onward `generic_grid_vector_components_rzphi` will be
            # replaced by `generic_grid_vector_components_rphiz`
            # Get vector and rzphi-vector array quantities
            elif subquantity.structure_reference in [
                "generic_grid_vector_components",
                "generic_grid_vector_components_rzphi",
            ]:
                vector_array_paths.append(subquantity.path)
            # Recursively search
            else:
                _recursive_path_search(
                    subquantity,
                    scalar_array_paths,
                    vector_array_paths,
                )

        # If subquantity is a structure
        elif subquantity.data_type == IDSDataType.STRUCTURE:
            # Skip "grid" quantity, this can occur if the grid is stored within the
            # GGD. e.g. in distribution_sources distributions IDSs
            if subquantity.name != "grid":
                _recursive_path_search(
                    subquantity,
                    scalar_array_paths,
                    vector_array_paths,
                )


def _get_nodes_from_path(node, path):
    return list(_iter_nodes_from_path(node, path.parts))


def _iter_nodes_from_path(node, path_parts):
    child_node = node[path_parts[0]]
    if len(path_parts) == 1:
        yield child_node
    elif isinstance(child_node, IDSStructArray):
        for structure in child_node:
            yield from _iter_nodes_from_path(structure, path_parts[1:])
    else:
        yield from _iter_nodes_from_path(child_node, path_parts[1:])
