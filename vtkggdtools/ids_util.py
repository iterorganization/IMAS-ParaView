from imaspy.ids_data_type import IDSDataType


def get_arrays_from_ids(ids, get_empty_arrays=False):
    """Fetches all GGD scalar and vector arrays that reside in the IDS.

    Args:
        get_arrays_from_ids (bool): Whether to return empty arrays

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
        # Only checkout subquantity if it is non-empty
        if not subquantity.has_value and not get_empty_arrays:
            continue

        metadata = subquantity.metadata
        # If subquantity is a struct array
        if metadata.data_type == IDSDataType.STRUCT_ARRAY:
            # Get scalar array quantities
            if metadata.structure_reference == "generic_grid_scalar":
                scalar_array_list.append(subquantity)
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
