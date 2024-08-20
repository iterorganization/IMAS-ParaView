from imaspy.ids_data_type import IDSDataType
from imaspy.util import is_lazy_loaded


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
    ids_is_lazy_loaded = is_lazy_loaded(ids)
    _recursive_array_search(
        ids, scalar_array_list, vector_array_list, get_empty_arrays, ids_is_lazy_loaded
    )
    return scalar_array_list, vector_array_list


def _recursive_array_search(
    quantity, scalar_array_list, vector_array_list, get_empty_arrays, ids_is_lazy_loaded
):
    """Recursively searches through the IDS node for scalar (real & complex) and
    vector arrays, and appends these to the scalar_array_list and vector_array_list
    respectively.

    Args:
        quantity: The IDS node to search from
        scalar_array_list: The GGD scalar arrays (real & complex)
        vector_array_list: The GGD vector arrays
        get_empty_arrays (bool): Whether to return empty arrays
        ids_is_lazy_loaded (bool): Whether the IDS was lazy loaded
    """
    for subquantity in quantity:
        if not ids_is_lazy_loaded:
            # Only checkout subquantity if it is non-empty
            if not subquantity.has_value and not get_empty_arrays:
                continue

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
                    ids_is_lazy_loaded,
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
                    ids_is_lazy_loaded,
                )
