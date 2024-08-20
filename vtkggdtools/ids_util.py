from vtkggdtools._ids_util import _get_nodes_from_path, _recursive_ggd_path_search


def get_arrays_from_ids(ids, get_empty_arrays=False):
    """Fetches all GGD scalar and vector arrays that reside in the IDS.

    Args:
        ids: The IDS from which to fetch GGD arrays

    Returns:
        scalar_array_list: The GGD scalar arrays (real & complex)
        vector_array_list: The GGD vector arrays (normal & rphiz)
    """
    # Recursively search the IDS for GGD paths
    scalar_array_paths = []
    vector_array_paths = []
    _recursive_ggd_path_search(
        ids.metadata,
        scalar_array_paths,
        vector_array_paths,
    )

    # Find scalar and vector GGD arrays in the IDS from the paths
    scalar_array_list = []
    vector_array_list = []
    for scalar_path in scalar_array_paths:
        scalar_array_list.extend(_get_nodes_from_path(ids, scalar_path))

    for vector_path in vector_array_paths:
        vector_array_list.extend(_get_nodes_from_path(ids, vector_path))

    # Remove empty arrays
    if not get_empty_arrays:
        scalar_array_list = [array for array in scalar_array_list if len(array) > 0]
        vector_array_list = [array for array in vector_array_list if len(array) > 0]

    return scalar_array_list, vector_array_list
