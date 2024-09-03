from imaspy.ids_data_type import IDSDataType

from vtkggdtools._ids_util import _get_nodes_from_path


def get_arrays_from_ids(
    ids,
    ggd_idx=None,
    get_empty_arrays=False,
    scalar_array_paths=None,
    vector_array_paths=None,
):
    """Fetches all GGD scalar and vector arrays that reside in the IDS.

    Args:
        ids: The IDS from which to fetch GGD arrays
        get_arrays_from_ids (bool): Whether to return empty GGD arrays
        ggd_idx: The GGD time step to load. Defaults to None, which corresponds with
        loading all timesteps.
    Returns:
        scalar_array_list: The GGD scalar arrays (real & complex)
        vector_array_list: The GGD vector arrays (normal & rphiz)
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
            _get_nodes_from_path(ids, scalar_path, get_empty_arrays, ggd_idx)
        )

    for vector_path in vector_array_paths:
        vector_array_list.extend(
            _get_nodes_from_path(ids, vector_path, get_empty_arrays, ggd_idx)
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
