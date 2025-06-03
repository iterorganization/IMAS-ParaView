from imas.ids_struct_array import IDSStructArray


def _get_nodes_from_path(
    node, path, get_empty_arrays, ggd_idx=None, create_empty_structs=False
):
    """Retrieve a list of nodes from a given IDSPath.

    Args:
        node: The starting node to navigate.
        path: An IDSPath to traverse.
        get_empty_arrays (bool): Whether to return empty GGD arrays
        ggd_idx: The GGD time step to load. Defaults to None, which corresponds with
        loading all timesteps.
        create_empty_structs: If required to traverse through an empty structure to
            reach a GGD array, resize this structure to have length 1.

    Returns:
        A list of nodes obtained from the specified path.
    """
    return list(
        _iter_nodes_from_path(
            node, path.parts, get_empty_arrays, ggd_idx, create_empty_structs
        )
    )


def _iter_nodes_from_path(
    node, path_parts, get_empty_arrays, ggd_idx, create_empty_structs
):
    """Recursively iterate through nodes of an IDS node based on path parts.

    Args:
        node: The current node being traversed.
        path_parts: A list of IDSPath segments (parts).
        get_empty_arrays (bool): Whether to return empty GGD arrays
        ggd_idx: The GGD time step to load.
        create_empty_structs: If required to traverse through an empty structure to
            reach a GGD array, resize this structure to have length 1.

    Yields:
        The next node in the structure corresponding to the current path part.
    """
    child_node = node[path_parts[0]]
    if len(path_parts) == 1:
        # The path_parts refer to nodes that have a defined length, such as struct
        # arrays
        if len(child_node) >= 1 or get_empty_arrays:
            yield child_node
    elif isinstance(child_node, IDSStructArray):
        # Only load specific timeidx from ggd node
        if ggd_idx is not None and path_parts[0] == "ggd":
            if len(child_node) > ggd_idx:
                structure = child_node[ggd_idx]
                yield from _iter_nodes_from_path(
                    structure,
                    path_parts[1:],
                    get_empty_arrays,
                    ggd_idx,
                    create_empty_structs,
                )
        else:
            # If the structure is empty, we cannot traverse to the GGD arrays which are
            # (grand)children of this node. Therefore, we need to resize it.
            if len(child_node) == 0 and create_empty_structs:
                child_node.resize(1)

            for structure in child_node:
                yield from _iter_nodes_from_path(
                    structure,
                    path_parts[1:],
                    get_empty_arrays,
                    ggd_idx,
                    create_empty_structs,
                )
    else:
        yield from _iter_nodes_from_path(
            child_node, path_parts[1:], get_empty_arrays, ggd_idx, create_empty_structs
        )
