from typing import Optional

import numpy as np


def iter_metadata_tree(meta):
    for child in meta._children.values():
        yield child
        yield from iter_metadata_tree(child)


def get_ggd_grid_path(ids_metadata) -> Optional[str]:
    # Find the DD node defining the GGD grid:
    for node in iter_metadata_tree(ids_metadata):
        structure_reference = getattr(node, "structure_reference", None)
        if structure_reference in ["generic_grid_dynamic", "generic_grid_aos3_root"]:
            if getattr(node, "lifecycle_status", None) != "obsolescent":
                return node.path_string
    return None  # There are no GGD grids inside this IDS


def get_ggd_path(ids_metadata) -> Optional[str]:
    for node in iter_metadata_tree(ids_metadata):
        metadata_name = getattr(node, "name", None)

        # Check for ggd in name metadata
        if metadata_name == "ggd":
            return node.path_string
    return None


def get_first_grid(ids):
    grid_path = get_ggd_grid_path(ids.metadata)
    if grid_path is None:
        return None

    node = ids
    for path in grid_path.split("/"):
        node = node[path]
        try:
            node = node[0]  # get first element of AoS
        except (LookupError, ValueError):
            pass  # apparently this was not an AoS :)

    return node


def create_first_grid(ids):
    grid_path = get_ggd_grid_path(ids.metadata)
    if grid_path is None:
        return None

    node = ids
    for path in grid_path.split("/"):
        node = node[path]
        try:
            if len(node) == 0:
                node.resize(1)
                node = node[0]
        except TypeError:
            pass  # This was a structure and not an AoS

    return node


def create_first_ggd(ids):
    grid_path = get_ggd_path(ids.metadata)
    if grid_path is None:
        return None

    node = ids
    for path in grid_path.split("/"):
        node = node[path]
        try:
            if len(node) == 0:
                node.resize(1)
                node = node[0]
            # Node was already resized
            elif len(node) == 1:
                node = node[0]
        except TypeError:
            pass  # This was a structure and not an AoS

    return node


def int64_to_int32(int_list):
    """
    Converts a list of type int64 to a numpy array containing values of
    type int32.
    """
    int32_array = np.array(int_list, dtype=np.int32)
    return int32_array
