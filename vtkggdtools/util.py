from typing import Optional


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


def get_first_grid(ids):
    grid_path = get_ggd_grid_path(ids.metadata)
    if grid_path is None:
        return None

    node = ids
    for path in grid_path.split("/"):
        node = node[path]
        try:
            node = node[0]  # get first element of AoS
        except LookupError:
            pass  # apparently this was not an AoS :)

    return node
