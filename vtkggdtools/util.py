import logging
from typing import Optional

import numpy as np

logger = logging.getLogger("vtkggdtools")

SUPPORTED_IDS_NAMES = [
    "edge_profiles",
    "edge_sources",
    "edge_transport",
    "equilibrium",
    "mhd",
    "radiation",
    "runaway_electrons",
    "wall",
]

# FIXME: Some IDSs do not have the structure as described by the GGD guidelines. e.g.
# waves has its grid defined under coherent_wave/full_wave/grid instead of in a
# dedicated grid_ggd. They are for now denoted as experimental IDS, as they are
# currently not covered by the unit tests.
EXPERIMENTAL_IDS_NAMES = [
    "distribution_sources",
    "distributions",
    "tf",
    "transport_solver_numerics",
    "waves",
]


class FauxIndexMap:
    def __getitem__(self, item):
        return 0

    def get(self, name, default=None):
        return 0


def format_units(node) -> str:
    """Return the unit of the node surrounded by square brackets."""
    return f"[{node.metadata.units}]"


def iter_metadata_tree(meta):
    for child in meta._children.values():
        yield child
        yield from iter_metadata_tree(child)


def get_ggd_grid_path(ids_metadata) -> Optional[str]:
    """Finds the path of the GGD grid node within IDS metadata.

    Args:
        ids_metadata: The metadata of an IDS.

    Returns:
        The path string of the GGD grid node, or None if not found.
    """
    # Find the DD node defining the GGD grid:
    for node in iter_metadata_tree(ids_metadata):
        structure_reference = getattr(node, "structure_reference", None)
        if structure_reference in ["generic_grid_dynamic", "generic_grid_aos3_root"]:
            if getattr(node, "lifecycle_status", None) != "obsolescent":
                return node.path_string
    return None  # There are no GGD grids inside this IDS


def get_ggd_path(ids_metadata) -> Optional[str]:
    """Finds the path of the GGD node within IDS metadata.

    Args:
        ids_metadata: The metadata of an IDS.

    Returns:
        The path string of the GGD node, or None if not found.
    """
    for node in iter_metadata_tree(ids_metadata):
        metadata_name = getattr(node, "name", None)

        # Check for ggd in name metadata
        if metadata_name == "ggd":
            return node.path_string
    return None


def get_grid_ggd(ids, ggd_idx=0):
    """Finds and returns the first grid_ggd within IDS.

    Args:
        ids: The IDS for which to return the first grid_gdd.
        ggd_idx: Time index for which to load the grid. Defaults to 0.

    Returns:
        The first grid_ggd node found, or None if not found.
    """
    grid_path = get_ggd_grid_path(ids.metadata)
    if grid_path is None:
        return None

    node = ids
    for path in grid_path.split("/"):

        try:
            node = node[path]
        except ValueError:
            logger.warning(
                "Could not find a valid grid_ggd to load, because node "
                f"{node.metadata.name} does not have a {path}."
            )
            return None

        try:
            node = node[ggd_idx]
        except (LookupError, ValueError):
            # if node at ggd_idx does not exist, instead try at index 0
            try:
                node = node[0]
                logger.warning(
                    f"The GGD grid was not found at time index {ggd_idx}, so first "
                    "grid was loaded instead."
                )
            except (LookupError, ValueError):
                pass  # apparently this was not an AoS :)

    return node


def create_first_grid(ids):
    """Creates and returns the first grid_ggd within IDS.

    Args:
        ids: The IDS for which to create and return the first grid_gdd.

    Returns:
        The created first grid_gdd, or None if grid_path is None.
    """
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
    """Creates and returns the first GGD within IDS.

    Args:
        ids: The IDS for which to create and return the first GGD.

    Returns:
        The created first GGD, or None if ggd_path is None.
    """
    ggd_path = get_ggd_path(ids.metadata)
    if ggd_path is None:
        return None

    node = ids
    for path in ggd_path.split("/"):
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


def int32array(int_list):
    """Converts a list of ints to a numpy array containing values of type int32.

    Args:
        int_list: List of int64s to be converted

    Returns:
        List of converted int32s
    """
    return np.array(int_list, dtype=np.int32)
