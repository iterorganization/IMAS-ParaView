import logging
from typing import Optional

import numpy as np
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkLine, vtkPolyData

logger = logging.getLogger("imas_paraview")


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
    """Finds and returns the grid_ggd within IDS at the time index ggd_idx. If the
    grid_ggd at ggd_idx time index does not exist, it tries to return the first
    grid_ggd. If this does not exist, it returns None.

    Args:
        ids: The IDS for which to return the grid_gdd.
        ggd_idx: Time index for which to load the grid.

    Returns:
        The grid_ggd node found, or None if not found.
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


def find_closest_indices(values_to_extract, source_array):
    """Find indices of the closest values in source_array that are less than or equal
    to each value in values_to_extract.

    Args:
        values_to_extract: Values to find in the source array.
        source_array: Array to search for closest values.

    Returns:
        Indices of closest values in source_array for each value in values_to_extract.
    """
    indices = np.searchsorted(source_array, values_to_extract, side="right") - 1
    return list(indices[indices >= 0])


def pol_to_cart(rho, phi):
    """Convert from polar (or cylindrical) coordinates to cartesian.

    Args:
        rho: the radial distance
        phi: the azimuth angle in radians

    Returns:
        Tuple containing the x and y coordinates
    """
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x, y)


def points_to_vtkpoly(points, is_closed=False):
    """Convert a list of 3D points to a vtkPoints and vtkCellArray, which are combined
    into a single VtkPolyData object. The expected format of the points is:
    [(x1,y1,z1),(x2,y2,z2),...].

    Args:
        points: A list of 3D points, containing the X,Y,Z-coordinates in tuples of the
            form (x,y,z).
        is_closed: Boolean flag whether to close the contour. If set to true, the first
            and last point are connected by a vtkLine.

    Returns:
        vtkPolyData containing the vtkPoints and vtkLines.
    """
    vtk_points = vtkPoints()
    vtk_lines = vtkCellArray()
    for i, point in enumerate(points):
        vtk_points.InsertNextPoint(*point)
        line = vtkLine()

        if i != len(points) - 1:
            line.GetPointIds().SetId(0, i)
            line.GetPointIds().SetId(1, i + 1)
            vtk_lines.InsertNextCell(line)

    # Close loop if the points are closed
    if is_closed:
        line = vtkLine()
        line.GetPointIds().SetId(0, len(points) - 1)
        line.GetPointIds().SetId(1, 0)
        vtk_lines.InsertNextCell(line)

    vtk_poly = vtkPolyData()
    vtk_poly.SetPoints(vtk_points)
    vtk_poly.SetLines(vtk_lines)
    return vtk_poly
