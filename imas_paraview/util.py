import logging
from typing import Optional

import numpy as np
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData

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
        if (
            structure_reference in ["generic_grid_dynamic", "generic_grid_aos3_root"]
            and getattr(node, "lifecycle_status", None) != "obsolescent"
        ):
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


def get_grid_ggd(ids, time=0, parent_idx=0):
    """Finds and returns the grid_ggd within IDS at the time index ggd_idx. If the
    grid_ggd at ggd_idx time index does not exist, it tries to return the first
    grid_ggd. If this does not exist, it returns None.

    Args:
        ids: The IDS for which to return the grid_gdd.
        time: Time value for which to load the grid.
        parent_idx: Index for any non-time-dependent parent Array of Structures.
            For example ``description_ggd[parent_idx].grid_ggd[ggd_idx]`` in the wall
            IDS.

    Returns:
        The grid_ggd node found, or None if not found.
    """
    grid_path = get_ggd_grid_path(ids.metadata)
    if grid_path is None:
        return None

    node = ids
    for path in grid_path.split("/"):
        node = node[path]
        if node.metadata.ndim == 0:
            pass  # Current node is a structure
        elif node.metadata.coordinate1.is_time_coordinate:
            # Time dependent array of structure
            # Let IMAS-Python handle the time mode (homogeneous/heterogeneous):
            time_array = node.coordinates[0]
            # Load closest previous time index
            ggd_idx = np.searchsorted(time_array, time, side="right") - 1

            if 0 <= ggd_idx < len(node):
                node = node[ggd_idx]
            else:
                node = node[0]
                logger.warning(
                    "The GGD grid was not found at time index %d, so first "
                    "grid was loaded instead.",
                    ggd_idx,
                )
        else:
            node = node[parent_idx]

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


def points_to_vtkpoly(points, is_closed=False, is_filled=False):
    """Convert a list of 3D points to a vtkPoints and vtkCellArray, which are combined
    into a single VtkPolyData object. The expected format of the points is:
    [(x1,y1,z1),(x2,y2,z2),...].

    Args:
        points: A list of 3D points, containing the X,Y,Z-coordinates in tuples of the
            form (x,y,z).
        is_closed: Boolean flag whether to close the contour. If set to true, the first
            and last point are connected by a vtkLine.
        is_filled: Boolean flag whether to fill a closed contour. If True, the
            vtkPolyData will contain a single polygon, otherwise it will contain the
            outline of line segments. Only allowed when is_closed=True.

    Returns:
        vtkPolyData containing the vtkPoints and vtkLines.
    """
    if not is_closed and is_filled:
        raise ValueError("Converting points to filled polygons requires is_closed=True")

    points = np.asarray(points)
    n = len(points)

    vtk_points = vtkPoints()
    vtk_points.SetData(numpy_to_vtk(points))

    poly = vtkPolyData()
    poly.SetPoints(vtk_points)

    if is_filled:
        # Create a single polygon cell
        polys = vtkCellArray()
        polys.InsertNextCell(n)
        for i in range(n):
            polys.InsertCellPoint(i)
        poly.SetPolys(polys)
    else:
        # Create outline from line segments
        lines = vtkCellArray()
        n_lines = n - 1 + (1 if is_closed else 0)
        for i in range(n_lines):
            lines.InsertNextCell(2)
            lines.InsertCellPoint(i % n)
            lines.InsertCellPoint((i + 1) % n)
        poly.SetLines(lines)

    return poly
