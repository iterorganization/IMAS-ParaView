import logging

import numpy as np
from imas import identifiers
from imas.ids_structure import IDSStructure
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkStructuredGrid

from imas_paraview.util import pol_to_cart

logger = logging.getLogger("imas_paraview")


def convert_structured_grid_subset_to_vtk(grid_ggd, subset_idx) -> vtkStructuredGrid:
    """Convert a grid subset of a structured GGD grid to a vtkStructuredGrid.

    Args:
        grid_ggd: A grid_ggd IDS node containing a structured GGD grid.
        subset_idx: Index of the grid subset.

    Returns:
        A vtkStructuredGrid for the requested subset.
    """

    space0, space1 = grid_ggd.space[0], grid_ggd.space[1]
    coords0 = _extract_1d_coords(space0)
    coords1 = _extract_1d_coords(space1)
    coord0_type = _get_coord_type(space0)
    coord1_type = _get_coord_type(space1)
    if len(grid_ggd.grid_subset[subset_idx].element) != 0:
        subset = grid_ggd.grid_subset[subset_idx]
        coords0, coords1 = _build_subset(subset, coords0, coords1)

    return _build_sgrid(coords0, coords1, coord0_type, coord1_type)


def _build_sgrid(coords0, coords1, coord0_type, coord1_type) -> vtkStructuredGrid:
    """Build a vtkStructuredGrid from the 1D coordinate arrays.

    Args:
        coords0: The 1D coordinate array for the first space.
        coords1: The 1D coordinate array for the second space.
        coord0_type: Coordinate type for the first space.
        coord1_type: Coordinate type for the second space.

    Returns:
        The created vtkStructuredGrid.
    """
    grid0, grid1 = np.meshgrid(coords0, coords1)
    grids = {coord0_type: grid0, coord1_type: grid1}

    cid = identifiers.coordinate_identifier
    R, PHI, Z, X, Y = cid.r.index, cid.phi.index, cid.z.index, cid.x.index, cid.y.index

    x_arr = np.zeros_like(grid0)
    y_arr = np.zeros_like(grid0)
    z_arr = np.zeros_like(grid0)

    if R in grids and PHI in grids:
        x_arr, y_arr = pol_to_cart(grids[R], grids[PHI])
    elif any(k in grids for k in {X, Y, Z, R}):
        # R is mapped to X by default in R-Z plots
        x_arr = grids.get(X, grids.get(R, x_arr))
        y_arr = grids.get(Y, y_arr)
        z_arr = grids.get(Z, z_arr)
    else:
        logger.warning("No spatial coordinates found. Showing grid on X-Y plane.")
        x_arr, y_arr = grid0, grid1

    xyz = np.stack([x_arr.flat, y_arr.flat, z_arr.flat], axis=1)
    vtk_pts = vtkPoints()
    vtk_pts.SetData(numpy_to_vtk(xyz, deep=True))

    sgrid = vtkStructuredGrid()
    sgrid.SetDimensions(len(coords0), len(coords1), 1)
    sgrid.SetPoints(vtk_pts)
    return sgrid


def _build_subset(subset, coords0, coords1):
    """Extract the specific 1D coordinate arrays spanned by a grid subset.

    Args:
        subset: The grid subset node containing element coordinate references.
        coords0: The full 1D coordinate array for the first space.
        coords1: The full 1D coordinate array for the second space.

    Returns:
        Tuple containing the subset coordinate arrays for both spaces.
    """

    idx0, idx1 = set(), set()
    for element in subset.element:
        for obj in element.object:
            if obj.space == 1:
                idx0.add(obj.index - 1)
            else:
                idx1.add(obj.index - 1)

    subset_coords0 = coords0[sorted(idx0)]
    subset_coords1 = coords1[sorted(idx1)]

    return subset_coords0, subset_coords1


def _get_coord_type(space) -> int:
    """Return the coordinate identifier for a space."""
    coordinates_type = space.coordinates_type
    if isinstance(coordinates_type[0], IDSStructure):
        return int(coordinates_type[0].index)
    return int(coordinates_type[0])


def _extract_1d_coords(space) -> np.ndarray:
    """Extract the coordinate geometry values from a 1D space.

    Args:
        space: The space node containing objects_per_dimension.

    Returns:
        NumPy array containing the coordinate geometry values.
    """
    objects = space.objects_per_dimension[0].object
    return np.array([obj.geometry[0] for obj in objects])
