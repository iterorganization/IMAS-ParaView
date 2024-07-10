from vtkmodules.vtkCommonCore import vtkPoints

from vtkggdtools.io.read_geom import (
    convert_grid_subset_geometry_to_unstructured_grid,
    fill_vtk_points,
)
from vtkggdtools.io.read_ps import read_plasma_state
from vtkggdtools.util import FauxIndexMap, get_first_grid


def test_read_plasma_state(ids_name, dummy_ids):
    """Tests reading plasma state from the IDS."""

    ids = dummy_ids
    subset_idx = 0
    space_idx = 0
    grid_ggd = get_first_grid(ids)
    points = vtkPoints()
    _aos_index_values = FauxIndexMap()  # FIX
    fill_vtk_points(grid_ggd, space_idx, points, ids_name)
    ugrid = convert_grid_subset_geometry_to_unstructured_grid(
        grid_ggd, subset_idx, points
    )
    read_plasma_state(ids_name, ids, _aos_index_values, subset_idx, ugrid)
