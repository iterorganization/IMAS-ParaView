import pytest
from vtkmodules.vtkCommonCore import vtkPoints

from vtkggdtools.io.read_geom import (
    convert_grid_subset_geometry_to_unstructured_grid,
    fill_vtk_points,
)
from vtkggdtools.io.read_ps import read_plasma_state
from vtkggdtools.util import get_first_grid


def test_read_plasma_state(ids_name, dummy_ids):
    """Tests reading plasma state from the IDS."""
    if ids_name == "transport_solver_numerics":
        pytest.skip(reason="boundary_conditions_ggd are not filled yet")
    elif ids_name == "runaway_electrons":
        pytest.skip(reason="runaway_electrons is not implemented")

    ids = dummy_ids
    subset_idx = 0
    space_idx = 0
    grid_ggd = get_first_grid(ids)
    points = vtkPoints()
    fill_vtk_points(grid_ggd, space_idx, points, ids_name)
    ugrid = convert_grid_subset_geometry_to_unstructured_grid(
        grid_ggd, subset_idx, points
    )
    read_plasma_state(ids, subset_idx, ugrid)
