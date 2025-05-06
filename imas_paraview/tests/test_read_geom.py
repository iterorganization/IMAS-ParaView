from vtkmodules.vtkCommonCore import vtkPoints

from imas_paraview.io.read_geom import (
    convert_grid_subset_geometry_to_unstructured_grid,
    fill_vtk_points,
)
from imas_paraview.util import get_grid_ggd


def test_fill_vtk_points(ids_name, dummy_ids):
    """Tests filling VTK points from the IDS grid."""

    space_idx = 0
    points = vtkPoints()
    grid_ggd = get_grid_ggd(dummy_ids)
    fill_vtk_points(grid_ggd, space_idx, points, ids_name)
    assert points.GetNumberOfPoints() > 0


def test_convert_grid_subset_geometry_to_unstructured_grid(ids_name, dummy_ids):
    """Tests grid subset geometry conversion to unstructured grid."""

    subset_idx = 0
    space_idx = 0
    grid_ggd = get_grid_ggd(dummy_ids)
    points = vtkPoints()
    fill_vtk_points(grid_ggd, space_idx, points, ids_name)
    grid = convert_grid_subset_geometry_to_unstructured_grid(
        grid_ggd, subset_idx, points
    )
    assert grid.GetPoints() is not None
