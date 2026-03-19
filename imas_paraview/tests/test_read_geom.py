import pytest
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import (
    VTK_HEXAHEDRON,
    VTK_LINE,
    VTK_POLY_LINE,
    VTK_POLYGON,
    VTK_POLYHEDRON,
    VTK_PYRAMID,
    VTK_QUAD,
    VTK_TETRA,
    VTK_TRIANGLE,
    VTK_VERTEX,
    VTK_WEDGE,
)

from imas_paraview.io.read_geom import (
    _get_vtk_cell_type,
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


def test_convert_grid_subset_geometry_wall(ids_name, dummy_ids):
    grid_ggd = get_grid_ggd(dummy_ids)
    points = vtkPoints()
    fill_vtk_points(grid_ggd, 0, points, ids_name)
    grid = convert_grid_subset_geometry_to_unstructured_grid(grid_ggd, -1, points)
    assert grid.GetPoints() is not None


@pytest.mark.parametrize(
    "dimension, num_points, expected",
    [
        (0, 1, VTK_VERTEX),
        (1, 2, VTK_LINE),
        (1, 5, VTK_POLY_LINE),
        (2, 3, VTK_TRIANGLE),
        (2, 4, VTK_QUAD),
        (2, 7, VTK_POLYGON),
        (3, 4, VTK_TETRA),
        (3, 5, VTK_PYRAMID),
        (3, 6, VTK_WEDGE),
        (3, 7, VTK_POLYHEDRON),
        (3, 8, VTK_HEXAHEDRON),
    ],
)
def test_get_vtk_cell_type(dimension, num_points, expected):
    assert _get_vtk_cell_type(dimension, num_points) == expected
