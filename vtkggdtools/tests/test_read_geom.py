import imaspy
import pytest
from vtkmodules.vtkCommonCore import vtkPoints

from vtkggdtools.io.read_geom import (
    convert_grid_subset_geometry_to_unstructured_grid,
    fill_vtk_points,
)
from vtkggdtools.tests.fill_ggd import fill_ids
from vtkggdtools.util import get_first_grid


@pytest.fixture(
    params=[
        "distribution_sources",
        "edge_profiles",
        "edge_transport",
        "mhd",
        "runaway_electrons",
        "transport_solver_numerics",
        "waves",
        "distributions",
        "edge_sources",
        "equilibrium",
        "radiation",
        "tf",
        "wall",
    ]
)
def ids_name(request):
    """Provides IDS names for testing."""

    return request.param


@pytest.fixture
def dummy_ids(ids_name):
    """Creates a dummy IDS object with a dummy grid and random GGD values for
    testing purposes."""

    ids = imaspy.IDSFactory().new(ids_name)
    fill_ids(ids)
    return ids


def test_validate_dummy_ids(dummy_ids):
    """Validates the dummy IDS object created by the fixture."""

    ids = dummy_ids
    ids.validate()


def test_fill_vtk_points(ids_name, dummy_ids):
    """Tests filling VTK points from the IDS grid."""

    ids = dummy_ids
    space_idx = 0
    points = vtkPoints()
    grid_ggd = get_first_grid(ids)
    fill_vtk_points(grid_ggd, space_idx, points, ids_name)
    assert points.GetNumberOfPoints() > 0


def test_convert_grid_subset_geometry_to_unstructured_grid(ids_name, dummy_ids):
    """Tests grid subset geometry conversion to unstructured grid."""

    ids = dummy_ids
    subset_idx = 0
    space_idx = 0
    grid_ggd = get_first_grid(ids)
    points = vtkPoints()
    fill_vtk_points(grid_ggd, space_idx, points, ids_name)
    grid = convert_grid_subset_geometry_to_unstructured_grid(
        grid_ggd, subset_idx, points
    )
    assert grid.GetPoints() is not None
