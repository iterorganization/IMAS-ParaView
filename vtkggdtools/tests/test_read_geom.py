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
    """Provides IDS names for testing.

    Args:
        request (FixtureRequest): The pytest request object.

    Returns:
        str: The name of the IDS to test with.
    """
    return request.param


@pytest.fixture
def create_dummy_ids(ids_name):
    """Creates a dummy IDS object with a dummy grid and random GGD values for
    testing purposes.

    Args:
        ids_name (str): The name of the IDS to create.

    Returns:
        The IDS object, or None if the IDS name is not valid.
    """
    factory = imaspy.IDSFactory()
    method = getattr(factory, ids_name, None)

    if method is not None:
        ids = method()
        fill_ids(ids)
        return ids
    else:
        return None


def test_validate_dummy_ids(create_dummy_ids):
    """Validates the dummy IDS object created by the fixture.

    Args:
        create_dummy_ids (fixture): The fixture to create the dummy IDS.

    Asserts:
        The created IDS object is valid.
    """
    ids = create_dummy_ids
    ids.validate()


def test_fill_vtk_points(ids_name, create_dummy_ids):
    """Tests filling VTK points from the IDS grid.

    Args:
        ids_name (str): The name of the IDS.
        create_dummy_ids (fixture): The fixture to create the dummy IDS.

    Asserts:
        The number of points in the VTK points object is greater than 0.
    """
    ids = create_dummy_ids
    space_idx = 0
    points = vtkPoints()
    grid_ggd = get_first_grid(ids)
    fill_vtk_points(grid_ggd, space_idx, points, ids_name)
    assert points.GetNumberOfPoints() > 0


def test_convert_grid_subset_geometry_to_unstructured_grid(ids_name, create_dummy_ids):
    """Tests grid subset geometry conversion to unstructured grid.

    Args:
        ids_name (str): The name of the IDS.
        create_dummy_ids (fixture): The fixture to create the dummy IDS.

    Asserts:
        The resulting grid has points.
    """
    ids = create_dummy_ids
    subset_idx = 0
    space_idx = 0
    grid_ggd = get_first_grid(ids)
    points = vtkPoints()
    fill_vtk_points(grid_ggd, space_idx, points, ids_name)
    grid = convert_grid_subset_geometry_to_unstructured_grid(
        grid_ggd, subset_idx, points
    )
    assert grid.GetPoints() is not None
