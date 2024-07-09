import imaspy
import pytest
from vtkmodules.vtkCommonCore import vtkPoints

from vtkggdtools.io.read_geom import fill_vtk_points
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
        "wave",
        "distributions",
        "edge_sources",
        "equilibrium",
        "radiation",
        "tf",
        "wall",
    ]
)
def ids_name(request):
    return request.param


def test_fill_vtk_points(ids_name):
    factory = imaspy.IDSFactory()
    method = getattr(factory, ids_name, None)

    if method is not None:
        ids = method()
        fill_ids(ids)
        space_idx = 0
        points = vtkPoints()
        grid_ggd = get_first_grid(ids)
        fill_vtk_points(grid_ggd, space_idx, points, ids_name)
