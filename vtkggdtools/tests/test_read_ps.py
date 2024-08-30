import imaspy
import pytest
from vtkmodules.vtkCommonCore import vtkPoints

from vtkggdtools.io.read_geom import (
    convert_grid_subset_geometry_to_unstructured_grid,
    fill_vtk_points,
)
from vtkggdtools.io.read_ps import PlasmaStateReader
from vtkggdtools.util import get_grid_ggd


def test_read_plasma_state(ids_name, dummy_ids):
    """Tests reading plasma state from the IDS."""
    if ids_name == "transport_solver_numerics":
        pytest.skip(reason="boundary_conditions_ggd are not filled yet")
    elif ids_name == "runaway_electrons":
        pytest.skip(reason="runaway_electrons is not implemented")

    subset_idx = 0
    space_idx = 0
    grid_ggd = get_grid_ggd(dummy_ids)
    points = vtkPoints()
    fill_vtk_points(grid_ggd, space_idx, points, ids_name)
    ugrid = convert_grid_subset_geometry_to_unstructured_grid(
        grid_ggd, subset_idx, points
    )
    ps_reader = PlasmaStateReader(dummy_ids)
    ps_reader.read_plasma_state(subset_idx, ugrid)


def test_create_name_with_units():
    ids = imaspy.IDSFactory(version="3.40.0").new("edge_sources")
    ps_reader = PlasmaStateReader(ids)

    # Resize relevant structures
    ids.source.resize(1)
    ids.source[0].ggd.resize(1)
    ids.source[0].ggd[0].ion.resize(1)
    ids.source[0].ggd[0].ion[0].state.resize(1)
    array = ids.source[0].ggd[0].ion[0].state[0].energy

    # Set identifier name and labels
    ids.source[0].identifier.name = "Charge exchange"
    ids.source[0].ggd[0].ion[0].label = "Ar"
    ids.source[0].ggd[0].ion[0].state[0].label = "Ar+14"

    name = ps_reader._create_name_with_units(array)
    assert name == "source (Charge exchange) ion (Ar) state (Ar+14) energy [W.m^-3]"
