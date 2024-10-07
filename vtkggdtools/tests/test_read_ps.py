import imaspy
import pytest
from imaspy.ids_path import IDSPath
from vtkmodules.vtkCommonCore import vtkPoints

from vtkggdtools.io.read_geom import (
    convert_grid_subset_geometry_to_unstructured_grid,
    fill_vtk_points,
)
from vtkggdtools.io.read_ps import PlasmaStateReader
from vtkggdtools.tests.fill_ggd import fill_scalar_quantity, fill_vector_quantity
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
    assert name == "Source (Charge exchange) Ion (Ar) State (Ar+14) Energy [W.m^-3]"


def test_load_paths_from_ids_empty():
    """Test if load_paths_from_ids returns empty lists if no GGDs are filled"""
    ids = imaspy.IDSFactory(version="3.41.0").new("edge_sources")
    ps_reader = PlasmaStateReader(ids)

    _, _, scalar_array_paths, vector_array_paths = ps_reader.load_paths_from_ids()
    assert scalar_array_paths == []
    assert vector_array_paths == []


def test_load_paths_from_ids_filled():
    """Test if load_paths_from_ids returns only filled GGD paths"""
    ids = imaspy.IDSFactory(version="3.41.0").new("edge_sources")
    ps_reader = PlasmaStateReader(ids)

    ids.source.resize(1)
    ids.source[0].ggd.resize(1)
    ids.source[0].ggd[0].ion.resize(1)
    ids.source[0].ggd[0].ion[0].state.resize(1)
    energy = ids.source[0].ggd[0].ion[0].state[0].energy
    momentum = ids.source[0].ggd[0].ion[0].state[0].momentum
    _, _, scalar_array_paths, vector_array_paths = ps_reader.load_paths_from_ids()
    assert scalar_array_paths == []
    assert vector_array_paths == []

    fill_scalar_quantity(energy, 1, 1, 1)
    fill_vector_quantity(momentum, 1, 1, 1)
    _, _, scalar_array_paths, vector_array_paths = ps_reader.load_paths_from_ids()
    assert scalar_array_paths == [energy.metadata.path]
    assert vector_array_paths == [momentum.metadata.path]


def test_load_paths_from_ids_empty_first():
    """Test if load_paths_from_ids handles first entry being empty in IDSStructure"""
    ids = imaspy.IDSFactory(version="3.41.0").new("edge_sources")
    ps_reader = PlasmaStateReader(ids)

    ids.source.resize(1)
    ids.source[0].ggd.resize(1)
    ids.source[0].ggd[0].ion.resize(2)
    ids.source[0].ggd[0].ion[1].state.resize(1)
    energy = ids.source[0].ggd[0].ion[1].state[0].energy
    momentum = ids.source[0].ggd[0].ion[1].state[0].momentum

    _, _, scalar_array_paths, vector_array_paths = ps_reader.load_paths_from_ids()
    assert scalar_array_paths == []
    assert vector_array_paths == []

    fill_scalar_quantity(energy, 1, 1, 1)
    fill_scalar_quantity(momentum, 1, 1, 1)
    _, _, scalar_array_paths, vector_array_paths = ps_reader.load_paths_from_ids()
    assert scalar_array_paths == [energy.metadata.path]
    assert vector_array_paths == [momentum.metadata.path]


def test_load_paths_from_ids_all():
    """Test if load_paths_from_ids returns all GGD paths"""
    ids = imaspy.IDSFactory(version="3.41.0").new("edge_sources")
    ps_reader = PlasmaStateReader(ids)

    # Only load filled GGD array paths
    scalar_array_paths, vector_array_paths, _, _ = ps_reader.load_paths_from_ids(
        return_empty=True
    )

    es_scalar_v3_41 = [
        IDSPath("source/ggd/electrons/particles"),
        IDSPath("source/ggd/electrons/energy"),
        IDSPath("source/ggd/ion/energy"),
        IDSPath("source/ggd/ion/particles"),
        IDSPath("source/ggd/ion/state/particles"),
        IDSPath("source/ggd/ion/state/energy"),
        IDSPath("source/ggd/neutral/particles"),
        IDSPath("source/ggd/neutral/energy"),
        IDSPath("source/ggd/neutral/state/particles"),
        IDSPath("source/ggd/neutral/state/energy"),
        IDSPath("source/ggd/total_ion_energy"),
        IDSPath("source/ggd/current"),
    ]
    es_vector_v3_41 = [
        IDSPath("source/ggd/ion/momentum"),
        IDSPath("source/ggd/ion/state/momentum"),
        IDSPath("source/ggd/neutral/momentum"),
        IDSPath("source/ggd/neutral/state/momentum"),
        IDSPath("source/ggd/momentum"),
    ]

    assert set(scalar_array_paths) == set(es_scalar_v3_41)
    assert set(vector_array_paths) == set(es_vector_v3_41)


def test_load_paths_from_ids_multiple_timesteps():
    """Test if load_paths_from_ids when GGD contains multiple time steps"""
    ids = imaspy.IDSFactory(version="3.41.0").new("edge_sources")
    ps_reader = PlasmaStateReader(ids)

    ids.source.resize(1)
    ids.source[0].ggd.resize(2)
    ids.source[0].ggd[0].ion.resize(1)
    ids.source[0].ggd[0].ion[0].state.resize(1)
    energy = ids.source[0].ggd[0].ion[0].state[0].energy
    ids.source[0].ggd[1].ion.resize(1)
    ids.source[0].ggd[1].ion[0].state.resize(1)
    momentum = ids.source[0].ggd[1].ion[0].state[0].momentum
    _, _, scalar_array_paths, vector_array_paths = ps_reader.load_paths_from_ids(
        ggd_idx=0
    )
    assert scalar_array_paths == []
    assert vector_array_paths == []

    _, _, scalar_array_paths, vector_array_paths = ps_reader.load_paths_from_ids(
        ggd_idx=1
    )
    assert scalar_array_paths == []
    assert vector_array_paths == []

    fill_scalar_quantity(energy, 1, 1, 1)
    fill_vector_quantity(momentum, 1, 1, 1)
    _, _, scalar_array_paths, vector_array_paths = ps_reader.load_paths_from_ids(
        ggd_idx=0
    )
    assert scalar_array_paths == [energy.metadata.path]
    assert vector_array_paths == []

    _, _, scalar_array_paths, vector_array_paths = ps_reader.load_paths_from_ids(
        ggd_idx=1
    )
    assert scalar_array_paths == []
    assert vector_array_paths == [momentum.metadata.path]
