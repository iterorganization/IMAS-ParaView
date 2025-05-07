import imas
import pytest
from imas.ids_path import IDSPath

from imas_paraview.ids_util import (
    create_name_recursive,
    get_arrays_from_ids,
    recursive_ggd_path_search,
)
from imas_paraview.tests.fill_ggd import fill_scalar_quantity, fill_vector_quantity


@pytest.fixture
def ids_fixture():
    """Creates a dummy IDS with 4 time steps"""
    ids = imas.IDSFactory().new("edge_sources")
    ids.time = [1.1, 2.2, 3.3, 4.4]
    num_timesteps = len(ids.time)
    ids.source.resize(1)
    ids.source[0].ggd.resize(num_timesteps)
    input_scalars = [None] * num_timesteps
    input_vectors = [None] * num_timesteps
    for ggd_idx in range(num_timesteps):
        ids.source[0].ggd[ggd_idx].ion.resize(1)
        ids.source[0].ggd[ggd_idx].ion[0].state.resize(1)
        input_scalars[ggd_idx] = ids.source[0].ggd[ggd_idx].ion[0].energy
        input_vectors[ggd_idx] = ids.source[0].ggd[ggd_idx].ion[0].momentum
    return ids, input_scalars, input_vectors


def fill_quantities(ids, input_scalars, input_vectors):
    """Fills the input scalar and vector GGD arrays with dummy values"""
    num_vertices = 6
    num_edges = 7
    num_faces = 2
    for ggd_idx in range(len(ids.time)):
        fill_scalar_quantity(input_scalars[ggd_idx], num_vertices, num_edges, num_faces)
        fill_vector_quantity(input_vectors[ggd_idx], num_vertices, num_edges, num_faces)


def test_get_arrays_from_ids_empty(ids_fixture):
    """Tests get_arrays_from_ids when IDS has no filled GGD arrays"""
    ids, _, _ = ids_fixture
    retrieved_scalar_arrays, retrieved_vector_arrays = get_arrays_from_ids(ids)
    assert retrieved_scalar_arrays == []
    assert retrieved_vector_arrays == []


def test_get_arrays_from_ids_all_timesteps(ids_fixture):
    """Tests get_arrays_from_ids when requesting all time steps in the IDS"""

    ids, input_scalars, input_vectors = ids_fixture
    fill_quantities(ids, input_scalars, input_vectors)

    retrieved_scalar_arrays, retrieved_vector_arrays = get_arrays_from_ids(ids)
    assert len(retrieved_scalar_arrays) == len(ids.time)
    assert len(retrieved_vector_arrays) == len(ids.time)
    assert all(scalar in retrieved_scalar_arrays for scalar in input_scalars)
    assert all(vector in retrieved_vector_arrays for vector in input_vectors)


def test_get_arrays_from_ids_single_timesteps(ids_fixture):
    """Tests get_arrays_from_ids when requesting only a single time step in the IDS"""

    ids, input_scalars, input_vectors = ids_fixture
    fill_quantities(ids, input_scalars, input_vectors)

    for ggd_idx in range(len(ids.time)):
        retrieved_scalar_arrays, retrieved_vector_arrays = get_arrays_from_ids(
            ids, ggd_idx=ggd_idx
        )
        assert len(retrieved_scalar_arrays) == 1
        assert len(retrieved_vector_arrays) == 1
        assert input_scalars[ggd_idx] == retrieved_scalar_arrays[0]
        assert input_vectors[ggd_idx] == retrieved_vector_arrays[0]


def test_get_arrays_from_ids_not_defined_timesteps(ids_fixture):
    """Tests get_arrays_from_ids when requesting a single time step that does not exist
    in the IDS"""
    ids, input_scalars, input_vectors = ids_fixture
    fill_quantities(ids, input_scalars, input_vectors)

    # Get timestep not in IDS
    ggd_idx = len(ids.time)
    retrieved_scalar_arrays, retrieved_vector_arrays = get_arrays_from_ids(
        ids, ggd_idx=ggd_idx
    )
    assert retrieved_scalar_arrays == []
    assert retrieved_vector_arrays == []


def test_recursive_ggd_path_search():
    """Test if recursive_ggd_path_search returns all GGD array paths"""
    ids = imas.IDSFactory(version="3.41.0").new("edge_sources")

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

    scalar_arrays = []
    vector_arrays = []
    recursive_ggd_path_search(ids.metadata, scalar_arrays, vector_arrays)
    assert set(scalar_arrays) == set(es_scalar_v3_41)
    assert set(vector_arrays) == set(es_vector_v3_41)


def test_create_name_recursive_ggd():
    ids = imas.IDSFactory(version="3.40.0").new("edge_sources")

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

    name = create_name_recursive(array)
    assert name == "Source (Charge exchange) Ion (Ar) State (Ar+14) Energy"


def test_create_name_recursive_profiles():
    ids = imas.IDSFactory(version="3.40.0").new("core_profiles")

    ids.profiles_1d.resize(1)
    phi_potential = ids.profiles_1d[0].phi_potential
    electrons_pressure = ids.profiles_1d[0].electrons.pressure
    ids.profiles_1d[0].ion.resize(2)
    ids.profiles_1d[0].ion[0].label = "D"
    temperature = ids.profiles_1d[0].ion[0].temperature
    ids.profiles_1d[0].ion[1].label = "T"
    density = ids.profiles_1d[0].ion[1].density

    assert create_name_recursive(phi_potential) == "Phi_potential"
    assert create_name_recursive(electrons_pressure) == "Electrons Pressure"
    assert create_name_recursive(temperature) == "Ion (D) Temperature"
    assert create_name_recursive(density) == "Ion (T) Density"
