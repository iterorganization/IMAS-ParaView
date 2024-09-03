import imaspy
import pytest

from vtkggdtools.ids_util import get_arrays_from_ids
from vtkggdtools.tests.fill_ggd import fill_scalar_quantity, fill_vector_quantity


@pytest.fixture
def ids_fixture():
    # Create dummy IDS with 4 time steps
    ids = imaspy.IDSFactory().new("edge_sources")
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
    num_vertices = 6
    num_edges = 7
    num_faces = 2
    for ggd_idx in range(len(ids.time)):
        fill_scalar_quantity(input_scalars[ggd_idx], num_vertices, num_edges, num_faces)
        fill_vector_quantity(input_vectors[ggd_idx], num_vertices, num_edges, num_faces)


def test_get_arrays_from_ids_empty(ids_fixture):
    ids, _, _ = ids_fixture
    retrieved_scalar_arrays, retrieved_vector_arrays = get_arrays_from_ids(ids)
    assert retrieved_scalar_arrays == []
    assert retrieved_vector_arrays == []


def test_get_arrays_from_ids_all_timesteps(ids_fixture):
    ids, input_scalars, input_vectors = ids_fixture
    fill_quantities(ids, input_scalars, input_vectors)

    retrieved_scalar_arrays, retrieved_vector_arrays = get_arrays_from_ids(ids)
    assert len(retrieved_scalar_arrays) == len(ids.time)
    assert len(retrieved_vector_arrays) == len(ids.time)
    assert all(scalar in retrieved_scalar_arrays for scalar in input_scalars)
    assert all(vector in retrieved_vector_arrays for vector in input_vectors)


def test_get_arrays_from_ids_single_timesteps(ids_fixture):
    ids, input_scalars, input_vectors = ids_fixture
    fill_quantities(ids, input_scalars, input_vectors)

    ggd_idx = 1
    retrieved_scalar_arrays, retrieved_vector_arrays = get_arrays_from_ids(
        ids, ggd_idx=ggd_idx
    )
    assert len(retrieved_scalar_arrays) == 1
    assert len(retrieved_vector_arrays) == 1
    assert input_scalars[ggd_idx] == retrieved_scalar_arrays[0]
    assert input_vectors[ggd_idx] == retrieved_vector_arrays[0]


def test_get_arrays_from_ids_not_defined_timesteps(ids_fixture):
    ids, input_scalars, input_vectors = ids_fixture
    fill_quantities(ids, input_scalars, input_vectors)

    # Get timestep not in IDS
    ggd_idx = len(ids.time)
    retrieved_scalar_arrays, retrieved_vector_arrays = get_arrays_from_ids(
        ids, ggd_idx=ggd_idx
    )
    assert retrieved_scalar_arrays == []
    assert retrieved_vector_arrays == []
