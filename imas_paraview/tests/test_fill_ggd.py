import imas
import numpy as np

from imas_paraview.tests.fill_ggd import fill_ids


def test_validate_dummy_ids(dummy_ids):
    """Validates the dummy IDS object created by the fixture."""

    ids = dummy_ids
    ids.validate()


def test_fill_ggd():
    ids_name = "edge_profiles"
    dd_version = "4.0.0"
    ids = imas.IDSFactory(version=dd_version).new(ids_name)
    fill_ids(ids)

    space = ids.grid_ggd[0].space[0]

    expected_vertices = [
        [0.0, 0.0],
        [0.5, 0.0],
        [0.0, 1.0],
        [0.5, 1.0],
    ]
    for i, expected in enumerate(expected_vertices):
        assert np.allclose(space.objects_per_dimension[0].object[i].geometry, expected)

    expected_edges = [
        [1, 2],
        [3, 4],
        [1, 3],
        [2, 4],
    ]
    for i, expected in enumerate(expected_edges):
        assert np.array_equal(space.objects_per_dimension[1].object[i].nodes, expected)

    assert np.array_equal(space.objects_per_dimension[2].object[0].nodes, [1, 2, 4, 3])
