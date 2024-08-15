import imaspy

from vtkggdtools.ids_util import get_arrays_from_ids
from vtkggdtools.tests.fill_ggd import fill_scalar_quantity, fill_vector_quantity


def test_get_arrays_from_ids():
    ids = imaspy.IDSFactory().new("edge_sources")
    ids.source.resize(1)
    ids.source[0].ggd.resize(1)
    ids.source[0].ggd[0].ion.resize(1)
    ids.source[0].ggd[0].ion[0].state.resize(1)
    scalar_quantity1 = ids.source[0].ggd[0].ion[0].energy
    scalar_quantity2 = ids.source[0].ggd[0].ion[0].state[0].energy
    vector_quantity1 = ids.source[0].ggd[0].ion[0].momentum
    vector_quantity2 = ids.source[0].ggd[0].ion[0].state[0].momentum

    scalar_array_list, vector_array_list = get_arrays_from_ids(ids)
    assert scalar_array_list == []
    assert vector_array_list == []

    # Fill scalar and vector quantity
    num_vertices = 6
    num_edges = 7
    num_faces = 2
    fill_scalar_quantity(scalar_quantity1, num_vertices, num_edges, num_faces)
    fill_scalar_quantity(scalar_quantity2, num_vertices, num_edges, num_faces)
    fill_vector_quantity(vector_quantity1, num_vertices, num_edges, num_faces)
    fill_vector_quantity(vector_quantity2, num_vertices, num_edges, num_faces)

    scalar_array_list, vector_array_list = get_arrays_from_ids(ids)
    assert scalar_array_list == [scalar_quantity1, scalar_quantity2]
    assert vector_array_list == [vector_quantity1, vector_quantity2]
