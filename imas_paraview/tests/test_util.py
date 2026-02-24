import imas
import numpy as np
import pytest
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from imas_paraview.util import find_closest_indices, get_grid_ggd, points_to_vtkpoly


@pytest.fixture
def points():
    return [(0, 0, 0), (1, 0, 0), (1, 1, 0)]


def extract_indices(cell_array):
    """Extract all point indices as a list of tuples from vtkCellArray."""
    indices = []
    cell_array.InitTraversal()
    id_list = vtk.vtkIdList()
    while cell_array.GetNextCell(id_list):
        indices.append(tuple(id_list.GetId(i) for i in range(id_list.GetNumberOfIds())))
    return indices


@pytest.mark.parametrize(
    "request_times, expected_indices",
    [
        ([2.2, 3.3], [1, 2]),
        ([4.4, 5.5], [3, 4]),
        ([1.5, 3.8, 5.9], [0, 2, 4]),
        ([10], [4]),
        ([0.1, 0.5, 0.9, 1.01], []),
        ([0.5, 1.01, 2.5, 4.7], [1, 3]),
        ([2.2, 2.2, 4.4, 4.4], [1, 1, 3, 3]),
    ],
)
def test_find_closest_indices(request_times, expected_indices):
    time_array = np.array([1.1, 2.2, 3.3, 4.4, 5.5])
    time_indices = find_closest_indices(request_times, time_array)
    assert time_indices == expected_indices


def test_points_to_vtkpoly_line_open(points):
    poly = points_to_vtkpoly(points, is_closed=False, is_filled=False)
    assert poly.GetNumberOfPoints() == 3
    assert poly.GetNumberOfLines() == 2
    assert poly.GetNumberOfPolys() == 0
    assert np.array_equal(vtk_to_numpy(poly.GetPoints().GetData()), np.array(points))
    assert extract_indices(poly.GetLines()) == [(0, 1), (1, 2)]


def test_points_to_vtkpoly_line_closed(points):
    poly = points_to_vtkpoly(points, is_closed=True, is_filled=False)
    assert poly.GetNumberOfPoints() == 3
    assert poly.GetNumberOfLines() == 3
    assert poly.GetNumberOfPolys() == 0
    assert np.array_equal(vtk_to_numpy(poly.GetPoints().GetData()), np.array(points))
    assert extract_indices(poly.GetLines()) == [(0, 1), (1, 2), (2, 0)]


def test_points_to_vtkpoly_filled_polygon(points):
    poly = points_to_vtkpoly(points, is_closed=True, is_filled=True)

    np_points = vtk_to_numpy(poly.GetPoints().GetData())
    np.testing.assert_array_equal(np_points, np.array(points))
    assert poly.GetNumberOfPoints() == 3
    assert poly.GetNumberOfPolys() == 1
    assert poly.GetNumberOfLines() == 0
    assert np.array_equal(vtk_to_numpy(poly.GetPoints().GetData()), np.array(points))
    assert extract_indices(poly.GetPolys()) == [(0, 1, 2)]


def test_points_to_vtkpoly_filled_requires_closed(points):
    with pytest.raises(ValueError):
        points_to_vtkpoly(points, is_closed=False, is_filled=True)


def test_get_grid_ggd():
    dd_version = "4.0.0"
    ids = imas.IDSFactory(version=dd_version).new("edge_profiles")
    ids.time = [0.0, 1.1, 2.2]
    ids.ids_properties.homogeneous_time = imas.ids_defs.IDS_TIME_MODE_HOMOGENEOUS
    assert get_grid_ggd(ids) is None
    ids.grid_ggd.resize(3)
    # Test with id() otherwise IDSStructure will compare contents of structure
    assert id(get_grid_ggd(ids)) == id(ids.grid_ggd[0])
    assert id(get_grid_ggd(ids, time=0.5)) == id(ids.grid_ggd[0])
    assert id(get_grid_ggd(ids, time=-0.5)) == id(ids.grid_ggd[0])
    assert id(get_grid_ggd(ids, time=1.1)) == id(ids.grid_ggd[1])
    assert id(get_grid_ggd(ids, time=2.0)) == id(ids.grid_ggd[1])
    assert id(get_grid_ggd(ids, time=2.2)) == id(ids.grid_ggd[2])
    assert id(get_grid_ggd(ids, time=3.3)) == id(ids.grid_ggd[2])

    ids = imas.IDSFactory(version=dd_version).new("wall")
    ids.time = [0.0]
    ids.ids_properties.homogeneous_time = imas.ids_defs.IDS_TIME_MODE_HOMOGENEOUS
    assert get_grid_ggd(ids) is None
    ids.description_ggd.resize(2)
    ids.description_ggd[0].grid_ggd.resize(1)
    ids.description_ggd[1].grid_ggd.resize(1)
    assert id(get_grid_ggd(ids)) == id(ids.description_ggd[0].grid_ggd[0])
    assert id(get_grid_ggd(ids, parent_idx=1)) == id(ids.description_ggd[1].grid_ggd[0])
    assert get_grid_ggd(ids, parent_idx=2) is None

    ids = imas.IDSFactory(version=dd_version).new("equilibrium")
    ids.time = [0.0]
    ids.ids_properties.homogeneous_time = imas.ids_defs.IDS_TIME_MODE_HOMOGENEOUS
    assert get_grid_ggd(ids) is None
    ids.grids_ggd.resize(1)
    assert get_grid_ggd(ids) is None
    ids.grids_ggd[0].grid.resize(1)
    assert id(get_grid_ggd(ids)) == id(ids.grids_ggd[0].grid[0])

    ids = imas.IDSFactory(version=dd_version).new("camera_ir")  # IDS without GGD grid
    ids.time = [0.0]
    ids.ids_properties.homogeneous_time = imas.ids_defs.IDS_TIME_MODE_HOMOGENEOUS
    assert get_grid_ggd(ids) is None
