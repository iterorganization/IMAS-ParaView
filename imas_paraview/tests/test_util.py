import imas
import numpy as np
import pytest
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from imas_paraview.util import (
    angles_to_vectors,
    create_vtk_spheres,
    ensure_unique_name,
    find_closest_indices,
    get_grid_ggd,
    points_to_vtkpoly,
    vel_pol_to_cart,
)


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


def test_vel_pol_to_cart():
    v_r = np.array([1.0, 1.0])
    v_phi = np.array([0.0, 0.0])
    v_z = np.array([2.0, 0.0])
    phi = np.array([0.0, np.pi / 2])

    result = vel_pol_to_cart(v_r, v_phi, v_z, phi)
    expected = np.array([[1.0, 0.0, 2.0], [0.0, 1.0, 0.0]])
    assert np.allclose(result, expected)


def test_create_vtk_spheres():
    p1 = [1.0, 2.0, 3.0]
    p2 = [5.0, 4.0, 3.0]
    r1 = 1.2
    r2 = 3.4
    positions = np.array([p1, p2])
    radii = np.array([r1, r2])
    spheres = create_vtk_spheres(positions, radii)
    points = vtk_to_numpy(spheres.GetPoints().GetData())

    n_per_sphere = points.shape[0] // 2
    d1 = np.linalg.norm(points[:n_per_sphere] - p1, axis=1)
    d2 = np.linalg.norm(points[n_per_sphere:] - p2, axis=1)

    assert np.allclose(d1, r1)
    assert np.allclose(d2, r2)


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


def test_ensure_unique_name():
    existing = ["coil1", "coil1 #1", "coil2", "coil3"]
    assert ensure_unique_name("coil4", existing) == "coil4"
    assert ensure_unique_name("coil2", existing) == "coil2 #1"
    assert ensure_unique_name("coil1", existing) == "coil1 #2"


def test_angles_to_vectors():
    r, phi, z = 1.0, 0.0, 0.0
    pol, tor = 0.0, 0.0
    position, direction = angles_to_vectors(r, phi, z, pol, tor)
    assert np.allclose(position, np.array([[1.0, 0.0, 0.0]]))
    assert np.allclose(direction, np.array([[1.0, 0.0, 0.0]]))

    r, phi, z = 2.0, np.pi / 2, 1.0
    pol, tor = 0.0, 0.0
    position, direction = angles_to_vectors(r, phi, z, pol, tor)
    assert np.allclose(position, np.array([[0.0, 2.0, 1.0]]))
    assert np.allclose(direction, np.array([[0.0, 1.0, 0.0]]))

    r, phi, z = 1.0, 0.0, 0.0
    pol, tor = np.pi / 2, 0.0
    position, direction = angles_to_vectors(r, phi, z, pol, tor)
    assert np.allclose(position, np.array([[1.0, 0.0, 0.0]]))
    assert np.allclose(direction, np.array([[0.0, 0.0, -1.0]]))

    r, phi, z = 5.0, np.pi, -1.0
    pol, tor = np.pi / 2, 0.0
    position, direction = angles_to_vectors(r, phi, z, pol, tor)
    assert np.allclose(position, np.array([[-5.0, 0.0, -1.0]]))
    assert np.allclose(direction, np.array([[0.0, 0.0, -1.0]]))

    r, phi, z = 3.0, np.pi / 2, 1.5
    pol, tor = 0.0, np.pi / 2
    position, direction = angles_to_vectors(r, phi, z, pol, tor)
    assert np.allclose(position, np.array([[0.0, 3.0, 1.5]]))
    assert np.allclose(direction, np.array([[-1.0, 0.0, 0.0]]))

    r, phi, z = 4.0, 3 * np.pi / 4, 5.0
    pol, tor = np.pi / 4, np.pi / 4
    position, direction = angles_to_vectors(r, phi, z, pol, tor)
    assert np.allclose(position, np.array([[-2 * np.sqrt(2), 2 * np.sqrt(2), 5.0]]))
    assert np.allclose(direction, np.array([[-np.sqrt(2) / 2, 0.0, -np.sqrt(2) / 2]])
