import imas
import numpy as np
import pytest
from conftest import DD_VERSION

from imas_paraview.plugins.coils_non_axisymmetric import CoilsNonAxisymmetricReader


@pytest.fixture
def elements():
    """Returns a elements IDS node of a coils_non_axisymmetric IDS"""
    ids = imas.IDSFactory(version=DD_VERSION).new("coils_non_axisymmetric")
    ids.coil.resize(1)
    ids.coil[0].conductor.resize(1)
    return ids.coil[0].conductor[0].elements


@pytest.fixture
def reader():
    return CoilsNonAxisymmetricReader()


@pytest.fixture
def base_points():
    """Returns a standard valid arc points setup."""
    p_start = np.array([5.0, 0.0, 0.0])
    p_intermediate = np.array([0.0, 5.0, 0.0])
    p_end = np.array([-5.0, 0.0, 0.0])
    p_centre = np.array([0.0, 0.0, 0.0])
    return p_start, p_intermediate, p_end, p_centre


def test_create_line_segment(elements, reader):
    elements.start_points.r = [10.0]
    elements.start_points.phi = [0.0]
    elements.start_points.z = [5.0]

    elements.end_points.r = [10.0]
    elements.end_points.phi = [np.pi / 2]
    elements.end_points.z = [3.0]

    points = reader._create_line_segment(elements, 0)

    assert len(points) == 2
    assert np.allclose(points[0], [10.0, 0.0, 5.0])
    assert np.allclose(points[1], [0.0, 10.0, 3.0])


def test_arc_of_circle_horizontal(elements, reader):
    elements.start_points.r = [10.0]
    elements.start_points.phi = [0.0]
    elements.start_points.z = [5.0]

    elements.intermediate_points.r = [10.0]
    elements.intermediate_points.phi = [np.pi / 2]
    elements.intermediate_points.z = [5.0]

    elements.end_points.r = [10.0]
    elements.end_points.phi = [np.pi]
    elements.end_points.z = [5.0]

    elements.centres.r = [0.0]
    elements.centres.phi = [0.0]
    elements.centres.z = [5.0]

    resolution = 10
    reader.resolution = resolution
    points = reader._create_circular_geometry(elements, 0, is_full_circle=False)

    assert np.allclose(points[0], [10.0, 0.0, 5.0])
    assert np.allclose(points[-1], [-10.0, 0.0, 5.0])
    assert np.all(points[:, 2] == 5.0)
    center = np.array([0.0, 0.0, 5.0])
    for p in points:  # All points should lie on a circle
        assert np.linalg.norm(np.array(p) - center) == pytest.approx(10.0)


def test_arc_of_circle_vertical(elements, reader):
    elements.start_points.r = [0.0]
    elements.start_points.phi = [0.0]
    elements.start_points.z = [-6.0]

    elements.intermediate_points.r = [5.0]
    elements.intermediate_points.phi = [0.0]
    elements.intermediate_points.z = [-1.0]

    elements.end_points.r = [0.0]
    elements.end_points.phi = [0.0]
    elements.end_points.z = [4.0]

    elements.centres.r = [0.0]
    elements.centres.phi = [0.0]
    elements.centres.z = [-1.0]

    resolution = 10
    reader.resolution = resolution
    points = reader._create_circular_geometry(elements, 0, is_full_circle=False)

    assert np.allclose(points[0], [0.0, 0.0, -6.0])
    assert np.allclose(points[-1], [0.0, 0.0, 4.0])
    assert np.all(points[:, 1] == 0.0)
    center = np.array([0.0, 0.0, -1.0])
    for p in points:  # All points should lie on a circle
        assert np.linalg.norm(np.array(p) - center) == pytest.approx(5.0)


def test_arc_of_circle_diagonal(elements, reader):
    elements.start_points.r = [5.0]
    elements.start_points.phi = [0.0]
    elements.start_points.z = [5.0]

    elements.intermediate_points.r = [np.sqrt(2) * 5.0]
    elements.intermediate_points.phi = [np.pi / 2]
    elements.intermediate_points.z = [0.0]

    elements.end_points.r = [5.0]
    elements.end_points.phi = [np.pi]
    elements.end_points.z = [-5.0]

    elements.centres.r = [0.0]
    elements.centres.phi = [0.0]
    elements.centres.z = [0.0]

    resolution = 10
    reader.resolution = resolution
    points = reader._create_circular_geometry(elements, 0, is_full_circle=False)

    assert np.allclose(points[0], [5.0, 0.0, 5.0])
    assert np.allclose(points[-1], [-5.0, 0.0, -5.0])
    center = np.array([0.0, 0.0, 0.0])
    for p in points:  # All points should lie on a circle
        assert np.linalg.norm(np.array(p) - center) == pytest.approx(np.sqrt(2) * 5.0)


def test_full_circle(elements, reader):
    elements.start_points.r = [3.0]
    elements.start_points.phi = [np.pi / 2]
    elements.start_points.z = [1.0]

    elements.intermediate_points.r = [3.0]
    elements.intermediate_points.phi = [np.pi / 4]
    elements.intermediate_points.z = [1.0]

    elements.centres.r = [0.0]
    elements.centres.phi = [0.0]
    elements.centres.z = [1.0]

    resolution = 10
    reader.resolution = resolution
    points = reader._create_circular_geometry(elements, 0, is_full_circle=True)

    assert np.allclose(points[0], [0.0, 3.0, 1.0])
    assert np.allclose(points[-1], [0.0, 3.0, 1.0])
    center = np.array([0.0, 0.0, 1.0])
    for p in points:  # All points should lie on a circle
        assert np.linalg.norm(np.array(p) - center) == pytest.approx(3.0)


def test_arc_valid(reader, base_points):
    p_start, p_intermediate, p_end, p_centre = base_points
    assert reader._are_circular_points_valid(
        p_start, p_intermediate, p_end, p_centre, False
    )


def test_arc_zero_radius(reader, base_points):
    _, p_intermediate, p_end, p_centre = base_points
    p_start = np.array([0.0, 0.0, 0.0])  # coincides with centre
    assert not reader._are_circular_points_valid(
        p_start, p_intermediate, p_end, p_centre, False
    )


def test_arc_intermediate_radius_not_equal(reader, base_points):
    p_start, _, p_end, p_centre = base_points
    p_intermediate = np.array([0.0, 6.0, 0.0])  # different radius
    assert not reader._are_circular_points_valid(
        p_start, p_intermediate, p_end, p_centre, False
    )


def test_arc_end_radius_not_equal(reader, base_points):
    p_start, p_intermediate, _, p_centre = base_points
    p_end = np.array([-6.0, 0.0, 0.0])  # different radius
    assert not reader._are_circular_points_valid(
        p_start, p_intermediate, p_end, p_centre, False
    )


def test_arc_degenerate_plane(reader, base_points):
    p_start, p_intermediate, p_end, p_centre = base_points
    p_intermediate = np.array([10.0, 0.0, 0.0])  # start and intermediate collinear
    assert not reader._are_circular_points_valid(
        p_start, p_intermediate, p_end, p_centre, False
    )


def test_arc_non_coplanar_end_point(reader, base_points):
    p_start, p_intermediate, p_end, p_centre = base_points
    p_end = np.array([-5.0, 0.0, 1.0])  # not coplanar
    assert not reader._are_circular_points_valid(
        p_start, p_intermediate, p_end, p_centre, False
    )


def test_full_circle_valid(reader, base_points):
    p_start, p_intermediate, p_end, p_centre = base_points
    p_end = None  # ignored for full circle
    assert reader._are_circular_points_valid(
        p_start, p_intermediate, p_end, p_centre, True
    )


def test_full_circle_radius_not_equal(reader, base_points):
    p_start, p_intermediate, p_end, p_centre = base_points
    p_end = None  # ignored for full circle
    p_intermediate = np.array([0.0, 6.0, 0.0])  # different radius
    assert not reader._are_circular_points_valid(
        p_start, p_intermediate, p_end, p_centre, True
    )
