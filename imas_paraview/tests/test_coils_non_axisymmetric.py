import imas
import numpy as np
import pytest
from conftest import DD_VERSION

from imas_paraview.plugins.coils_non_axisymmetric import CoilsNonAxisymmetricReader


@pytest.fixture
def conductor():
    """Returns a conductor IDS node of a coils_non_axisymmetric IDS."""
    ids = imas.IDSFactory(version=DD_VERSION).new("coils_non_axisymmetric")
    ids.coil.resize(1)
    ids.coil[0].conductor.resize(1)
    return ids.coil[0].conductor[0]


@pytest.fixture
def elements(conductor):
    """Returns a elements IDS node of a coils_non_axisymmetric conductor"""
    return conductor.elements


@pytest.fixture
def reader():
    return CoilsNonAxisymmetricReader()


def test_create_line_frame_with_intermediate_point(elements, reader):
    elements.start_points.r = [10.0]
    elements.start_points.phi = [0.0]
    elements.start_points.z = [0.0]
    elements.end_points.r = [10.0]
    elements.end_points.phi = [0.0]
    elements.end_points.z = [5.0]
    elements.intermediate_points.r = [11.0]
    elements.intermediate_points.phi = [0.0]
    elements.intermediate_points.z = [0.0]

    normal, binormal, tangent = reader._create_line_frame(elements, 0, None)

    assert np.allclose(tangent, [0.0, 0.0, 1.0])
    assert np.allclose(normal, [1.0, 0.0, 0.0])
    assert np.allclose(binormal, [0.0, 1.0, 0.0])


def test_create_line_frame_missing_intermediate_no_prev(elements, reader):
    """If intermediate_points is not filled and there is no previous frame, an
    arbitrary but valid perpendicular frame should be returned."""
    elements.start_points.r = [10.0]
    elements.start_points.phi = [0.0]
    elements.start_points.z = [0.0]
    elements.end_points.r = [10.0]
    elements.end_points.phi = [0.0]
    elements.end_points.z = [5.0]

    normal, binormal, tangent = reader._create_line_frame(elements, 0, None)

    assert np.allclose(tangent, [0.0, 0.0, 1.0])
    assert np.isclose(np.dot(normal, tangent), 0.0)
    assert np.isclose(np.dot(binormal, tangent), 0.0)
    assert np.isclose(np.dot(normal, binormal), 0.0)
    assert np.isclose(np.linalg.norm(normal), 1.0)
    assert np.isclose(np.linalg.norm(binormal), 1.0)


def test_create_line_frame_missing_intermediate_uses_prev_frame(elements, reader):
    """If intermediate_points is not filled but a previous frame is available, the
    frame should be parallel-transported from the previous element instead."""
    elements.start_points.r = [10.0]
    elements.start_points.phi = [0.0]
    elements.start_points.z = [0.0]
    elements.end_points.r = [10.0]
    elements.end_points.phi = [0.0]
    elements.end_points.z = [5.0]

    prev_frame = (np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), [0.0, 0.0, 1.0])
    normal, binormal, tangent = reader._create_line_frame(elements, 0, prev_frame)

    # Tangent is unchanged, so the frame should be transported without any rotation
    assert np.allclose(tangent, [0.0, 0.0, 1.0])
    assert np.allclose(normal, [1.0, 0.0, 0.0])
    assert np.allclose(binormal, [0.0, 1.0, 0.0])


def test_parallel_transport_frame_rotation(reader):
    """Rotating the tangent by 90 degrees should rotate the frame consistently."""
    prev_tangent = np.array([0.0, 0.0, 1.0])
    prev_normal = np.array([1.0, 0.0, 0.0])
    prev_binormal = np.array([0.0, 1.0, 0.0])
    new_tangent = np.array([1.0, 0.0, 0.0])

    normal, binormal = reader._parallel_transport_frame(
        prev_normal, prev_binormal, prev_tangent, new_tangent
    )

    assert np.isclose(np.dot(normal, new_tangent), 0.0)
    assert np.isclose(np.dot(binormal, new_tangent), 0.0)
    assert np.isclose(np.linalg.norm(normal), 1.0)
    assert np.isclose(np.linalg.norm(binormal), 1.0)


def test_create_circular_frames(elements, reader):
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

    reader.resolution = 4
    normals, binormals, tangents = reader._create_circular_frames(
        elements, 0, is_full_circle=False
    )

    # Normal should always point from the path point towards the centre
    points = reader._create_circular_geometry(elements, 0, is_full_circle=False)
    centre = np.array([0.0, 0.0, 5.0])
    for point, normal in zip(points, normals):
        assert np.allclose(centre - point, 10.0 * normal)
    # Binormal is constant (the z-axis, since the arc is in the horizontal plane)
    assert np.allclose(binormals, [0.0, 0.0, 1.0])
    # Tangents are perpendicular to the (constant-radius) normal at every point
    for normal, tangent in zip(normals, tangents):
        assert np.isclose(np.dot(normal, tangent), 0.0)


def test_polygon_cross_section_on_line_segment(conductor, reader):
    """A rectangular polygon outline swept along a straight line segment should
    produce a box with 8 vertices and 6 quad faces."""
    conductor.elements.types = [1]
    conductor.elements.start_points.r = [10.0]
    conductor.elements.start_points.phi = [0.0]
    conductor.elements.start_points.z = [0.0]
    conductor.elements.end_points.r = [10.0]
    conductor.elements.end_points.phi = [0.0]
    conductor.elements.end_points.z = [5.0]
    conductor.elements.intermediate_points.r = [11.0]
    conductor.elements.intermediate_points.phi = [0.0]
    conductor.elements.intermediate_points.z = [0.0]

    conductor.cross_section.resize(1)
    cross_section = conductor.cross_section[0]
    cross_section.geometry_type.index = 1  # polygon
    cross_section.outline.normal = [0.1, -0.1, -0.1, 0.1]
    cross_section.outline.binormal = [0.05, 0.05, -0.05, -0.05]

    output = reader.create_conductor_geometry(conductor)

    assert output.GetNumberOfPoints() == 8
    assert output.GetNumberOfCells() == 6
    x_min, x_max, y_min, y_max, z_min, z_max = output.GetBounds()
    assert np.isclose(x_min, 9.9)
    assert np.isclose(x_max, 10.1)
    assert np.isclose(y_min, -0.05)
    assert np.isclose(y_max, 0.05)
    assert np.isclose(z_min, 0.0)
    assert np.isclose(z_max, 5.0)


def test_polygon_cross_section_without_intermediate_point(conductor, reader):
    """The polygon cross-section should still be produced (via the fallback
    orientation) even when the line segment does not provide an intermediate
    point."""
    conductor.elements.types = [1]
    conductor.elements.start_points.r = [10.0]
    conductor.elements.start_points.phi = [0.0]
    conductor.elements.start_points.z = [0.0]
    conductor.elements.end_points.r = [10.0]
    conductor.elements.end_points.phi = [0.0]
    conductor.elements.end_points.z = [5.0]

    conductor.cross_section.resize(1)
    cross_section = conductor.cross_section[0]
    cross_section.geometry_type.index = 1  # polygon
    cross_section.outline.normal = [0.1, -0.1, -0.1, 0.1]
    cross_section.outline.binormal = [0.05, 0.05, -0.05, -0.05]

    output = reader.create_conductor_geometry(conductor)

    assert output.GetNumberOfPoints() == 8
    assert output.GetNumberOfCells() == 6


def test_polygon_cross_section_on_arc_preserves_radius(conductor, reader):
    """Sweeping a polygon cross-section along an arc should keep the cross-section
    at a fixed radial offset from the arc's centreline at every path point."""
    conductor.elements.types = [2]
    conductor.elements.start_points.r = [10.0]
    conductor.elements.start_points.phi = [0.0]
    conductor.elements.start_points.z = [5.0]
    conductor.elements.intermediate_points.r = [10.0]
    conductor.elements.intermediate_points.phi = [np.pi / 2]
    conductor.elements.intermediate_points.z = [5.0]
    conductor.elements.end_points.r = [10.0]
    conductor.elements.end_points.phi = [np.pi]
    conductor.elements.end_points.z = [5.0]
    conductor.elements.centres.r = [0.0]
    conductor.elements.centres.phi = [0.0]
    conductor.elements.centres.z = [5.0]

    conductor.cross_section.resize(1)
    cross_section = conductor.cross_section[0]
    cross_section.geometry_type.index = 1  # polygon
    cross_section.outline.normal = [0.1, -0.1, -0.1, 0.1]
    cross_section.outline.binormal = [0.05, 0.05, -0.05, -0.05]

    reader.resolution = 8
    output = reader.create_conductor_geometry(conductor)

    for i in range(output.GetNumberOfPoints()):
        x, y, z = output.GetPoint(i)
        radius = np.hypot(x, y)
        assert radius == pytest.approx(9.9, abs=1e-6) or radius == pytest.approx(
            10.1, abs=1e-6
        )
        assert z == pytest.approx(4.95, abs=1e-6) or z == pytest.approx(5.05, abs=1e-6)


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
