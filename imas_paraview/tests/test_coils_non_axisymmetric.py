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


def test_create_line_segment(elements):
    elements.start_points.r = [10.0]
    elements.start_points.phi = [0.0]
    elements.start_points.z = [5.0]

    elements.end_points.r = [10.0]
    elements.end_points.phi = [np.pi / 2]
    elements.end_points.z = [3.0]

    points = CoilsNonAxisymmetricReader()._create_line_segment(elements, 0)

    assert len(points) == 2
    assert np.allclose(points[0], [10.0, 0.0, 5.0])
    assert np.allclose(points[1], [0.0, 10.0, 3.0])


def test_create_arc_segment(elements):
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
    reader = CoilsNonAxisymmetricReader()
    reader.resolution = resolution
    points = reader._create_circular_geometry(elements, 0, is_full_circle=False)

    assert len(points) == resolution
    assert np.allclose(points[0], [10.0, 0.0, 5.0])
    assert np.allclose(points[-1], [-10.0, 0.0, 5.0])
    center = np.array([0.0, 0.0, 5.0])
    for p in points:  # All points should lie on a circle
        assert np.linalg.norm(np.array(p) - center) == pytest.approx(10.0)


def test_create_full_circle(elements):
    elements.start_points.r = [3.0]
    elements.start_points.phi = [np.pi / 2]
    elements.start_points.z = [1.0]

    elements.intermediate_points.r = [5.0]
    elements.intermediate_points.phi = [np.pi / 4]
    elements.intermediate_points.z = [1.0]

    elements.centres.r = [0.0]
    elements.centres.phi = [0.0]
    elements.centres.z = [1.0]

    resolution = 10
    reader = CoilsNonAxisymmetricReader()
    reader.resolution = resolution
    points = reader._create_circular_geometry(elements, 0, is_full_circle=True)

    assert len(points) == resolution
    assert np.allclose(points[0], [0.0, 3.0, 1.0])
    assert np.allclose(points[-1], [0.0, 3.0, 1.0])
    center = np.array([0.0, 0.0, 1.0])
    for p in points:  # All points should lie on a circle
        assert np.linalg.norm(np.array(p) - center) == pytest.approx(3.0)
