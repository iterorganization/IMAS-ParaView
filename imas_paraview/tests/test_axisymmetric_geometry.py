import imas
import numpy as np
import pytest
import vtk
from conftest import DD_VERSION
from vtkmodules.util.numpy_support import vtk_to_numpy

from imas_paraview.plugins.axisymmetric_geometry import AxisymmetricGeometryReader


@pytest.fixture(params=["pf_active", "pf_passive"])
def ids_geometry(request):
    """Returns a geometry IDS node from a pf_active or pf_passive IDS"""
    ids = imas.IDSFactory(version=DD_VERSION).new(request.param)
    node = ids.coil if ids.metadata.name == "pf_active" else ids.loop
    node.resize(1)
    node[0].element.resize(1)
    return node[0].element[0].geometry


def get_points(poly):
    points = poly.GetPoints()
    return vtk_to_numpy(points.GetData())


def get_line_cells(poly):
    cells = poly.GetLines()
    ids = []
    cells.InitTraversal()

    idlist = vtk.vtkIdList()
    while cells.GetNextCell(idlist):
        ids.append([idlist.GetId(i) for i in range(idlist.GetNumberOfIds())])
    return ids


def test_outline_geometry(ids_geometry):
    outline = ids_geometry.outline
    outline.r = [0.0, 1.0, 2.0, 3.0]
    outline.z = [4.0, 5.0, 6.0, 7.0]
    result = AxisymmetricGeometryReader()._create_outline(outline)

    points = get_points(result)
    assert len(points) == 4
    assert np.array_equal(points[0], [0.0, 0.0, 4.0])
    assert np.array_equal(points[1], [1.0, 0.0, 5.0])
    assert np.array_equal(points[2], [2.0, 0.0, 6.0])
    assert np.array_equal(points[3], [3.0, 0.0, 7.0])

    cells = get_line_cells(result)
    assert len(cells) == 4


def test_rectangle_geometry(ids_geometry):
    rectangle = ids_geometry.rectangle
    rectangle.r = 1.0
    rectangle.z = 2.0
    rectangle.width = 3.0
    rectangle.height = 4.0
    result = AxisymmetricGeometryReader()._create_rectangle(rectangle)

    points = get_points(result)
    assert len(points) == 4
    assert np.array_equal(points[0], [-0.5, 0.0, 0.0])
    assert np.array_equal(points[1], [2.5, 0.0, 0.0])
    assert np.array_equal(points[2], [2.5, 0.0, 4.0])
    assert np.array_equal(points[3], [-0.5, 0.0, 4.0])

    cells = get_line_cells(result)
    assert len(cells) == 4


def test_oblique_geometry(ids_geometry):
    oblique = ids_geometry.oblique
    oblique.r = 1.0
    oblique.z = 2.0
    oblique.length_alpha = 3.0
    oblique.length_beta = 4.0
    oblique.alpha = np.pi / 6.0
    oblique.beta = -np.pi / 6.0
    result = AxisymmetricGeometryReader()._create_oblique(oblique)

    points = get_points(result)
    assert len(points) == 4
    expected_pts = [
        [1.0, 0.0, 2.0],
        [1.0 + 1.5 * np.sqrt(3), 0.0, 3.5],
        [3.0 + 1.5 * np.sqrt(3), 0.0, 3.5 + 2.0 * np.sqrt(3)],
        [3.0, 0.0, 2.0 + 2.0 * np.sqrt(3)],
    ]
    assert np.allclose(np.array(points), expected_pts)

    cells = get_line_cells(result)
    assert len(cells) == 4


def test_arcs_of_circle_geometry(ids_geometry):
    arcs = ids_geometry.arcs_of_circle
    arcs.r = [0.0, 1.0]
    arcs.z = [0.0, 0.0]
    arcs.curvature_radii = [0.5, 0.5]
    resolution = 10

    reader = AxisymmetricGeometryReader()
    reader.resolution = resolution
    result = reader._create_arcs_of_circle(arcs)
    points = get_points(result)

    # All points should lie on a circle
    center = np.array([0.5, 0.0, 0.0])
    for p in points:
        assert np.linalg.norm(np.array(p) - center) == pytest.approx(0.5)
    assert len(points) == 2 * resolution

    cells = get_line_cells(result)
    assert len(cells) == 2 * resolution


def test_annulus_geometry(ids_geometry):
    annulus = ids_geometry.annulus
    annulus.r = 1.0
    annulus.z = 2.0
    annulus.radius_inner = 3.0
    annulus.radius_outer = 4.0
    resolution = 10
    reader = AxisymmetricGeometryReader()
    reader.resolution = resolution
    result = reader._create_annulus(annulus)

    points = get_points(result)
    assert len(points) == resolution * 2

    # Check if radius is correct
    center = np.array([1.0, 0.0, 2.0])
    outer_points = points[resolution:]
    inner_points = points[:resolution]
    for p in outer_points:
        assert np.linalg.norm(np.array(p) - center) == pytest.approx(4.0)
    for p in inner_points:
        assert np.linalg.norm(np.array(p) - center) == pytest.approx(3.0)

    cells = get_line_cells(result)
    assert len(cells) == 2
    assert len(cells[0]) == resolution + 1
    assert len(cells[1]) == resolution + 1


def test_thick_line_geometry(ids_geometry):
    thick_line = ids_geometry.thick_line
    thick_line.first_point.r = 1.0
    thick_line.first_point.z = 2.0
    thick_line.second_point.r = 3.0
    thick_line.second_point.z = 4.0
    thick_line.thickness = 2 * np.sqrt(2.0)
    result = AxisymmetricGeometryReader()._create_thick_line(thick_line)

    points = get_points(result)
    assert len(points) == 4
    expected_pts = [
        [0.0, 0.0, 3.0],
        [2.0, 0.0, 5.0],
        [4.0, 0.0, 3.0],
        [2.0, 0.0, 1.0],
    ]
    assert np.allclose(np.array(points), expected_pts)

    cells = get_line_cells(result)
    assert len(cells) == 4
