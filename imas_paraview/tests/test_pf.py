import imas
import numpy as np
import pytest
import vtk
from conftest import DD_VERSION
from vtkmodules.vtkCommonDataModel import vtkPolyData

from imas_paraview.plugins.pf import PFReader


@pytest.fixture(params=["pf_active", "pf_passive"])
def ids_geometry(request):
    """Returns a geometry IDS node from a pf_active or pf_passive IDS"""
    ids = imas.IDSFactory(version=DD_VERSION).new(request.param)
    if ids.metadata.name == "pf_active":
        node = ids.coil
    else:
        node = ids.loop
    node.resize(1)
    node[0].element.resize(1)
    return node[0].element[0].geometry


def get_points(poly: vtkPolyData):
    pts = poly.GetPoints()
    return [pts.GetPoint(i) for i in range(pts.GetNumberOfPoints())]


def get_line_cells(poly: vtkPolyData):
    cells = poly.GetLines()
    ids = []
    cells.InitTraversal()

    idlist = vtk.vtkIdList()
    while cells.GetNextCell(idlist):
        ids.append([idlist.GetId(i) for i in range(idlist.GetNumberOfIds())])
    return ids


def test_rectangle_geometry(ids_geometry):
    rectangle = ids_geometry.rectangle
    rectangle.r = 1.0
    rectangle.z = 2.0
    rectangle.width = 3.0
    rectangle.height = 4.0
    result = PFReader()._create_rectangle(rectangle)

    pts = get_points(result)
    assert len(pts) == 5
    assert pts[0] == (-0.5, 0.0, 0.0)
    assert pts[1] == (2.5, 0.0, 0.0)
    assert pts[2] == (2.5, 0.0, 4.0)
    assert pts[3] == (-0.5, 0.0, 4.0)
    assert pts[4] == (-0.5, 0.0, 0.0)

    cells = get_line_cells(result)
    assert len(cells) == 1
    assert len(cells[0]) == 5


def test_oblique_geometry(ids_geometry):
    oblique = ids_geometry.oblique
    oblique.r = 1.0
    oblique.z = 2.0
    oblique.length_alpha = 3.0
    oblique.length_beta = 4.0
    oblique.alpha = np.pi / 6.0
    oblique.beta = -np.pi / 6.0
    result = PFReader()._create_oblique(oblique)

    pts = get_points(result)
    assert len(pts) == 5
    expected_pts = np.array(
        [
            [1.0, 0.0, 2.0],
            [1.0 + 1.5 * np.sqrt(3), 0.0, 3.5],
            [3.0 + 1.5 * np.sqrt(3), 0.0, 3.5 + 2.0 * np.sqrt(3)],
            [3.0, 0.0, 2.0 + 2.0 * np.sqrt(3)],
            [1.0, 0.0, 2.0],
        ]
    )
    assert np.allclose(np.array(pts), expected_pts)

    cells = get_line_cells(result)
    assert len(cells) == 1
    assert len(cells[0]) == 5


def test_annulus_geometry(ids_geometry):
    annulus = ids_geometry.annulus
    annulus.r = 1.0
    annulus.z = 2.0
    annulus.radius_inner = 3.0
    annulus.radius_outer = 4.0
    resolution = 10
    result = PFReader()._create_annulus(annulus, resolution=resolution)

    pts = get_points(result)
    assert len(pts) == resolution * 2

    # Check if radius is correct
    center = np.array([1.0, 0.0, 2.0])
    outer_pts = pts[0::2]
    inner_pts = pts[1::2]
    for p in outer_pts:
        assert np.linalg.norm(np.array(p) - center) == pytest.approx(4.0)
    for p in inner_pts:
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
    result = PFReader()._create_thick_line(thick_line)

    pts = get_points(result)
    assert len(pts) == 5
    expected_pts = np.array(
        [
            [0.0, 0.0, 3.0],
            [2.0, 0.0, 5.0],
            [4.0, 0.0, 3.0],
            [2.0, 0.0, 1.0],
            [0.0, 0.0, 3.0],
        ]
    )
    assert np.allclose(np.array(pts), expected_pts)

    cells = get_line_cells(result)
    assert len(cells) == 1
    assert len(cells[0]) == 5
