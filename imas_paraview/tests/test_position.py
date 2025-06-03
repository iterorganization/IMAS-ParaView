import imas
import numpy as np
import pytest
from imas import DBEntry
from vtkmodules.vtkCommonDataModel import vtkPolyData

from imas_paraview.plugins.position import PositionReader
from imas_paraview.util import pol_to_cart


@pytest.mark.skip(reason="no IMAS-Core available")
def test_load_position_magnetics():
    """Test if positions of magnetics barometry IDS are saved into vtkPolyData."""
    reader = PositionReader()
    entry = DBEntry(
        (
            "imas:hdf5?path=/work/imas/shared/imasdb/ITER_MACHINE_DESCRIPTION/"
            "/3/150100/5/"
        ),
        "r",
    )
    ids = entry.get("magnetics", lazy=True, autoconvert=False)

    reader._ids = ids
    reader.setup_ids()

    r1 = ids.b_field_pol_probe[0].position.r
    phi1 = ids.b_field_pol_probe[0].position.phi
    z1 = ids.b_field_pol_probe[0].position.z
    name1 = f"{ids.b_field_pol_probe[0].name} / {ids.b_field_pol_probe[0].identifier}"
    point1 = (*pol_to_cart(r1, phi1), z1)
    output = vtkPolyData()
    reader._selected = [name1]
    reader._load_position(output)
    assert output.GetNumberOfPoints() == 1
    assert np.all(np.isclose(point1, output.GetPoint(0)))


@pytest.mark.skip(reason="no IMAS-Core available")
def test_load_position_barometry():
    """Test if positions of gauges of barometry IDS are saved into vtkPolyData."""
    ids = imas.IDSFactory().barometry()
    ids.gauge.resize(2)

    r1 = 1.1
    phi1 = 2.2
    z1 = 3.3
    name1 = "gauge1"
    ids.gauge[0].name = name1
    ids.gauge[0].position.r = r1
    ids.gauge[0].position.phi = phi1
    ids.gauge[0].position.z = z1
    point1 = (*pol_to_cart(r1, phi1), z1)

    r2 = 4.4
    phi2 = 5.5
    z2 = 6.6
    name2 = "gauge2"
    ids.gauge[1].name = name2
    ids.gauge[1].position.r = r2
    ids.gauge[1].position.phi = phi2
    ids.gauge[1].position.z = z2
    point2 = (*pol_to_cart(r2, phi2), z2)

    reader = PositionReader()
    reader._ids = ids
    reader.setup_ids()

    output = vtkPolyData()
    reader._selected = [name1]
    reader._load_position(output)
    assert output.GetNumberOfPoints() == 1
    assert np.all(np.isclose(point1, output.GetPoint(0)))

    output = vtkPolyData()
    reader._selected = [name1, name2]
    reader._load_position(output)
    assert output.GetNumberOfPoints() == 2
    assert np.all(np.isclose(point1, output.GetPoint(0)))
    assert np.all(np.isclose(point2, output.GetPoint(1)))


def test_load_position_langmuir_probes():
    """Test if positions of embedded langmuir probes are saved into vtkPolyData."""
    ids = imas.IDSFactory().langmuir_probes()
    ids.embedded.resize(2)

    r1 = 1.1
    phi1 = 2.2
    z1 = 3.3
    name1 = "probe1"
    ids.embedded[0].name = name1
    ids.embedded[0].position.r = r1
    ids.embedded[0].position.phi = phi1
    ids.embedded[0].position.z = z1
    point1 = (*pol_to_cart(r1, phi1), z1)

    r2 = 4.4
    phi2 = 5.5
    z2 = 6.6
    name2 = "probe2"
    ids.embedded[1].name = name2
    ids.embedded[1].position.r = r2
    ids.embedded[1].position.phi = phi2
    ids.embedded[1].position.z = z2
    point2 = (*pol_to_cart(r2, phi2), z2)

    reader = PositionReader()
    reader._ids = ids
    reader.setup_ids()

    output = vtkPolyData()
    reader._selected = [name1]
    reader._load_position(output)
    assert output.GetNumberOfPoints() == 1
    assert np.all(np.isclose(point1, output.GetPoint(0)))

    output = vtkPolyData()
    reader._selected = [name1, name2]
    reader._load_position(output)
    assert output.GetNumberOfPoints() == 2
    assert np.all(np.isclose(point1, output.GetPoint(0)))
    assert np.all(np.isclose(point2, output.GetPoint(1)))
