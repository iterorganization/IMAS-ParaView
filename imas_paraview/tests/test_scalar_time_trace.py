import imas
import numpy as np
import pytest
from imas.ids_defs import IDS_TIME_MODE_HOMOGENEOUS
from vtk.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkTable

from imas_paraview.plugins.scalar_time_trace import ScalarTimeTraceReader


@pytest.fixture
def reader():
    return ScalarTimeTraceReader()


@pytest.fixture
def wall_ids():
    ids = imas.IDSFactory(version="4.1.0").new("wall")
    ids.ids_properties.homogeneous_time = IDS_TIME_MODE_HOMOGENEOUS
    n_time_points = 5
    ids.time = 1.1 * np.arange(n_time_points)
    ids.global_quantities.power_incident = 2.2 * np.arange(n_time_points)
    ids.global_quantities.neutral.resize(2)
    ids.global_quantities.neutral[0].name = "neutral 1"
    ids.global_quantities.neutral[0].gas_puff = 3.3 * np.arange(n_time_points)
    ids.global_quantities.neutral[1].name = "neutral 2"
    ids.global_quantities.neutral[1].gas_puff = 4.4 * np.arange(n_time_points)
    ids.global_quantities.temperature = 5.5 * np.arange(n_time_points)
    return ids


def test_time_array(wall_ids, reader):
    pi_name = "Global_quantities Power_incident [W]"
    gp1_name = "Global_quantities Neutral (Neutral 1) Gas_puff [s^-1]"
    gp2_name = "Global_quantities Neutral (Neutral 2) Gas_puff [s^-1]"
    reader._ids = wall_ids

    reader.setup_ids()

    reader._selected = [gp1_name, gp2_name, pi_name]
    for i in range(5):
        output = vtkTable()
        reader._load_time_dependent_data(output, 1.1 * i)

        assert output.GetNumberOfRows() == i + 1
        assert output.GetNumberOfColumns() == 4
        time_col = vtk_to_numpy(output.GetColumnByName("Time [s]"))
        pi_col = vtk_to_numpy(output.GetColumnByName(pi_name))
        gp1_col = vtk_to_numpy(output.GetColumnByName(gp1_name))
        gp2_col = vtk_to_numpy(output.GetColumnByName(gp2_name))

        assert np.all(time_col == np.array(1.1 * np.arange(i + 1)))
        assert np.all(pi_col == np.array(2.2 * np.arange(i + 1)))
        assert np.all(gp1_col == np.array(3.3 * np.arange(i + 1)))
        assert np.all(gp2_col == np.array(4.4 * np.arange(i + 1)))


def test_time_array_full_time_trace(wall_ids, reader):
    pi_name = "Global_quantities Power_incident [W]"
    gp1_name = "Global_quantities Neutral (Neutral 1) Gas_puff [s^-1]"
    gp2_name = "Global_quantities Neutral (Neutral 2) Gas_puff [s^-1]"
    reader._ids = wall_ids
    reader.setup_ids()
    reader._selected = [gp1_name, gp2_name, pi_name]
    reader._show_full_time_trace = True
    output = vtkTable()
    reader._load_time_dependent_data(output, 2.2)

    assert output.GetNumberOfRows() == 6  # 5 times + 1 duplicate
    assert output.GetNumberOfColumns() == 5  # time + marker + 3 selected columns
    time_col = vtk_to_numpy(output.GetColumnByName("Time [s]"))
    pi_col = vtk_to_numpy(output.GetColumnByName(pi_name))
    gp1_col = vtk_to_numpy(output.GetColumnByName(gp1_name))
    gp2_col = vtk_to_numpy(output.GetColumnByName(gp2_name))
    marker = vtk_to_numpy(output.GetColumnByName("Time Marker"))

    assert np.allclose(time_col, [0, 1.1, 2.2, 2.2, 3.3, 4.4])
    assert np.allclose(pi_col, [0.0, 2.2, 4.4, 4.4, 6.6, 8.8])
    assert np.allclose(gp1_col, [0.0, 3.3, 6.6, 6.6, 9.9, 13.2])
    assert np.allclose(gp2_col, [0.0, 4.4, 8.8, 8.8, 13.2, 17.6])
    assert np.allclose(
        marker,
        [np.nan, np.nan, 0.0, 17.6 * 1.01, np.nan, np.nan],
        equal_nan=True,
    )


def test_time_slice(reader):
    ids = imas.IDSFactory(version="4.1.0").new("equilibrium")
    ids.ids_properties.homogeneous_time = IDS_TIME_MODE_HOMOGENEOUS
    n_time_points = 5
    ids.time = 1.1 * np.arange(n_time_points)
    ids.time_slice.resize(n_time_points)

    for idx, time_slice in enumerate(ids.time_slice):
        time_slice.global_quantities.ip = 2.0 * idx
        time_slice.global_quantities.beta_tor_norm = 3.0 * idx
        time_slice.boundary.psi = 4.0 * idx

        time_slice.boundary.gap.resize(2)
        time_slice.boundary.gap[0].name = "gap 1"
        time_slice.boundary.gap[1].name = "gap 2"
        time_slice.boundary.gap[0].r = 5.0 * idx
        time_slice.boundary.gap[1].r = 6.0 * idx

    reader._ids = ids
    reader.setup_ids()

    psi_name = "Boundary Psi [Wb]"
    ip_name = "Global_quantities Ip [A]"
    gap1_name = "Boundary Gap (Gap 1) R [m]"
    gap2_name = "Boundary Gap (Gap 2) R [m]"
    reader._selected = [psi_name, ip_name, gap1_name, gap2_name]
    for i in range(n_time_points):
        output = vtkTable()
        reader._load_time_dependent_data(output, 1.1 * i)

        assert output.GetNumberOfRows() == i + 1
        assert output.GetNumberOfColumns() == 5
        time_col = vtk_to_numpy(output.GetColumnByName("Time [s]"))
        psi_col = vtk_to_numpy(output.GetColumnByName(psi_name))
        ip_col = vtk_to_numpy(output.GetColumnByName(ip_name))
        gap1_col = vtk_to_numpy(output.GetColumnByName(gap1_name))
        gap2_col = vtk_to_numpy(output.GetColumnByName(gap2_name))

        assert np.all(time_col == np.array(1.1 * np.arange(i + 1)))
        assert np.all(ip_col == np.array(2.0 * np.arange(i + 1)))
        assert np.all(psi_col == np.array(4.0 * np.arange(i + 1)))
        assert np.all(gap1_col == np.array(5.0 * np.arange(i + 1)))
        assert np.all(gap2_col == np.array(6.0 * np.arange(i + 1)))


def test_nested_aos_time_slice(reader):
    ids = imas.IDSFactory(version="4.1.0").new("camera_ir")
    ids.ids_properties.homogeneous_time = IDS_TIME_MODE_HOMOGENEOUS
    n_time_points = 5
    ids.time = 1.1 * np.arange(n_time_points)
    ids.channel.resize(2)
    ids.channel[0].name = "channel 1"
    ids.channel[1].name = "channel 2"
    for channel in ids.channel:
        channel.camera.resize(2)
        channel.camera[0].name = "camera 1"
        channel.camera[1].name = "camera 2"
        for camera in channel.camera:
            camera.frame.resize(n_time_points)
            for time_idx in range(n_time_points):
                camera.frame[time_idx].filter.resize(2)
                camera.frame[time_idx].filter[0].wavelength_central = 10
                camera.frame[time_idx].filter[1].wavelength_central = 20

    reader._ids = ids
    reader.setup_ids()

    f1_name = "Channel (Channel 1) Camera (Camera 1) Filter (#1) Wavelength_central [m]"
    f2_name = "Channel (Channel 2) Camera (Camera 2) Filter (#2) Wavelength_central [m]"
    reader._selected = [f1_name, f2_name]
    for i in range(n_time_points):
        output = vtkTable()
        reader._load_time_dependent_data(output, 1.1 * i)

        assert output.GetNumberOfRows() == i + 1
        assert output.GetNumberOfColumns() == 3
        time_col = vtk_to_numpy(output.GetColumnByName("Time [s]"))
        f1_col = vtk_to_numpy(output.GetColumnByName(f1_name))
        f2_col = vtk_to_numpy(output.GetColumnByName(f2_name))

        assert np.all(time_col == np.array(1.1 * np.arange(i + 1)))
        assert np.all(f1_col == 10)
        assert np.all(f2_col == 20)
