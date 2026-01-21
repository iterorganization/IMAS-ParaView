import imas
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkTable

from imas_paraview.plugins.time_dependent_0d import TimeDependent0DReader


def test_time_array():
    ids = imas.IDSFactory(version="4.0.0").new("wall")
    reader = TimeDependent0DReader()

    n_time_points = 5
    ids.time = 1.1 * np.arange(n_time_points)
    ids.global_quantities.power_incident = 2.2 * np.arange(n_time_points)
    ids.global_quantities.neutral.resize(2)
    ids.global_quantities.neutral[0].name = "neutral 1"
    ids.global_quantities.neutral[0].gas_puff = 3.3 * np.arange(n_time_points)
    ids.global_quantities.neutral[1].name = "neutral 2"
    ids.global_quantities.neutral[1].gas_puff = 4.4 * np.arange(n_time_points)
    ids.global_quantities.temperature = 5.5 * np.arange(n_time_points)

    pi_name = "Global_quantities Power_incident [W]"
    gp1_name = "Global_quantities Neutral (Neutral 1) Gas_puff [s^-1]"
    gp2_name = "Global_quantities Neutral (Neutral 2) Gas_puff [s^-1]"
    reader._ids = ids
    reader.setup_ids()

    reader._selected = [gp1_name, gp2_name, pi_name]
    for i in range(1, n_time_points):
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


def test_time_slice():
    ids = imas.IDSFactory(version="4.0.0").new("equilibrium")
    reader = TimeDependent0DReader()

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
    for i in range(1, n_time_points):
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
