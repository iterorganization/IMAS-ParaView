import imas
import numpy as np
import pytest
from imas import DBEntry
from imas.backends.imas_core.db_entry_helpers import IDS_TIME_MODE_HOMOGENEOUS
from vtk.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkPartitionedDataSetCollection

from imas_paraview.plugins.profiles_2d import Profiles2DReader


@pytest.mark.skip(reason="no IMAS-Core available")
def test_load_profiles():
    """Test if 2D profiles structures are loaded in the VTK
    PartitionedDatasetCollection."""
    reader = Profiles2DReader()

    with DBEntry(
        "imas:hdf5?path=/work/imas/shared/imasdb/ITER_SCENARIOS/3/110004/1",
        "r",
    ) as entry:
        ids = entry.get("equilibrium", lazy=True, autoconvert=False)
        reader._ids = ids
        reader.setup_ids()

        profile1 = ids.time_slice[0].profiles_2d[0].psi
        name1 = "Psi"
        profile2 = ids.time_slice[0].profiles_2d[0].j_tor
        name2 = "J_tor"

        # 1 selection
        output = vtkPartitionedDataSetCollection()
        reader._selected = [name1]
        reader._load_profiles(output)

        assert output.GetNumberOfPartitionedDataSets() == 1
        check_ugrid(output.GetPartitionedDataSet(0).GetPartition(0), profile1, reader)

        # 2 selections
        output = vtkPartitionedDataSetCollection()
        reader._selected = [name1, name2]
        reader._load_profiles(output)

        assert output.GetNumberOfPartitionedDataSets() == 2
        check_ugrid(output.GetPartitionedDataSet(0).GetPartition(0), profile1, reader)
        check_ugrid(output.GetPartitionedDataSet(1).GetPartition(0), profile2, reader)


def check_ugrid(ugrid, expected_profile, reader):
    """assert that the data in the ugrid matches the IDS profiles arrays."""
    num_points = ugrid.GetNumberOfPoints()
    vtk_scalars = ugrid.GetPointData().GetScalars()
    scalars_array = vtk_to_numpy(vtk_scalars)
    points_np = np.array([ugrid.GetPoint(i) for i in range(num_points)])
    r_vtk, z_vtk = points_np[:, 0], points_np[:, 2]

    assert np.array_equal(expected_profile.ravel(), scalars_array)
    assert np.array_equal(reader.r.ravel(), r_vtk)
    assert np.array_equal(reader.z.ravel(), z_vtk)


def test_load_profiles_dummy_equilibrium():
    ids = imas.IDSFactory(version="4.0.0").new("equilibrium")
    ids.time = [0]
    ids.ids_properties.homogeneous_time = IDS_TIME_MODE_HOMOGENEOUS
    ids.time_slice.resize(1)
    ids.time_slice[0].profiles_2d.resize(1)
    profiles_2d = ids.time_slice[0].profiles_2d[0]
    profiles_2d.r = np.array([[1.1, 2.2], [3.3, 4.4]])
    profiles_2d.z = np.array([[5.5, 6.6], [7.7, 8.8]])
    profiles_2d.psi = np.array([[10.1, 10.2], [10.3, 10.4]])

    reader = Profiles2DReader()
    reader._ids = ids
    reader.setup_ids()
    output = vtkPartitionedDataSetCollection()
    reader._selected = ["Psi"]
    reader._load_profiles(output)
    assert output.GetNumberOfPartitionedDataSets() == 1
    profile = ids.time_slice[0].profiles_2d[0].psi
    check_ugrid(output.GetPartitionedDataSet(0).GetPartition(0), profile, reader)


def test_load_profiles_dummy_core_profiles():
    ids = imas.IDSFactory(version="4.0.0").new("core_profiles")
    ids.time = [0]
    ids.ids_properties.homogeneous_time = IDS_TIME_MODE_HOMOGENEOUS
    ids.profiles_2d.resize(1)
    ids.profiles_2d[0].ion.resize(1)
    profiles_2d = ids.profiles_2d[0]
    profiles_2d.grid.dim1 = [1, 2]
    profiles_2d.grid.dim2 = [3, 4]
    profiles_2d.grid_type = 1
    profiles_2d.ion[0].pressure = np.array([[10.1, 10.2], [10.3, 10.4]])

    reader = Profiles2DReader()
    reader._ids = ids
    reader.setup_ids()
    output = vtkPartitionedDataSetCollection()
    reader._selected = ["Ion Pressure"]
    reader._load_profiles(output)
    assert output.GetNumberOfPartitionedDataSets() == 1
    profile = profiles_2d.ion[0].pressure
    check_ugrid(output.GetPartitionedDataSet(0).GetPartition(0), profile, reader)
