import numpy as np
from imaspy import DBEntry
from vtk.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkPartitionedDataSetCollection

from vtkggdtools.plugins.profiles_2d import IMASPyProfiles2DReader


def test_load_2d_profiles():
    """Test if 2D profiles structures are loaded in the VTK
    PartitionedDatasetCollection."""
    reader = IMASPyProfiles2DReader()

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
