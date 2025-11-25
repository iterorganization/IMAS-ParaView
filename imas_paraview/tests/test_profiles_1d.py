import numpy as np
import pytest
from imas import DBEntry
from vtk.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkTable

from imas_paraview.plugins.profiles_1d import Profiles1DReader


@pytest.mark.external_data
def test_load_profiles_1d(test_data_dir):
    """Test if 1D profile structures are loaded in the VTK Table."""
    reader = Profiles1DReader()

    with DBEntry(test_data_dir / "iter_scenario_53298_seq1_DD4.nc", "r") as entry:
        ids = entry.get("core_profiles", autoconvert=False)
        reader._ids = ids
        reader.setup_ids()

        profile1 = ids.profiles_1d[0].electrons.temperature
        name1 = "Electrons Temperature"
        profile2 = ids.profiles_1d[0].j_total
        name2 = "J_total"

        # 1 selection
        output = vtkTable()
        reader._selected = [name1]
        reader._load_profiles(output)

        check_table(output, [(name1, profile1)])

        # 2 selections
        output = vtkTable()
        reader._selected = [name1, name2]
        reader._load_profiles(output)

        check_table(output, [(name1, profile1), (name2, profile2)])


def check_table(vtk_table, expected_profiles):
    """Assert that the data in the vtkTable matches the expected 1D profiles."""
    assert vtk_table.GetNumberOfColumns() == len(expected_profiles) + 1

    for name, profile in expected_profiles:
        array = vtk_to_numpy(vtk_table.GetColumnByName(name))
        assert np.array_equal(profile.ravel(), array)
