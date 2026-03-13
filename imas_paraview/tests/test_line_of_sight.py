import numpy as np
import pytest
from imas import DBEntry
from vtk.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from imas_paraview.plugins.line_of_sight import LineOfSightReader
from imas_paraview.util import pol_to_cart


@pytest.mark.external_data
def test_load_los(test_data_dir):
    """Test if line_of_sight structures are loaded in the VTK Multiblock Dataset."""
    reader = LineOfSightReader()

    with DBEntry(test_data_dir / "iter_md_bolometer_150401_3.nc", "r") as entry:
        ids = entry.get("bolometer", autoconvert=False)
        reader._ids = ids
        reader.setup_ids()

        los1 = ids.channel[0].line_of_sight
        name1 = str(ids.channel[0].name)
        los2 = ids.channel[1].line_of_sight
        name2 = str(ids.channel[1].name)

        # 1 selection
        output = vtkMultiBlockDataSet()
        reader._selected = [name1]
        reader._load_los(output)
        assert output.GetNumberOfBlocks() == 1
        assert_values_match(los1, output.GetBlock(0))

        # 2 selections
        output = vtkMultiBlockDataSet()
        reader._selected = [name1, name2]
        reader._load_los(output)
        assert output.GetNumberOfBlocks() == 2
        assert_values_match(los1, output.GetBlock(0))
        assert_values_match(los2, output.GetBlock(1))

        # All selected
        output = vtkMultiBlockDataSet()
        reader._selected = [str(channel.name) for channel in ids.channel]
        reader._load_los(output)
        assert output.GetNumberOfBlocks() == len(ids.channel)
        for i, channel in enumerate(ids.channel):
            los = channel.line_of_sight
            assert_values_match(los, output.GetBlock(i))


@pytest.mark.external_data
def test_load_los_ece(test_data_dir):
    """Test if line_of_sight structures in ece are loaded in the VTK
    Multiblock Dataset."""
    reader = LineOfSightReader()
    with DBEntry(test_data_dir / "iter_md_ece_150601_22.nc", "r") as entry:
        ids = entry.get("ece", autoconvert=False)
        reader._ids = ids
        reader.setup_ids()

        # this ece IDS does not have line_of_sight structures for each channel,
        # instead the global line_of_sight is used
        los = ids.line_of_sight
        name1 = str(ids.channel[0].name)
        name2 = str(ids.channel[1].name)

        # 1 selection
        output = vtkMultiBlockDataSet()
        reader._selected = [name1]
        reader._load_los(output)
        assert output.GetNumberOfBlocks() == 1
        assert_values_match(los, output.GetBlock(0))

        # 2 selections
        output = vtkMultiBlockDataSet()
        reader._selected = [name1, name2]
        reader._load_los(output)
        assert output.GetNumberOfBlocks() == 2
        assert_values_match(los, output.GetBlock(0))
        assert_values_match(los, output.GetBlock(1))

        # All selected
        output = vtkMultiBlockDataSet()
        reader._selected = [str(channel.name) for channel in ids.channel]
        reader._load_los(output)
        assert output.GetNumberOfBlocks() == len(ids.channel)
        for i in range(len(ids.channel)):
            assert_values_match(los, output.GetBlock(i))


def assert_values_match(los, block):
    """Check if points in the line_of_sight match with the values in vtk block"""
    vtk_points = block.GetPoints()
    numpy_points = vtk_to_numpy(vtk_points.GetData())

    first_point = [
        *pol_to_cart(los.first_point.r, los.first_point.phi),
        los.first_point.z,
    ]
    second_point = [
        *pol_to_cart(los.second_point.r, los.second_point.phi),
        los.second_point.z,
    ]
    los_points = np.array([first_point, second_point])

    assert np.all(np.isclose(los_points, numpy_points))
