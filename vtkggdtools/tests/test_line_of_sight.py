import numpy as np
from imaspy import DBEntry
from vtk.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from vtkggdtools.plugins.line_of_sight import IMASPyLineOfSightReader
from vtkggdtools.util import pol_to_cart


def test_load_los():
    """Test if limiters are loaded in the VTK Multiblock Dataset."""
    reader = IMASPyLineOfSightReader()
    entry = DBEntry(
        "imas:hdf5?path=/work/imas/shared/imasdb/ITER_MACHINE_DESCRIPTION/3/150401/3/",
        "r",
    )
    ids = entry.get("bolometer", lazy=True, autoconvert=False)
    reader._ids = ids
    reader.setup_ids()

    los1 = ids.channel[0].line_of_sight
    name1 = ids.channel[0].name
    los2 = ids.channel[1].line_of_sight
    name2 = ids.channel[1].name

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
    reader._selected = []
    for channel in ids.channel:
        reader._selected.append(channel.name)
    reader._load_los(output)
    assert output.GetNumberOfBlocks() == len(ids.channel)
    for i, channel in enumerate(ids.channel):
        los = channel.line_of_sight
        assert_values_match(los, output.GetBlock(i))


def assert_values_match(los, block):
    """Check if values in limiter unit match with the values in vtk block"""
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
