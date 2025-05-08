import numpy as np
import pytest
from imas import DBEntry
from vtk.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from imas_paraview.plugins.wall_limiter import WallLimiterReader


@pytest.mark.skip(reason="no IMAS-Core available")
def test_load_limiters():
    """Test if limiters are loaded in the VTK Multiblock Dataset."""
    reader = WallLimiterReader()

    with DBEntry(
        "imas:hdf5?path=/work/imas/shared/imasdb/ITER_MACHINE_DESCRIPTION/3/116000/5/",
        "r",
    ) as entry:
        ids = entry.get("wall", lazy=True, autoconvert=False)
        reader._ids = ids
        reader.setup_ids()
        description_name = ids.description_2d[0].type.name
        limiter_name = ids.description_2d[0].limiter.type.name

        unit1 = ids.description_2d[0].limiter.unit[0]
        name1 = f"{description_name} / {limiter_name} / {unit1.name}"
        unit2 = ids.description_2d[0].limiter.unit[1]
        name2 = f"{description_name} / {limiter_name} / {unit2.name}"

        # 1 selection
        output = vtkMultiBlockDataSet()
        reader._selected = [name1]
        reader._load_limiters(output)
        assert output.GetNumberOfBlocks() == 1
        assert_values_match(unit1, output.GetBlock(0))

        # 2 selections
        output = vtkMultiBlockDataSet()
        reader._selected = [name1, name2]
        reader._load_limiters(output)
        assert output.GetNumberOfBlocks() == 2
        assert_values_match(unit1, output.GetBlock(0))
        assert_values_match(unit2, output.GetBlock(1))


def assert_values_match(unit, block):
    """Check if values in limiter unit match with the values in vtk block"""
    vtk_points = block.GetPoints()
    numpy_points = vtk_to_numpy(vtk_points.GetData())
    assert np.all(np.isclose(unit.outline.r, numpy_points[:, 0]))
    assert np.all(np.isclose(unit.outline.z, numpy_points[:, 2]))
