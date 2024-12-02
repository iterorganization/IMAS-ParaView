from imaspy import DBEntry
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from vtkggdtools.plugins.wall_limiter import IMASPyWallLimiterReader


def test_load_limiters():
    """Test if limiters are loaded in the VTK Multiblock Dataset."""
    reader = IMASPyWallLimiterReader()
    entry = DBEntry(
        "imas:hdf5?path=/work/imas/shared/imasdb/ITER_MACHINE_DESCRIPTION/3/116000/5/",
        "r",
    )
    ids = entry.get("wall", lazy=True, autoconvert=False)
    reader._ids = ids
    reader.setup_ids()
    description_name = ids.description_2d[0].type.name
    limiter_name = ids.description_2d[0].limiter.type.name
    unit_name1 = ids.description_2d[0].limiter.unit[0].name
    unit_name2 = ids.description_2d[0].limiter.unit[0].name

    output = vtkMultiBlockDataSet()
    assert output.GetNumberOfBlocks() == 0

    output = vtkMultiBlockDataSet()
    reader._selected = [f"{description_name} / {limiter_name} / {unit_name1}"]
    reader._load_limiters(output)
    assert output.GetNumberOfBlocks() == 1

    output = vtkMultiBlockDataSet()
    reader._selected = [
        f"{description_name} / {limiter_name} / {unit_name1}",
        f"{description_name} / {limiter_name} / {unit_name2}",
    ]
    reader._load_limiters(output)
    assert output.GetNumberOfBlocks() == 2
