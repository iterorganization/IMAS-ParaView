from imaspy import DBEntry
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from vtkggdtools.plugins.wall_limiter import IMASPyWallLimiterReader


def test_load_limiters():
    """Test if limiters are loaded in the VTK Multiblock Dataset."""
    reader = IMASPyWallLimiterReader()
    entry = DBEntry(
            "imas:hdf5?path=/work/imas/shared/imasdb/ITER_MACHINE_DESCRIPTION/3/116000/5/", "r"
    )
    ids = entry.get("wall", lazy=True, autoconvert=False)
    reader._ids = ids 
    reader.setup_ids()
    
    output = vtkMultiBlockDataSet()
    assert output.GetNumberOfBlocks() == 0

    output = vtkMultiBlockDataSet()
    reader._selected = ["Wall description at operating temperature (100 C) / Divertor and FW / First Wall"] 
    reader._load_limiters(output)
    assert output.GetNumberOfBlocks() == 1
    
    output = vtkMultiBlockDataSet()
    reader._selected = ["Wall description at operating temperature (100 C) / Divertor and FW / First Wall", 
    "Wall description at operating temperature (100 C) / Divertor and FW / Divertor"] 
    reader._load_limiters(output)
    assert output.GetNumberOfBlocks() == 2
