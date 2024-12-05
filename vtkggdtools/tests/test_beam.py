from imaspy import DBEntry
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from vtkggdtools.plugins.beam import IMASPyBeamReader


def test_load_limiters():
    """Test if limiters are loaded in the VTK Multiblock Dataset."""
    reader = IMASPyBeamReader()
    entry = DBEntry(
        (
            "imas:hdf5?path=/work/imas/shared/imasdb/ITER_MACHINE_DESCRIPTION/"
            "3/120000/1304/"
        ),
        "r",
    )
    ids = entry.get("wall", lazy=True, autoconvert=False)
    reader._ids = ids
    reader.setup_ids()
    name1 = ids.beam[0].name
    name2 = ids.beam[1].name

    # 1 selection
    output = vtkMultiBlockDataSet()
    reader._selected = [name1]
    reader._load_beam(output)
    assert output.GetNumberOfBlocks() == 1

    # 2 selections
    output = vtkMultiBlockDataSet()
    reader._selected = [name1, name2]
    reader._load_beam(output)
    assert output.GetNumberOfBlocks() == 2

    # All selected
    output = vtkMultiBlockDataSet()
    reader._selected = []
    for beam in ids.beam:
        reader._selected.append(beam.name)
    reader._load_beam(output)
    assert output.GetNumberOfBlocks() == len(ids.beam)
