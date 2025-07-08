from imas import DBEntry
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from imas_paraview.plugins.beam import BeamReader


def test_load_beam():
    """Test if limiters are loaded in the VTK Multiblock Dataset."""
    reader = BeamReader()
    with DBEntry(
        "/home/ITER/blokhus/public/imas_paraview_tests/iter_md-120000-1304.nc", "r"
    ) as entry:
        ids = entry.get("ec_launchers", autoconvert=False)
        reader._ids = ids
        reader.setup_ids()
        time_idx = 0
        name1 = str(ids.beam[0].name)
        name2 = str(ids.beam[1].name)

        # 1 selection
        output = vtkMultiBlockDataSet()
        reader._selected = [name1]
        reader._load_beam(output, time_idx)
        assert output.GetNumberOfBlocks() == 1

        # 2 selections
        output = vtkMultiBlockDataSet()
        reader._selected = [name1, name2]
        reader._load_beam(output, time_idx)
        assert output.GetNumberOfBlocks() == 2

        # All selected
        output = vtkMultiBlockDataSet()
        reader._selected = [str(beam.name) for beam in ids.beam]
        reader._load_beam(output, time_idx)
        assert output.GetNumberOfBlocks() == len(ids.beam)
