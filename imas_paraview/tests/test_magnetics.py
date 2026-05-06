import pytest
from imas import DBEntry
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from imas_paraview.plugins.magnetics import MagneticsReader


@pytest.mark.external_data
def test_load_magnetics(test_data_dir):
    """Test if all magnetics geometries are loaded correctly."""
    reader = MagneticsReader()

    FLUX_LOOP_NAME = "Flux loop (55.AD.00-MSA-1001)"
    POL_PROBE_NAME = "Poloidal field probe (55.A3.00-MLF-3001)"
    TOR_PROBE_NAME = "Toroidal field probe (55.AC Toroidal Coils)"
    ROGOWSKI_COIL_NAME = "Rogowski coil (55.AP.00-MRG-1201)"

    with DBEntry(test_data_dir / "iter_md_magnetics_150100_5.nc", "r") as entry:
        ids = entry.get("magnetics", autoconvert=False)
        reader._ids = ids
        reader.setup_ids()

        output = vtkMultiBlockDataSet()
        reader._selected = [
            FLUX_LOOP_NAME,
            POL_PROBE_NAME,
            TOR_PROBE_NAME,
            ROGOWSKI_COIL_NAME,
        ]
        reader._convert_to_vtk(output)

        assert output.GetNumberOfBlocks() == 4
        for i in range(4):
            block = output.GetBlock(i)
            assert block is not None
            assert block.GetNumberOfPoints() > 0
