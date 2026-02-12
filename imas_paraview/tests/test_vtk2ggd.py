import imas

from imas_paraview.convert import Converter
from imas_paraview.tests.fill_ggd import fill_ids
from imas_paraview.vtk2ggd import VTK2GGDConverter


def test_round_trip():
    ids_name = "edge_profiles"
    dd_version = "4.0.0"
    ids = imas.IDSFactory(version=dd_version).new(ids_name)
    fill_ids(ids, fill_ggd=False, create_3d_grid=True)

    # Convert IDS to VTK
    converter = Converter(ids)
    vtk_pdsc = converter.ggd_to_vtk()

    # Convert VTK back to IDS
    converter2 = VTK2GGDConverter([vtk_pdsc], ids_name, dd_version=dd_version)
    ids2 = converter2.convert()

    # Check if the two IDSs are the same
    assert list(imas.util.idsdiffgen(ids, ids2)) == []
