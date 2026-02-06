import imas

from imas_paraview.convert import Converter
from imas_paraview.tests.fill_ggd import fill_ids
from imas_paraview.vtk2ggd import VTK2GGDConverter


def test_round_trip():
    ids_name = "edge_profiles"
    dd_version = "4.0.0"
    ids = imas.IDSFactory(version=dd_version).new(ids_name)
    fill_ids(ids, fill_ggd=False, create_3d_grid=True)
    converter = Converter(ids)
    vtk_pdsc = converter.ggd_to_vtk()  # get Partitioned Dataset Collection
    converter2 = VTK2GGDConverter([vtk_pdsc], ids_name, dd_version)
    ids2 = converter2.convert()
    # Check if grids are the same
    assert imas.util.calc_hash(ids.grid_ggd) == imas.util.calc_hash(ids2.grid_ggd)
