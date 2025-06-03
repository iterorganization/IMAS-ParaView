from pathlib import Path

import imas
from asv_runner.benchmarks.mark import skip_benchmark

from imas_paraview.convert import Converter
from imas_paraview.tests.fill_ggd import fill_ids

SIZE_GRID = 10


def create_uri(backend):
    path = Path.cwd() / f"DB-{backend}"
    return f"imas:{backend.lower()}?path={path}"


class Convert:
    def setup(self):
        self.edge_profiles = imas.IDSFactory().edge_profiles()
        fill_ids(self.edge_profiles, grid_size=SIZE_GRID)

    def time_convert(self):
        converter = Converter(self.edge_profiles)
        converter.ggd_to_vtk(time_idx=0)


class ConvertNetCDF:
    def setup(self):
        es = imas.IDSFactory().edge_profiles()
        fill_ids(es, grid_size=SIZE_GRID)
        self.uri = Path.cwd() / "DB-netcdf.nc"
        with imas.DBEntry(self.uri, "w") as dbentry:
            dbentry.put(es)

    def time_load_and_convert(self):
        entry = imas.DBEntry(self.uri, "r")
        ids = entry.get("edge_profiles", lazy=False, autoconvert=False)
        converter = Converter(ids)
        converter.ggd_to_vtk(time_idx=0)

    def time_load_and_convert_lazy(self):
        entry = imas.DBEntry(self.uri, "r")
        ids = entry.get("edge_profiles", lazy=True, autoconvert=False)
        converter = Converter(ids)
        converter.ggd_to_vtk(time_idx=0)


# TODO: HDF5 and MDSPlus benchmarks are skipped as IMAS-Core is not yet available
class ConvertHDF5:
    def setup(self):
        es = imas.IDSFactory().edge_profiles()
        fill_ids(es, grid_size=SIZE_GRID)
        self.uri = create_uri("hdf5")
        with imas.DBEntry(self.uri, "w") as dbentry:
            dbentry.put(es)

    @skip_benchmark
    def time_load_and_convert(self):
        entry = imas.DBEntry(self.uri, "r")
        ids = entry.get("edge_profiles", lazy=False, autoconvert=False)
        converter = Converter(ids)
        converter.ggd_to_vtk(time_idx=0)

    @skip_benchmark
    def time_load_and_convert_lazy(self):
        entry = imas.DBEntry(self.uri, "r")
        ids = entry.get("edge_profiles", lazy=True, autoconvert=False)
        converter = Converter(ids)
        converter.ggd_to_vtk(time_idx=0)


class ConvertMDSPlus:
    def setup(self):
        es = imas.IDSFactory().edge_profiles()
        fill_ids(es, grid_size=SIZE_GRID)
        self.uri = create_uri("mdsplus")
        with imas.DBEntry(self.uri, "w") as dbentry:
            dbentry.put(es)

    @skip_benchmark
    def time_load_and_convert(self):
        entry = imas.DBEntry(self.uri, "r")
        ids = entry.get("edge_profiles", lazy=False, autoconvert=False)
        converter = Converter(ids)
        converter.ggd_to_vtk(time_idx=0)

    @skip_benchmark
    def time_load_and_convert_lazy(self):
        entry = imas.DBEntry(self.uri, "r")
        ids = entry.get("edge_profiles", lazy=True, autoconvert=False)
        converter = Converter(ids)
        converter.ggd_to_vtk(time_idx=0)
