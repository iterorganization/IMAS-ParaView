from pathlib import Path

import imaspy

from vtkggdtools.converter import ggd_to_vtk
from vtkggdtools.tests.fill_ggd import fill_ids

SIZE_GRID = 10


def create_uri(backend):
    path = Path.cwd() / f"DB-{backend}"
    return f"imas:{backend.lower()}?path={path}"


class Convert:
    def setup(self):
        self.edge_profiles = imaspy.IDSFactory().edge_profiles()
        fill_ids(self.edge_profiles, grid_size=SIZE_GRID)

    def time_convert(self):
        ggd_to_vtk(self.edge_profiles, time_idx=0)


class ConvertHDF5:
    def setup(self):
        es = imaspy.IDSFactory().edge_profiles()
        fill_ids(es, grid_size=SIZE_GRID)
        self.uri = create_uri("hdf5")
        with imaspy.DBEntry(self.uri, "w") as dbentry:
            dbentry.put(es)

    def time_load_and_convert(self):
        entry = imaspy.DBEntry(self.uri, "r")
        ids = entry.get("edge_profiles", lazy=False, autoconvert=False)
        ggd_to_vtk(ids, time_idx=0)

    def time_load_and_convert_lazy(self):
        entry = imaspy.DBEntry(self.uri, "r")
        ids = entry.get("edge_profiles", lazy=True, autoconvert=False)
        ggd_to_vtk(ids, time_idx=0)


class ConvertMDSPlus:
    def setup(self):
        es = imaspy.IDSFactory().edge_profiles()
        fill_ids(es, grid_size=SIZE_GRID)
        self.uri = create_uri("mdsplus")
        with imaspy.DBEntry(self.uri, "w") as dbentry:
            dbentry.put(es)

    def time_load_and_convert(self):
        entry = imaspy.DBEntry(self.uri, "r")
        ids = entry.get("edge_profiles", lazy=False, autoconvert=False)
        ggd_to_vtk(ids, time_idx=0)

    def time_load_and_convert_lazy(self):
        entry = imaspy.DBEntry(self.uri, "r")
        ids = entry.get("edge_profiles", lazy=True, autoconvert=False)
        ggd_to_vtk(ids, time_idx=0)
