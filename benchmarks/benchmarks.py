from pathlib import Path

import imaspy

from vtkggdtools.convert import ggd_to_vtk
from vtkggdtools.tests.fill_ggd import fill_ids


def create_uri(backend):
    path = Path.cwd() / f"DB-{backend}"
    return f"imas:{backend.lower()}?path={path}"


class Convert:
    def setup(self):
        self.edge_profiles = imaspy.IDSFactory().edge_profiles()
        fill_ids(self.edge_profiles, N=3)

    def time_convert(self):
        ggd_to_vtk(self.edge_profiles, 0)


class ConvertHDF5:
    def setup(self):
        es = imaspy.IDSFactory().edge_profiles()
        fill_ids(es, N=3)
        self.uri = create_uri("hdf5")
        with imaspy.DBEntry(self.uri, "w") as dbentry:
            dbentry.put(es)

    def time_load_and_convert(self):
        entry = imaspy.DBEntry(self.uri, "r")
        ids = entry.get("edge_profiles", lazy=False, autoconvert=False)
        ggd_to_vtk(ids, 0)

    def time_load_and_convert_lazy(self):
        entry = imaspy.DBEntry(self.uri, "r")
        ids = entry.get("edge_profiles", lazy=True, autoconvert=False)
        ggd_to_vtk(ids, 0)


class ConvertMDSPlus:
    def setup(self):
        es = imaspy.IDSFactory().edge_profiles()
        fill_ids(es, N=3)
        self.uri = create_uri("mdsplus")
        with imaspy.DBEntry(self.uri, "w") as dbentry:
            dbentry.put(es)

    def time_load_and_convert(self):
        entry = imaspy.DBEntry(self.uri, "r")
        ids = entry.get("edge_profiles", lazy=False, autoconvert=False)
        ggd_to_vtk(ids, 0)

    def time_load_and_convert_lazy(self):
        entry = imaspy.DBEntry(self.uri, "r")
        ids = entry.get("edge_profiles", lazy=True, autoconvert=False)
        ggd_to_vtk(ids, 0)
