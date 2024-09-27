from pathlib import Path

import imaspy

from vtkggdtools.convert import ggd_to_vtk
from vtkggdtools.tests.fill_ggd import fill_ids

HDF5 = "HDF5"
MDSPLUS = "MDSplus"
MEMORY = "memory"
ASCII = "ASCII"

all_backends = [
    HDF5,
    MDSPLUS,
    MEMORY,
    ASCII,
]


def create_uri(backend):
    path = Path.cwd() / f"DB-{backend}"
    return f"imas:{backend.lower()}?path={path}"


class Convert:
    params = all_backends

    def setup(self, backend):
        self.edge_profiles = imaspy.IDSFactory().edge_profiles()
        fill_ids(self.edge_profiles, N=10)
        uri = create_uri(backend)
        with imaspy.DBEntry(uri, "w") as dbentry:
            dbentry.put(self.edge_profiles)

    def time_convert(self, backend):
        ggd_to_vtk(self.edge_profiles, 0)

    def time_load_and_convert(self, backend):
        uri = create_uri(backend)
        entry = imaspy.DBEntry(uri, "r")
        ids = entry.get("edge_profiles", lazy=False, autoconvert=False)
        ggd_to_vtk(ids, 0)

    def time_load_and_convert_lazy(self, backend):
        uri = create_uri(backend)
        entry = imaspy.DBEntry(uri, "r")
        ids = entry.get("edge_profiles", lazy=True, autoconvert=False)
        ggd_to_vtk(ids, 0)
