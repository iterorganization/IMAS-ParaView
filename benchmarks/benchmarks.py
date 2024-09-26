from pathlib import Path

import imaspy

from vtkggdtools.convert import ggd_to_vtk
from vtkggdtools.tests.fill_ggd import fill_ids

from .utils import all_backends, create_dbentry, create_uri


class Get:
    params = all_backends
    param_names = ["backend"]

    def convert(self, backend):
        path = Path.cwd() / f"DB-{backend}"
        self.dbentry = create_dbentry(backend, path)
        edge_profiles = imaspy.IDSFactory().edge_profiles()
        fill_ids(edge_profiles, N=20)
        self.dbentry.put(edge_profiles)
        uri = create_uri(backend, path)
        ggd_to_vtk(uri, 0)
