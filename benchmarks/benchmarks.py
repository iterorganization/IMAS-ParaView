import imaspy

from vtkggdtools.convert import ggd_to_vtk
from vtkggdtools.tests.fill_ggd import fill_ids


class Convert:

    def setup(self):
        self.edge_profiles = imaspy.IDSFactory().edge_profiles()
        fill_ids(self.edge_profiles, N=25)

    def time_convert(self):
        ggd_to_vtk(self.edge_profiles, 0)
