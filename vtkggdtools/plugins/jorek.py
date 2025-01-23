"""IMASPy plugin to view Jorek GGD data."""

import logging

from paraview.util.vtkAlgorithm import smhint, smproxy

from vtkggdtools.plugins.ggd_base_reader import GGDBaseReader

logger = logging.getLogger("vtkggdtools")


SUPPORTED_IDS_NAMES = [
    "mhd",
    "radiation",
]


@smproxy.source(label="JOREK Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class JorekGGDReader(
    GGDBaseReader,
    use_bezier=True,
    is_time_dependent=True,
):
    def __init__(self):
        supported_ids = SUPPORTED_IDS_NAMES
        super().__init__(supported_ids)
        self._n_plane = 1
        self._phi_start = 0
        self._phi_end = 0
