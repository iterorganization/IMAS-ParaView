"""Plugin to view Jorek GGD data."""

import logging

from paraview.util.vtkAlgorithm import smhint, smproxy

from imas_paraview.plugins.ggd_base_reader import GGDBaseReader

logger = logging.getLogger("imas_paraview")


SUPPORTED_IDS_NAMES = [
    "mhd",
    "radiation",
    "plasma_profiles",
]


@smproxy.source(label="JOREK Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class JorekGGDReader(
    GGDBaseReader,
    use_bezier=True,
    is_time_dependent=True,
):
    """ParaView plugin to view JOREK Grids in ParaView.

    JOREK uses Cubic Bézier elements in the poloidal plane and Fourier interpolation in
    the toroidal direction, which requires its own plugin.
    """

    def __init__(self):
        supported_ids = SUPPORTED_IDS_NAMES
        super().__init__(supported_ids)
        self._MIN_N_PLANE = self._n_plane = 1
