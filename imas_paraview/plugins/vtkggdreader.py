import logging

from paraview.util.vtkAlgorithm import smhint, smproxy

from imas_paraview.paraview_support.servermanager_tools import checkbox, propertygroup
from imas_paraview.plugins.ggd_base_reader import GGDBaseReader

logger = logging.getLogger("imas_paraview")


SUPPORTED_IDS_NAMES = [
    "edge_profiles",
    "edge_sources",
    "edge_transport",
    "mhd",
    "radiation",
    "runaway_electrons",
    "wall",
    "plasma_profiles",
    "plasma_sources",
    "plasma_transport",
]

# FIXME: Some IDSs do not have the grid structure in a separate `grid_ggd` object, as
# described by the GGD guidelines. They are for now denoted as experimental IDS, as they
# are currently not covered by the unit tests. If the DD for these stays like this, they
# will need to be handled separately. GGD grid locations for each IDS:
# equilibrium (grids_ggd/grid)
# distribution_sources (source/ggd/grid)
# distributions (distribution/ggd/grid)
# tf (field_map/grid)
# transport_solver_numerics (boundary_conditions_ggd/grid)
# waves (coherent_wave/full_wave/grid)

EXPERIMENTAL_IDS_NAMES = [
    "equilibrium",
    "distribution_sources",
    "distributions",
    "tf",
    "transport_solver_numerics",
    "waves",
]


@smproxy.source(label="GGD Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class GGDReader(GGDBaseReader, is_time_dependent=True):
    """ParaView plugin to view (spatial) GGD structures grids and associated data."""

    _show_parent_indices_dropdown = True

    def __init__(self):
        supported_ids = EXPERIMENTAL_IDS_NAMES + SUPPORTED_IDS_NAMES
        super().__init__(supported_ids)
        self._MIN_N_PLANE = self._n_plane = 0

    @propertygroup("Coordinate Transformation", ["ConvertCylToCart"])
    def PG4_CoordinateTransformation(self):
        """Dummy function to define a PropertyGroup."""

    @checkbox(
        name="ConvertCylToCart",
        label="Convert cylindrical to Cartesian components",
        default_values="0",
    )
    def P15_SetConvertCylToCart(self, val):
        """When enabled, cylindrical vector components (r, phi, z) are
        automatically converted to Cartesian (x, y, z) components. Disable to
        keep the original cylindrical components."""
        self.convert_cyl_to_cart = bool(val)
        self.Modified()
