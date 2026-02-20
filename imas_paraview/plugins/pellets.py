"""ParaView plugin to visualize SPI fragment positions"""

import logging
from dataclasses import dataclass

import numpy as np
from imas.ids_primitive import IDSNumericArray
from imas.ids_struct_array import IDSStructArray
from imas.ids_structure import IDSStructure
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet, vtkPolyData
from vtkmodules.vtkFiltersCore import vtkAppendPolyData

from imas_paraview.paraview_support.servermanager_tools import (
    doublevector,
    propertygroup,
)
from imas_paraview.plugins.base_class import GGDVTKPluginBase
from imas_paraview.util import (
    create_vtk_arrows,
    create_vtk_spheres,
    find_closest_indices,
    pol_to_cart,
    vel_pol_to_cart,
)


@dataclass
class Fragments:
    fragments: IDSStructArray  # injector.fragment


@dataclass
class VelMassCentres:
    vel_r: IDSNumericArray  # injector.velocity_mass_centre_fragments_r
    vel_phi: IDSNumericArray  # injector.velocity_mass_centre_fragments_phi
    vel_z: IDSNumericArray  # injector.velocity_mass_centre_fragments_z
    shatter_cone_origin: IDSStructure  # injector.shatter_cone.origin


logger = logging.getLogger("imas_paraview")

# TODO: Also visualize 'pellets' IDS
# TODO: Currently only the shattered pellet fragments are visualized, we can also
# - spi.injector(i1).pellet position/velocity over time
# - injector(i1).injection_direction
# - injector(i1).shatter_cone
# - injector(i1)/pellet/core/species(i2).name

SUPPORTED_IDS_NAMES = ["spi"]


@smproxy.source(label="(Shattered) Pellets Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class PelletReader(GGDVTKPluginBase, is_time_dependent=True):
    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", SUPPORTED_IDS_NAMES)
        self.selectable_map = {}
        self.frag_vel_scaling_factor = 1.0
        """Scaling factor for the velocity arrow glyphs of shattered pellet fragments"""
        self.frag_scaling_factor = 1.0
        """Scaling factor for the spheres of shattered pellet fragments"""

    @doublevector(
        label="VelocityScalingFactor",
        name="frag_vel_scaling_factor",
        default_values=1.0,
    )
    def P98_SetVelocityScalingFactor(self, val):
        self._update_property("frag_vel_scaling_factor", val)

    @doublevector(label="ScalingFactor", name="frag_scaling_factor", default_values=1.0)
    def P99_SetScalingFactor(self, val):
        self._update_property("frag_scaling_factor", val)

    @propertygroup(
        "Pellets Reader Settings", ["frag_scaling_factor", "frag_vel_scaling_factor"]
    )
    def PG3_PelletReaderGroup(self):
        """Dummy function to define a PropertyGroup."""

    def setup_ids(self):
        """Populate the array selector with injector names."""
        assert self._ids is not None, "IDS cannot be empty during setup."

        self.selectable_map = {}

        for i, injector in enumerate(self._ids.injector):
            injector_name = injector.name
            if not injector_name:
                injector_name = f"injector {i}"
                logger.warning(
                    "Found an injector without a name, loading as '%s'.", injector_name
                )
            self.populate_fragments(injector, injector_name)
            self.populate_vel_mass_centre(injector, injector_name)
        self._selectable = list(self.selectable_map.keys())

    def populate_fragments(self, injector, injector_name):
        if len(injector.fragment) > 0:
            injector_name = f"Shattered fragments ({injector_name})"
            self.selectable_map[injector_name] = Fragments(fragments=injector.fragment)
        else:
            logger.warning("'%s' has no fragments, skipping.", injector_name)

    def populate_vel_mass_centre(self, injector, injector_name):
        vel_r = injector.velocity_mass_centre_fragments_r

        try:
            vel_phi = injector.velocity_mass_centre_fragments_phi
        except AttributeError:
            vel_phi = injector.velocity_mass_centre_fragments_tor
        vel_z = injector.velocity_mass_centre_fragments_z
        origin = injector.shatter_cone.origin

        if (
            vel_r.has_value
            and vel_phi.has_value
            and vel_z.has_value
            and origin.r.has_value
            and origin.phi.has_value
            and origin.z.has_value
        ):
            injector_name = f"Fragment centre of mass velocity ({injector_name})"
            self.selectable_map[injector_name] = VelMassCentres(
                vel_r=vel_r,
                vel_phi=vel_phi,
                vel_z=vel_z,
                shatter_cone_origin=origin,
            )
        else:
            logger.warning(
                "'%s' does not have both the 'velocity_mass_centre' and "
                "'shatter_cone.origin' defined. The velocity centre of mass of this "
                "injector will not be loaded.",
                injector_name,
            )

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        time = self._get_selected_time_step(outInfo)
        if time is None:
            return 1

        time_idx = find_closest_indices([time], self._ids.time)[0]
        output = vtkMultiBlockDataSet.GetData(outInfo)

        for block_id, selection_name in enumerate(self._selected):
            selected = self.selectable_map[selection_name]
            vtk_object = None
            if isinstance(selected, Fragments):
                spheres_poly, arrows_poly = self._create_fragments_geom(
                    selected, time_idx
                )
                appended = vtkAppendPolyData()
                appended.AddInputData(spheres_poly)
                appended.AddInputData(arrows_poly)
                appended.Update()
                vtk_object = appended.GetOutput()
            elif isinstance(selected, VelMassCentres):
                vtk_object = self._build_vel_mass_centre_geom(selected)

            if vtk_object is not None:
                output.GetMetaData(block_id).Set(selection_name)
                output.SetBlock(block_id, vtk_object)
        return 1

    def _build_vel_mass_centre_geom(self, selectable: VelMassCentres) -> vtkPolyData:
        try:
            vel_phi = selectable.vel_phi
        except AttributeError:
            vel_phi = selectable.vel_tor

        pos_phi = selectable.shattering_position.phi
        pos_x, pos_y = pol_to_cart(selectable.shattering_position.r, pos_phi)

        position = np.array([[pos_x, pos_y, selectable.shattering_position.z]])
        velocity = vel_pol_to_cart(
            np.array([selectable.vel_r]),
            np.array([vel_phi]),
            np.array([selectable.vel_z]),
            np.array([pos_phi]),
        )

        return create_vtk_arrows(position, velocity)

    def _create_fragments_geom(self, selectable, time_idx):
        """Create sphere and arrow polydata for all fragments of one injector.

        Args:
            injector: The injector IDS node.
            time_idx: The time index to query.

        Returns:
            Tuple of (spheres_polydata, arrows_polydata).
        """
        fragments = selectable.fragments
        num_fragments = len(fragments)

        r = np.empty(num_fragments)
        phi = np.empty(num_fragments)
        z = np.empty(num_fragments)
        volumes = np.empty(num_fragments)
        v_r = np.empty(num_fragments)
        v_phi = np.empty(num_fragments)
        v_z = np.empty(num_fragments)

        for i, f in enumerate(fragments):
            pos = f.position
            r[i] = pos.r[time_idx]
            phi[i] = pos.phi[time_idx]
            z[i] = pos.z[time_idx]
            volumes[i] = f.volume[time_idx]
            v_r[i] = f.velocity_r[time_idx]
            v_z[i] = f.velocity_z[time_idx]
            try:
                v_phi[i] = f.velocity_phi[time_idx]
            except AttributeError:
                v_phi[i] = f.velocity_tor[time_idx]

        pos_x, pos_y = pol_to_cart(r, phi)
        positions = np.stack([pos_x, pos_y, z], axis=1)
        radii = (3.0 * volumes / (4.0 * np.pi)) ** (1.0 / 3.0)
        velocities = vel_pol_to_cart(v_r, v_phi, v_z, phi)

        spheres = create_vtk_spheres(positions, radii)
        arrows = create_vtk_arrows(positions, velocities)
        return spheres, arrows
