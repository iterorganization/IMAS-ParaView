"""ParaView plugin to visualize pellets"""

import logging
from dataclasses import dataclass

import numpy as np
from imas.ids_primitive import IDSNumericArray
from imas.ids_struct_array import IDSStructArray
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet
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
    fragment: IDSStructArray  # spi.injector.fragment


@dataclass
class VelocityMassCentre:
    vel_r: IDSNumericArray  # spi.injector.velocity_mass_centre_fragments_r
    vel_phi: IDSNumericArray  # spi.injector.velocity_mass_centre_fragments_phi
    vel_z: IDSNumericArray  # spi.injector.velocity_mass_centre_fragments_z
    origin_r: IDSNumericArray  # spi.injector.shatter_cone.origin.r
    origin_phi: IDSNumericArray  # spi.injector.shatter_cone.origin.phi
    origin_z: IDSNumericArray  # spi.injector.shatter_cone.origin.z


logger = logging.getLogger("imas_paraview")

# TODO: Currently only the shattered pellet fragments are visualized. This plugin
# can be expanded to also visualize the following:
# - pellet trajectory path in 'pellets' IDS
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
        """Scaling factor for the size of the velocity arrows of shattered pellet 
        fragments"""
        self.frag_scaling_factor = 1.0
        """Scaling factor for the size of spheres of shattered pellet fragments"""
        self.com_vel_scaling_factor = 1.0
        """Scaling factor for the size of the velocity arrows of the centre of mass of
        the fragments"""

    @doublevector(
        label="Fragment Sphere Scaling Factor",
        name="frag_scaling_factor",
        default_values=1.0,
    )
    def P97_SetFragmentScalingFactor(self, val):
        self._update_property("frag_scaling_factor", val)

    @doublevector(
        label="Fragment Velocity Scaling Factor",
        name="frag_vel_scaling_factor",
        default_values=1.0,
    )
    def P98_SetFragmentVelocityScalingFactor(self, val):
        self._update_property("frag_vel_scaling_factor", val)

    @doublevector(
        label="Velocity Centre of Mass Scaling Factor",
        name="frag_vel_scaling_factor",
        default_values=1.0,
    )
    def P99_SetCoMVelocityScalingFactor(self, val):
        self._update_property("com_vel_scaling_factor", val)

    @propertygroup(
        "Pellets Reader Settings",
        ["frag_scaling_factor", "frag_vel_scaling_factor", ""],
    )
    def PG3_PelletReaderGroup(self):
        """Dummy function to define a PropertyGroup."""

    def setup_ids(self):
        assert self._ids is not None, "IDS cannot be empty during setup."

        self.selectable_map = {}

        for i, injector in enumerate(self._ids.injector):
            injector_name = injector.name
            if not injector_name:
                injector_name = f"injector {i}"
                logger.warning(
                    "Found an injector without a name, loading as '%s'.", injector_name
                )
            self._populate_fragments(injector, injector_name)
            self._populate_vel_mass_centre(injector, injector_name)
        self._selectable = list(self.selectable_map.keys())

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
            if isinstance(selected, Fragments):
                vtk_object = self._create_fragments_geom(selected, time_idx)
            else:  # selected is a VelocityMassCentre
                vtk_object = self._build_vel_mass_centre_geom(selected)

            output.SetBlock(block_id, vtk_object)
        return 1

    def _populate_fragments(self, injector, injector_name):
        if len(injector.fragment) > 0:
            injector_name = f"Shattered fragments ({injector_name})"
            self.selectable_map[injector_name] = Fragments(injector.fragment)
        else:
            logger.warning("'%s' has no fragments, skipping.", injector_name)

    def _populate_vel_mass_centre(self, injector, injector_name):
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
            self.selectable_map[injector_name] = VelocityMassCentre(
                vel_r, vel_phi, vel_z, origin.r, origin.phi, origin.z
            )
        else:
            logger.warning(
                "'%s' does not have both the 'velocity_mass_centre' and "
                "'shatter_cone.origin' defined. The velocity centre of mass of this "
                "injector will not be loaded.",
                injector_name,
            )

    def _build_vel_mass_centre_geom(self, vel_mass_centre: VelocityMassCentre):
        vel_r = vel_mass_centre.vel_r
        vel_phi = vel_mass_centre.vel_phi
        vel_z = vel_mass_centre.vel_z
        pos_r = vel_mass_centre.origin_r
        pos_phi = vel_mass_centre.origin_phi
        pos_z = vel_mass_centre.origin_z

        pos_x, pos_y = pol_to_cart(pos_r, pos_phi)
        position = np.array([[pos_x, pos_y, pos_z]])
        velocity = vel_pol_to_cart(
            np.array([vel_r]),
            np.array([vel_phi]),
            np.array([vel_z]),
            np.array([pos_phi]),
        )
        return create_vtk_arrows(position, velocity)

    def _create_fragments_geom(self, frags: Fragments, time_idx):
        fragments = frags.fragment
        num_fragments = len(fragments)

        r = np.empty(num_fragments)
        phi = np.empty(num_fragments)
        z = np.empty(num_fragments)
        volumes = np.empty(num_fragments)
        vel_r = np.empty(num_fragments)
        vel_phi = np.empty(num_fragments)
        vel_z = np.empty(num_fragments)

        for i, f in enumerate(fragments):
            position = f.position
            r[i] = position.r[time_idx]
            phi[i] = position.phi[time_idx]
            z[i] = position.z[time_idx]
            volumes[i] = f.volume[time_idx]
            vel_r[i] = f.velocity_r[time_idx]
            try:
                vel_phi[i] = f.velocity_phi[time_idx]
            except AttributeError:
                vel_phi[i] = f.velocity_tor[time_idx]
            vel_z[i] = f.velocity_z[time_idx]

        pos_x, pos_y = pol_to_cart(r, phi)
        positions = np.stack([pos_x, pos_y, z], axis=1)

        # Volume to radius conversion of spheres
        radii = (3.0 * volumes / (4.0 * np.pi)) ** (1.0 / 3.0)
        velocities = vel_pol_to_cart(vel_r, vel_phi, vel_z, phi)

        spheres = create_vtk_spheres(
            positions, radii, scaling_factor=self.frag_scaling_factor
        )
        arrows = create_vtk_arrows(
            positions, velocities, scaling_factor=self.frag_vel_scaling_factor
        )

        vtk_fragments = vtkAppendPolyData()
        vtk_fragments.AddInputData(spheres)
        vtk_fragments.AddInputData(arrows)
        vtk_fragments.Update()
        return vtk_fragments.GetOutput()
