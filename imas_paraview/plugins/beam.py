"""Plugin to view beam structures of ec_launchers IDS in Paraview."""

import logging
from dataclasses import dataclass

import numpy as np
from imas.ids_structure import IDSStructure
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from imas_paraview.ids_util import get_object_by_name
from imas_paraview.paraview_support.servermanager_tools import (
    doublevector,
    propertygroup,
)
from imas_paraview.plugins.base_class import GGDVTKPluginBase
from imas_paraview.util import find_closest_indices, points_to_vtkpoly, pol_to_cart

logger = logging.getLogger("imas_paraview")


@dataclass
class Beam:
    """Data class that stores ``beam`` IDS structures, along with its name."""

    name: str
    beam: IDSStructure


@smproxy.source(label="Beam Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class BeamReader(GGDVTKPluginBase, is_time_dependent=True):
    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", ["ec_launchers"])
        self.distance = 10

    @doublevector(label="Beam Distance", name="beam_distance", default_values=10)
    def P99_SetDistance(self, val):
        """Sets the beam's distance from the launching position."""
        self._update_property("distance", val)

    @propertygroup("Beam settings", ["beam_distance"])
    def PG3_BeamGroup(self):
        """Dummy function to define a PropertyGroup."""

    def GetAttributeArrayName(self, idx) -> str:
        return self._selectable[idx].name

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        # Retrieve the selected time step and profiles
        time = self._get_selected_time_step(outInfo)
        if time is None:
            logger.warning("Selected invalid time step")
            return 1

        index_list = find_closest_indices([time], self._ids.time)
        time_idx = index_list[0]

        if len(self._selected) > 0:
            output = vtkMultiBlockDataSet.GetData(outInfo)
            self._load_beam(output, time_idx)
        return 1

    def setup_ids(self):
        """
        Select which beams to show in the array domain selector. The IDS is searched
        for beam structures and their names are added to the array domain selector.
        """
        assert self._ids is not None, "IDS cannot be empty during setup."

        self._selectable = []

        beams = self._ids.beam

        for i, beam in enumerate(beams):
            beam_name = beam.name
            if beam_name == "":
                beam_name = f"beam {i}"
                logger.warning(
                    f"Found a channel without a name, it will be loaded as {beam_name}."
                )

            selectable = Beam(str(beam_name), beam)
            self._selectable.append(selectable)

    def _load_beam(self, output, time_idx):
        """Go through the list of selected beams, and load each of them in a
        separate vtkPolyData object, which are all combined into a vtkMultiBlockDataSet.

        Args:
            output: The vtkMultiBlockDataSet containing the beams.
        """
        for i, beam_name in enumerate(self._selected):
            beam = get_object_by_name(self._selectable, beam_name)
            if beam is None:
                raise ValueError(f"Could not find {beam_name}")

            logger.info(f"Selected {beam.name}")
            vtk_poly = self._create_vtk_beam(beam.beam, time_idx)
            output.SetBlock(i, vtk_poly)

    def _create_vtk_beam(self, beam, time_idx):
        """Create a vtkPolyData containing two points from the beam's
        launching_position to a point in the direction that the steering angles point
        to, at a distance from the launching position. For this, the points and
        direction are converted to cartesian coordinates, and a vtkLine is created
        between these two points.

        Args:
            beam: the beam structure to create.

        Returns:
            vtkPolyData containing the beams
        """

        first_point, second_point = self.transform_points(
            beam.launching_position.r[time_idx],
            beam.launching_position.phi[time_idx],
            beam.launching_position.z[time_idx],
            beam.steering_angle_pol[time_idx],
            beam.steering_angle_tor[time_idx],
        )

        points = [first_point, second_point]
        vtk_poly = points_to_vtkpoly(points)
        return vtk_poly

    def transform_points(
        self,
        launch_pos_r,
        launch_pos_phi,
        launch_pos_z,
        steering_angle_pol,
        steering_angle_tor,
    ):
        """Transform the launching position and launching direction to cartesian
        coordinates.

        Args:
            launch_pos_r: Radial distance of the launching position
            launch_pos_phi: Azimuth angle of the launching position
            launch_pos_z: Height of the launching position
            steering_angle_pol: Steering angle of the beam in the R,Z plane
            steering_angle_tor: Steering angle of the beam away from the poloidal plane
            distance: Distance along direction vector at which to place the second point

        Returns:
            tuple containing the launching position and a point into the direction of
            the launch, in cartesian coordinates: ((x1,y1,z1),(x2,y2,z2))
        """

        first_point = (*pol_to_cart(launch_pos_r, launch_pos_phi), launch_pos_z)

        # Calculate radial components of the wave vector (taking |k|=1)
        k_r = -np.cos(steering_angle_pol) * np.cos(steering_angle_tor)
        k_phi = np.cos(steering_angle_pol) * np.sin(steering_angle_tor)

        # Convert to direction in cartesian coordinates
        k_x = k_r * np.cos(launch_pos_phi) - k_phi * np.sin(launch_pos_phi)
        k_y = k_r * np.sin(launch_pos_phi) + k_phi * np.cos(launch_pos_phi)
        k_z = -np.sin(steering_angle_pol)

        # Calculate the end point from the launching position into the direction of the
        # steering angles
        second_point = (
            first_point[0] + self.distance * k_x,
            first_point[1] + self.distance * k_y,
            first_point[2] + self.distance * k_z,
        )

        return first_point, second_point
