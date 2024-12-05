"""IMASPy plugin to view beam structures of ec_launchers IDS in Paraview."""

import logging
from dataclasses import dataclass

import numpy as np
from imaspy.ids_structure import IDSStructure
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkLine,
    vtkMultiBlockDataSet,
    vtkPolyData,
)

from vtkggdtools.ids_util import get_object_by_name
from vtkggdtools.plugins.base_class import GGDVTKPluginBase
from vtkggdtools.util import pol_to_cart

logger = logging.getLogger("vtkggdtools")


@dataclass
class Beam:
    """Data class that stores ``beam`` IDS structures, along with its name."""

    name: str
    beam: IDSStructure


@smproxy.source(label="Beam Reader")
@smhint.xml("""<ShowInMenu category="VTKGGDTools" />""")
class IMASPyBeamReader(GGDVTKPluginBase):
    """Beam reader based on IMASPy"""

    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", ["ec_launchers"])

    def GetAttributeArrayName(self, idx) -> str:
        return self._selectable[idx].name

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        if len(self._selected) > 0:
            output = vtkMultiBlockDataSet.GetData(outInfo)
            self._load_los(output)
        return 1

    def request_information(self):
        """
        Abstract method which is called during the RequestInformation pipeline step.
        Intentionally left empty.
        """

    def setup_ids(self):
        """
        Select which line_of_sights to show in the array domain selector. The
        channel AoS is search through for line_of_sight structures. Their names are
        generated based on the channel name, and added to the array domain selector.
        """
        assert self._ids is not None, "IDS cannot be empty during setup."

        self._selectable = []

        beams = self._ids.beam

        for i, beam in enumerate(beams):
            beam_name = beam.name
            if beam_name == "":
                beam_name = f"beam {i}"
                logger.warning(
                    "Found a channel without a name, "
                    f"it will be loaded as {beam_name}."
                )

            selectable = Beam(str(beam_name), beam)
            self._selectable.append(selectable)

    def _load_los(self, output):
        """Go through the list of selected line_of_sights, and load each of them in a
        separate vtkPolyData object, which are all combined into a vtkMultiBlockDataSet.

        Args:
            output: The vtkMultiBlockDataSet containing the line_of_sights.
        """
        for i, beam_name in enumerate(self._selected):
            beam = get_object_by_name(self._selectable, beam_name)
            if beam is None:
                raise ValueError(f"Could not find {beam_name}")

            logger.info(f"Selected {beam.name}")
            vtk_poly = self._create_vtk_beam(beam.beam)
            output.SetBlock(i, vtk_poly)

    def _create_vtk_beam(self, beam):
        """Create a vtkPolyData containing a line, based on the r, phi, and z
        coordinates in the line of sight structures. The r,phi,z-coordinates are
        converted to cartesian and stored as vtkPoints and connected using vtkLines,
        which both stored in a vtkPolyData object.

        Args:
            channel containing a line_of_sight structure

        Returns:
            vtkPolyData containing line_of_sight data
        """

        first_point, second_point = get_points(
            beam.launching_position,
            beam.steering_angle_pol,
            beam.steering_angle_tor,
            10,
        )

        vtk_points = vtkPoints()
        vtk_points.InsertNextPoint(*first_point)
        vtk_points.InsertNextPoint(*second_point)
        vtk_lines = vtkCellArray()
        line = vtkLine()
        line.GetPointIds().SetId(0, 0)
        line.GetPointIds().SetId(1, 1)
        vtk_lines.InsertNextCell(line)
        vtk_poly = vtkPolyData()
        vtk_poly.SetPoints(vtk_points)
        vtk_poly.SetLines(vtk_lines)
        return vtk_poly


def get_points(launch_pos, steering_angle_pol, steering_angle_tor, distance):
    """Apply the steering angles to the beam's position and calculate a point further along the direction.

    Args:
        r: Radial distance of the beam's position
        z: Z height of the beam's position
        phi: Azimuthal angle of the beam's position
        steering_angle_pol: Steering angle in the R,Z plane (radians)
        steering_angle_tor: Steering angle away from the poloidal plane (radians)
        distance: The distance along the beam's direction to the new point

    Returns:
        Adjusted Cartesian coordinates of the beam and the new position further along the direction.
    """
    first_point = (*pol_to_cart(launch_pos.r[0], launch_pos.phi[0]), launch_pos.z[0])

    x_pol = launch_pos.r[0] * np.cos(steering_angle_pol)
    y_pol = launch_pos.r[0] * np.sin(steering_angle_pol)

    x_tor = np.sin(steering_angle_tor)
    y_tor = 0
    z_tor = np.cos(steering_angle_tor)

    direction_x = x_pol + x_tor
    direction_y = y_pol + y_tor
    direction_z = z_tor

    x_end = first_point[0] + distance * direction_x
    y_end = first_point[1] + distance * direction_y
    z_end = first_point[2] + distance * direction_z

    return first_point, (x_end, y_end, z_end)
