"""IMASPy plugin to view line of sight nodes."""

import logging
from dataclasses import dataclass

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

SUPPORTED_LINE_OF_SIGHT_IDS = [
    "bolometer",
    "bremsstrahlung_visible",
    "ece",
    "hard_x_rays",
    "interferometer",
    "mse",
    "polarimeter",
    "refractometer",
    "soft_x_rays",
    "spectrometer_uv",
    "spectrometer_visible",
]


@dataclass
class LineOfSight:
    """Data class that stores 1d profiles, along with its name and coordinate array."""

    name: str
    line_of_sight: IDSStructure


@smproxy.source(label="line_of_sight Reader")
@smhint.xml("""<ShowInMenu category="VTKGGDTools" />""")
class IMASPyLineOfSightReader(GGDVTKPluginBase):
    """Line of sight reader based on IMASPy"""

    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", SUPPORTED_LINE_OF_SIGHT_IDS)

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
        Select which limiters to show in the array domain selector. The description_2d
        AoS is searched through for limiter structures. Their names are generated based
        on the limiter type name and the limiter unit name, and added to the array
        domain selector.
        """
        assert self._ids is not None, "IDS cannot be empty during setup."

        self._selectable = []

        channels = self._ids.channel

        for i, channel in enumerate(channels):
            channel_name = channel.name
            if channel_name == "":
                channel_name = f"channel {i}"
                logger.warning(
                    "Found a channel without a name, "
                    f"it will be loaded as {channel_name}."
                )

            # If the channel does not have a line of sight, or it is empty, use the
            # "global" line of sight instead
            if self._ids.metadata.name == "ece" and (
                not hasattr(channel, "line_of_sight")
                or not channel.line_of_sight.first_point.r.has_value
            ):
                selectable = LineOfSight(str(channel_name), self._ids.line_of_sight)
            else:
                selectable = LineOfSight(str(channel_name), channel.line_of_sight)
            self._selectable.append(selectable)

    def _load_los(self, output):
        """Go through the list of selected limiters, and load each of them in a
        separate vtkPolyData object, which are all combined into a vtkMultiBlockDataSet.

        Args:
            output: The vtkMultiBlockDataSet containing the limiter contours.
        """
        for i, channel_name in enumerate(self._selected):
            channel = get_object_by_name(self._selectable, channel_name)
            if channel is None:
                raise ValueError(f"Could not find {channel_name}")

            logger.info(f"Selected {channel.name}")
            vtk_poly = self._create_vtk_los(channel)
            output.SetBlock(i, vtk_poly)

    def _create_vtk_los(self, channel):
        """Create a contour based on the r,z coordinates in the limiter.
        The r,z-coordinates are stored as vtkPoints, and connected using vtkLines, which
        are both stored in a vtkPolyData object. If the contour is closed, the start
        and end points are connected.

        Args:
            limiter: limiter object containing contour

        Returns:
            vtkPolyData containing contour data
        """
        los = channel.line_of_sight

        first_point = los.first_point
        second_point = los.second_point

        first_point_cart = (*pol_to_cart(first_point.r, first_point.phi), first_point.z)
        second_point_cart = (
            *pol_to_cart(second_point.r, second_point.phi),
            second_point.z,
        )

        # Fill points
        vtk_points = vtkPoints()
        vtk_points.InsertNextPoint(
            first_point_cart[0], first_point_cart[1], float(first_point_cart[2])
        )
        vtk_points.InsertNextPoint(
            second_point_cart[0], second_point_cart[1], float(second_point_cart[2])
        )

        # Fill lines
        vtk_lines = vtkCellArray()
        line = vtkLine()
        line.GetPointIds().SetId(0, 0)
        line.GetPointIds().SetId(1, 1)
        vtk_lines.InsertNextCell(line)

        if hasattr(los, "third_point"):
            third_point = los.third_point
            if (
                third_point.r.has_value
                and third_point.phi.has_value
                and third_point.z.has_value
            ):
                third_point_cart = (
                    *pol_to_cart(third_point.r, third_point.phi),
                    third_point.z,
                )
                vtk_points.InsertNextPoint(*third_point_cart)
                line.GetPointIds().SetId(0, 1)
                line.GetPointIds().SetId(1, 2)
                vtk_lines.InsertNextCell(line)

        vtk_poly = vtkPolyData()
        vtk_poly.SetPoints(vtk_points)
        vtk_poly.SetLines(vtk_lines)
        return vtk_poly
