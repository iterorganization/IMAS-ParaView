"""Plugin to view line of sight IDS structures in Paraview."""

import logging
from dataclasses import dataclass

from imas.ids_structure import IDSStructure
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkLine,
    vtkMultiBlockDataSet,
    vtkPolyData,
)

from imas_paraview.ids_util import get_object_by_name
from imas_paraview.paraview_support.servermanager_tools import (
    doublevector,
    propertygroup,
)
from imas_paraview.plugins.base_class import GGDVTKPluginBase
from imas_paraview.util import pol_to_cart

logger = logging.getLogger("imas_paraview")

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
    """Data class that stores ``line_of_sight`` IDS structures, along with its name."""

    name: str
    line_of_sight: IDSStructure


@smproxy.source(label="Line of Sight Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class LineOfSightReader(GGDVTKPluginBase):
    """Line of sight reader based on IMAS-Python"""

    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", SUPPORTED_LINE_OF_SIGHT_IDS)
        self.scaling_factor = 1

    @doublevector(label="Scaling Factor", name="scaling_factor", default_values=1)
    def P99_SetScalingFactor(self, val):
        """Sets the scaling factor for the line of sight's length."""
        self._update_property("scaling_factor", val)

    @propertygroup("Line of Sight Settings", ["scaling_factor"])
    def PG3_LineOfSightGroup(self):
        """Dummy function to define a PropertyGroup."""

    def GetAttributeArrayName(self, idx) -> str:
        return self._selectable[idx].name

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        if len(self._selected) > 0:
            output = vtkMultiBlockDataSet.GetData(outInfo)
            self._load_los(output)
        return 1

    def setup_ids(self):
        """
        Select which line_of_sights to show in the array domain selector. The
        channel AoS is search through for line_of_sight structures. Their names are
        generated based on the channel name, and added to the array domain selector.
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

            # For ece, if the channel does not have a line of sight, or it is empty,
            # use the "global" line_of_sight instead
            if self._ids.metadata.name == "ece" and (
                not hasattr(channel, "line_of_sight")
                or not channel.line_of_sight.first_point.r.has_value
            ):
                selectable = LineOfSight(str(channel_name), self._ids.line_of_sight)
            else:
                selectable = LineOfSight(str(channel_name), channel.line_of_sight)
            self._selectable.append(selectable)

    def _load_los(self, output):
        """Go through the list of selected line_of_sights, and load each of them in a
        separate vtkPolyData object, which are all combined into a vtkMultiBlockDataSet.

        Args:
            output: The vtkMultiBlockDataSet containing the line_of_sights.
        """
        for i, channel_name in enumerate(self._selected):
            channel = get_object_by_name(self._selectable, channel_name)
            if channel is None:
                raise ValueError(f"Could not find {channel_name}")

            logger.info(f"Selected {channel.name}")
            vtk_poly = self._create_vtk_los(channel)
            output.SetBlock(i, vtk_poly)

    def _extend_line(self, first_point, second_point):
        """Extends a line in the direction from the first_point to the second_point. The
        distance from the first_point to the new point is determined by the original
        distance between the first and second point, scaled by the scaling factor.

        Args:
            first_point: Tuple containing the (x,y,z) coordinates of the first point.
            second_point: Tuple containing the (x,y,z) coordinates of the second point.

        Returns:
            Tuple containing a point (x,y,z) along the line from the first to the second
            point, at a distance from the first point defined by the scaling factor.
        """
        line1_direction = tuple(x - y for x, y in zip(second_point, first_point))
        return tuple(
            x + self.scaling_factor * y for x, y in zip(first_point, line1_direction)
        )

    def _create_vtk_los(self, channel):
        """Create a vtkPolyData containing a line, based on the r, phi, and z
        coordinates in the line of sight structures. The r-, phi-, and z-coordinates are
        converted to cartesian and stored as vtkPoints and connected using vtkLines,
        which are both stored in a vtkPolyData object. The lines can be scaled by a
        scaling factor.

        The connecting lines are created based on the graph below. Here, the points
        0 -> 1 describe the main LoS (los.first_point -> los.second_point). The
        reflection is described by points 1 -> 2 (los.second_point -> los.third_point).
        The extensions are visualized by the plusses in the graph, from 0 -> 3 for the
        main LoS, and from 1 -> 4, for the reflection. Only the lines from 0 -> 3 (and
        1 -> 4, if there is a reflection) are transformed into vtkLines.

        0----1+++++3
              \
               \
                2
                 +
                  +
                   4
        Args:
            channel: containing a line_of_sight structure.

        Returns:
            vtkPolyData containing line_of_sight data.
        """
        los = channel.line_of_sight
        points = [None] * 5
        points[0] = los.first_point
        points[1] = los.second_point

        # Add reflection point
        has_reflection = False
        if hasattr(los, "third_point"):
            third_point = los.third_point
            if (
                third_point.r.has_value
                and third_point.phi.has_value
                and third_point.z.has_value
            ):
                has_reflection = True
                points[2] = third_point

        # Convert from cylindrical to cartesian coordinates
        points = [
            (*pol_to_cart(point.r, point.phi), point.z) if point is not None else None
            for point in points
        ]

        # Add line extension
        points[3] = self._extend_line(points[0], points[1])
        if has_reflection:
            points[4] = self._extend_line(points[1], points[2])

        vtk_poly = self._create_vtk_poly(points, has_reflection)
        return vtk_poly

    def _create_vtk_poly(self, points, has_reflection):
        """Create vtkPolyData containing the vtkPoints and vtkLines in a vtkCellArray.

        Args:
            points: List of points to store in the vtkPolyData.
            has_reflection: boolean whether LoS contains a reflection point.

        Returns:
            vtkPolyData containing vtkPoints and vtkCellArray containing the vtkLines
            corresponding to the scaled LoS and optionally its reflection.
        """
        vtk_points = vtkPoints()
        vtk_lines = vtkCellArray()

        vtk_points.InsertNextPoint(*points[0])
        vtk_points.InsertNextPoint(*points[3])
        line = vtkLine()
        line.GetPointIds().SetId(0, 0)
        line.GetPointIds().SetId(1, 1)
        vtk_lines.InsertNextCell(line)
        if has_reflection:
            vtk_points.InsertNextPoint(*points[1])
            vtk_points.InsertNextPoint(*points[4])
            line = vtkLine()
            line.GetPointIds().SetId(0, 2)
            line.GetPointIds().SetId(1, 3)
            vtk_lines.InsertNextCell(line)

        vtk_poly = vtkPolyData()
        vtk_poly.SetPoints(vtk_points)
        vtk_poly.SetLines(vtk_lines)
        return vtk_poly
