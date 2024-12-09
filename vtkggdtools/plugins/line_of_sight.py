"""IMASPy plugin to view line of sight IDS structures in Paraview."""

import logging
from dataclasses import dataclass

from imaspy.ids_structure import IDSStructure
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from vtkggdtools.ids_util import get_object_by_name
from vtkggdtools.plugins.base_class import GGDVTKPluginBase
from vtkggdtools.util import points_to_vtkpoly, pol_to_cart

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
    """Data class that stores ``line_of_sight`` IDS structures, along with its name."""

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

    def _create_vtk_los(self, channel):
        """Create a vtkPolyData containing a line, based on the r, phi, and z
        coordinates in the line of sight structures. The r,phi,z-coordinates are
        converted to cartesian and stored as vtkPoints and connected using vtkLines,
        which both stored in a vtkPolyData object.

        Args:
            channel containing a line_of_sight structure

        Returns:
            vtkPolyData containing line_of_sight data
        """
        los = channel.line_of_sight

        points = [los.first_point, los.second_point]
        if hasattr(los, "third_point"):
            third_point = los.third_point
            if (
                third_point.r.has_value
                and third_point.phi.has_value
                and third_point.z.has_value
            ):
                points.append(third_point)

        converted_points = [
            (*pol_to_cart(point.r, point.phi), point.z) for point in points
        ]
        vtk_poly = points_to_vtkpoly(converted_points)
        return vtk_poly
