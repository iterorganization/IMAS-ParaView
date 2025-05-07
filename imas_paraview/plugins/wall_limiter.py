"""Plugin to view limiter structures of 2d wall descriptions nodes"""

import logging
from dataclasses import dataclass

from imas.ids_structure import IDSStructure
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from imas_paraview.ids_util import get_object_by_name
from imas_paraview.plugins.base_class import GGDVTKPluginBase
from imas_paraview.util import points_to_vtkpoly

logger = logging.getLogger("imas_paraview")


@dataclass
class Limiter:
    """Data class that stores limiter units, along with its name."""

    name: str
    unit: IDSStructure


@smproxy.source(label="Wall Limiter Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class WallLimiterReader(GGDVTKPluginBase):
    """IMAS-Python based reader that can load the limiter structures of 2D wall
    descriptions"""

    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", ["wall"])

    def GetAttributeArrayName(self, idx) -> str:
        return self._selectable[idx].name

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        if len(self._selected) > 0:
            output = vtkMultiBlockDataSet.GetData(outInfo)
            self._load_limiters(output)
        return 1

    def setup_ids(self):
        """
        Select which limiters to show in the array domain selector. The description_2d
        AoS is searched through for limiter structures. Their names are generated based
        on the limiter type name and the limiter unit name, and added to the array
        domain selector.
        """
        assert self._ids is not None, "IDS cannot be empty during setup."

        if not self._ids.description_2d.has_value:
            raise ValueError("This IDS does not contain a description_2d node")

        descriptions = self._ids.description_2d
        self._selectable = []

        for i, description in enumerate(descriptions):
            description_name = description.type.name
            if description_name == "":
                description_name = f"description {i}"
                logger.warning(
                    "Found a limiter without a description name, "
                    f"it will be loaded as {description_name}"
                )

            limiter = description.limiter
            type_name = limiter.type.name
            if type_name == "":
                type_name = f"type {i}"
                logger.warning(
                    "Found a limiter without a type name, "
                    f"it will be loaded as {type_name}"
                )

            for j, unit in enumerate(limiter.unit):
                unit_name = unit.name
                if unit_name == "":
                    unit_name = f"unit {j}"
                    logger.warning(
                        "Found a limiter unit without a name, "
                        f"it will be loaded as {unit_name}"
                    )
                name = " / ".join(
                    [str(description_name), str(type_name), str(unit_name)]
                )
                selectable = Limiter(name, unit)
                self._selectable.append(selectable)

    def _load_limiters(self, output):
        """Go through the list of selected limiters, and load each of them in a
        separate vtkPolyData object, which are all combined into a vtkMultiBlockDataSet.

        Args:
            output: The vtkMultiBlockDataSet containing the limiter contours.
        """
        for i, limiter_name in enumerate(self._selected):
            limiter = get_object_by_name(self._selectable, limiter_name)
            if limiter is None:
                raise ValueError(f"Could not find {limiter_name}")

            logger.info(f"Selected {limiter.name}")
            vtk_poly = self._create_contour(limiter)
            output.SetBlock(i, vtk_poly)

    def _create_contour(self, limiter):
        """Create a contour based on the r,z coordinates in the limiter.
        The r,z-coordinates are stored as vtkPoints, and connected using vtkLines, which
        are both stored in a vtkPolyData object. If the contour is closed, the start
        and end points are connected.

        Args:
            limiter: limiter object containing contour

        Returns:
            vtkPolyData containing contour data
        """
        # closed was removed in DD4 - data providers need to repeat the first point
        # for closed outlines, which Just Works with our is_closed=False logic
        if getattr(limiter.unit, "closed", 0) == 0:
            is_closed = False
        else:
            is_closed = True

        r = limiter.unit.outline.r
        z = limiter.unit.outline.z
        points = [(ri, 0.0, zi) for ri, zi in zip(r, z)]

        assert len(r) == len(z), "r and z must have the same length."

        vtk_poly = points_to_vtkpoly(points, is_closed)
        return vtk_poly
