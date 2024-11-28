"""IMASPy plugin to view wall description nodes"""

import logging
from dataclasses import dataclass

from imaspy.ids_structure import IDSStructure
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkDataObject,
    vtkLine,
    vtkMultiBlockDataSet,
    vtkPartitionedDataSetCollection,
    vtkPolyData,
    vtkPolyDataCollection,
)

from vtkggdtools.plugins.base_class import GGDVTKPluginBase

logger = logging.getLogger("vtkggdtools")


@dataclass
class Limiter:
    """Data class that stores limiter units, along with its name."""

    name: str
    unit: IDSStructure


@smproxy.source(label="description_2d Reader")
@smhint.xml("""<ShowInMenu category="VTKGGDTools" />""")
class IMASPyDescription2DReader(GGDVTKPluginBase):
    """Wall descriptions reader based on IMASPy"""

    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", ["wall"])

    def GetAttributeArrayName(self, idx) -> str:
        return self._selectable[idx].name

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        if len(self._selected) > 0:
            print(self._selected)
            output = vtkMultiBlockDataSet.GetData(outInfo)
            self._load_limiters(output)
        return 1

    def request_information(self):
        """
        AbstracSelect which profiles to show in the array domain selector, based
        on whether the "Show All" checkbox is enabled.
        """

    def setup_ids(self):
        """
        Select which profiles to show in the array domain selector, based
        on whether the "Show All" checkbox is enabled.
        """
        assert self._ids is not None, "IDS cannot be empty during setup."

        if not self._ids.description_2d.has_value:
            raise ValueError("This IDS does not contain a description_2d node")

        descriptions = self._ids.description_2d
        self._selectable = []

        for description in descriptions:
            limiter = description.limiter
            type_name = limiter.type.name
            for unit in limiter.unit:
                name = " / ".join([str(type_name), str(unit.name)])
                selectable = Limiter(name, unit)
                self._selectable.append(selectable)

    def _load_limiters(self, output):
        for i, limiter_name in enumerate(self._selected):
            limiter = self._get_limiter_by_name(limiter_name)
            if limiter is None:
                raise ValueError(f"Could not find {limiter_name}")

            logger.info(f"Selected {limiter.name}")
            vtk_poly = self._create_vtk_array(limiter)
            output.SetBlock(i, vtk_poly)

    def _get_limiter_by_name(self, limiter_name):
        for limiter in self._selectable:
            if limiter_name == limiter.name:
                return limiter
        return None

    def _create_vtk_array(self, limiter):
        if limiter.unit.closed == 0:
            is_closed = False
        else:
            is_closed = True

        r = limiter.unit.outline.r
        z = limiter.unit.outline.z
        assert len(r) == len(z), "r and z must have the same length."

        # Fill points
        vtk_points = vtkPoints()
        for i in range(len(r)):
            vtk_points.InsertNextPoint(r[i], 0.0, z[i])

        # Fill lines
        vtk_lines = vtkCellArray()
        for i in range(len(r) - 1):
            line = vtkLine()
            line.GetPointIds().SetId(0, i)
            line.GetPointIds().SetId(1, i + 1)
            vtk_lines.InsertNextCell(line)

        # Close loop if the unit is closed
        if is_closed:
            line = vtkLine()
            line.GetPointIds().SetId(0, len(r) - 1)
            line.GetPointIds().SetId(1, 0)
            vtk_lines.InsertNextCell(line)

        vtk_poly = vtkPolyData()
        vtk_poly.SetPoints(vtk_points)
        vtk_poly.SetLines(vtk_lines)
        return vtk_poly
