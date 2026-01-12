"""Plugin to view coils and their conductors"""

import logging

import numpy as np
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkMultiBlockDataSet,
    vtkPolyData,
)

from imas_paraview.plugins.base_class import GGDVTKPluginBase

logger = logging.getLogger("imas_paraview")

SUPPORTED_IDS_NAMES = ["pf_active", "pf_passive", "coils_non_axisymmetric"]


@smproxy.source(label="Coils Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class CoilsReader(GGDVTKPluginBase):
    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", SUPPORTED_IDS_NAMES)
        self.selectable_map = {}

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        if len(self._selected) > 0:
            output = vtkMultiBlockDataSet.GetData(outInfo)
            self._create_coil_vtk_objects(output)
        return 1

    def setup_ids(self):
        assert self._ids is not None, "IDS cannot be empty during setup."

        self.selectable_map = {}

        if self._ids.metadata.name == "pf_active":
            self._load_pf_active_coils()
        elif self._ids.metadata.name == "pf_passive":
            raise NotImplementedError("not implemented")
        elif self._ids.metadata.name == "coils_non_axisymmetric":
            raise NotImplementedError("not implemented")

    def _load_pf_active_coils(self):
        for i, coil in enumerate(self._ids.coil):
            coil_name = coil.name
            if not coil_name:
                coil_name = f"coil {i}"
                logger.warning(f"PF active coil without name found. Using {coil_name}")

            if len(coil.element) == 0:
                logger.warning(
                    f"{coil_name} has no elements, skipping selectable entries"
                )
                continue

            self.selectable_map[str(coil_name)] = coil

        self._selectable = list(self.selectable_map.keys())

    def _create_rectangular_coil(self, element):
        rect = element.geometry.rectangle
        r, z, width, height = rect.r, rect.z, rect.width, rect.height

        r0 = r - width / 2.0
        r1 = r + width / 2.0
        z0 = z - height / 2.0
        z1 = z + height / 2.0

        pts = vtkPoints()
        cells = vtkCellArray()

        rectangle_pts = [
            (r0, z0),
            (r1, z0),
            (r1, z1),
            (r0, z1),
            (r0, z0),  # close loop
        ]

        cells.InsertNextCell(len(rectangle_pts))
        for r_i, z_i in rectangle_pts:
            pts.InsertNextPoint(r_i, 0.0, z_i)
            cells.InsertCellPoint(pts.GetNumberOfPoints() - 1)

        poly = vtkPolyData()
        poly.SetPoints(pts)
        poly.SetLines(cells)
        return poly

    def _create_annulus_coil(self, element, resolution=10):
        ann = element.geometry.annulus
        r0, z0, r_in, r_out = ann.r, ann.z, ann.radius_inner, ann.radius_outer

        outer_ids = []
        inner_ids = []
        pts = vtkPoints()

        for i in range(resolution):
            theta = 2.0 * np.pi * i / resolution
            c = np.cos(theta)
            s = np.sin(theta)
            outer_ids.append(pts.InsertNextPoint(r0 + r_out * c, 0.0, z0 + r_out * s))
            inner_ids.append(pts.InsertNextPoint(r0 + r_in * c, 0.0, z0 + r_in * s))

        def fill_cells(cells, point_ids):
            cells.InsertNextCell(resolution + 1)
            for point_id in point_ids:
                cells.InsertCellPoint(point_id)
            cells.InsertCellPoint(point_ids[0])

        cells = vtkCellArray()
        fill_cells(cells, outer_ids)
        fill_cells(cells, inner_ids)

        poly = vtkPolyData()
        poly.SetPoints(pts)
        poly.SetLines(cells)
        return poly

    def _create_coil_vtk_objects(self, output: vtkMultiBlockDataSet):
        block_index = 0

        for coil_name in self._selected:
            coil = self.selectable_map[coil_name]
            for element in coil.element:
                geom_type = element.geometry.geometry_type

                if geom_type == 2:  # rectangle
                    vtk_coil = self._create_rectangular_coil(element)
                elif geom_type == 5:  # annulus
                    vtk_coil = self._create_annulus_coil(element)
                else:
                    logger.warning(
                        f"{coil_name} has unsupported geometry type: {geom_type}, skipping"
                    )
                    continue

                output.SetBlock(block_index, vtk_coil)
                output.GetMetaData(block_index).Set(output.NAME(), str(coil_name))
                block_index += 1

            logger.info(
                f"Loaded PF coil {coil_name} with {len(coil.element)} element(s)"
            )
