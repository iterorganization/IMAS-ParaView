"""Plugin to visualize coils and loops from the pf_active, pf_passive and
coils_non_axisymmetric IDSs"""

import logging

import numpy as np
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkMultiBlockDataSet,
    vtkPolyData,
)

from imas_paraview.paraview_support.servermanager_tools import intvector, propertygroup
from imas_paraview.plugins.base_class import GGDVTKPluginBase
from imas_paraview.util import pol_to_cart

logger = logging.getLogger("imas_paraview")

SUPPORTED_IDS_NAMES = ["coils_non_axisymmetric"]

# TODO: add tests
# TODO: add docs


@smproxy.source(label="Non-Axisymmetric Coils Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class CoilsNonAxisymmetricReader(GGDVTKPluginBase):
    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", SUPPORTED_IDS_NAMES)
        self.selectable_map = {}
        self.resolution = 10

    @intvector(label="Resolution", name="resolution", default_values=10)
    def P99_SetResolution(self, val):
        """Sets the number of points for 'arcs_of_circle' geometry type, if it is
        available in the loaded IDS."""
        self._update_property("resolution", val)

    @propertygroup("Non-Axisymmetric Coils Reader Settings", ["resolution"])
    def PG3_CoilsNonAxisymmetricReaderGroup(self):
        """Dummy function to define a PropertyGroup."""

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        if len(self._selected) > 0:
            output = vtkMultiBlockDataSet.GetData(outInfo)
            self._convert_to_vtk(output)

        return 1

    def setup_ids(self):
        assert self._ids is not None, "IDS cannot be empty during setup."

        self.selectable_map = {}
        self._load_coils(self._ids.coil)

    def _load_coils(self, ids_quantity):
        for i, coil in enumerate(ids_quantity):
            coil_name = coil.name
            if not coil_name:
                coil_name = f"coil {i}"
                logger.warning(
                    "Non-axisymmetric coil without name found. Using %s", coil_name
                )

            if len(coil.conductor) == 0:
                logger.warning("%s has no conductors, skipping it.", coil_name)
                continue

            has_elements = False
            for conductor in coil.conductor:
                if len(conductor.elements.types) > 0:
                    has_elements = True
                    break

            if not has_elements:
                logger.warning(
                    "%s has no elements in any conductor, skipping it.", coil_name
                )
                continue

            self.selectable_map[str(coil_name)] = coil
        self._selectable = list(self.selectable_map.keys())

    def _convert_to_vtk(self, output: vtkMultiBlockDataSet):
        block_index = 0
        for coil_name in self._selected:
            coil = self.selectable_map[coil_name]
            for conductor_idx, conductor in enumerate(coil.conductor):
                elements = conductor.elements
                conductor_points = []
                lines = vtkCellArray()

                if len(elements.types) == 0:
                    logger.warning(
                        "elements of conductor %d does not have types, skipping",
                        conductor_idx,
                    )
                    continue

                point_id = 0
                num_elements_added = 0

                for elem_idx, elem_type in enumerate(elements.types):
                    elem_points = np.array([])

                    # TODO: implement cross section
                    if elem_type == 1:
                        elem_points = self._create_line_segment(elements, elem_idx)
                    elif elem_type == 2:
                        elem_points = self._create_arc_of_circle(elements, elem_idx)
                    else:
                        logger.warning(
                            "%s conductor %d element %d has unsupported element "
                            "type %d, skipping",
                            coil_name,
                            conductor_idx,
                            elem_idx,
                            elem_type,
                        )

                    if len(elem_points) != 0:
                        conductor_points.append(elem_points)
                        num_pts = len(elem_points)
                        lines.InsertNextCell(num_pts)
                        for i in range(num_pts):
                            lines.InsertCellPoint(point_id + i)
                        point_id += num_pts
                        num_elements_added += 1

                if not conductor_points:
                    continue

                vtk_pts = vtkPoints()
                vtk_pts.SetData(numpy_to_vtk(np.vstack(conductor_points)))
                poly = vtkPolyData()
                poly.SetPoints(vtk_pts)
                poly.SetLines(lines)
                output.SetBlock(block_index, poly)
                block_index += 1
                logger.info(
                    "Loaded conductor %d from coil %r with %d element(s)",
                    conductor_idx,
                    coil_name,
                    num_elements_added,
                )
            logger.info(
                "Loaded non-axisymmetric coil %s with %d conductor(s)",
                coil_name,
                len(coil.conductor),
            )

    def _create_line_segment(self, elements, idx):
        r_start = elements.start_points.r[idx]
        phi_start = elements.start_points.phi[idx]
        z_start = elements.start_points.z[idx]

        r_end = elements.end_points.r[idx]
        phi_end = elements.end_points.phi[idx]
        z_end = elements.end_points.z[idx]

        x_start, y_start = pol_to_cart(r_start, phi_start)
        x_end, y_end = pol_to_cart(r_end, phi_end)

        return np.array([[x_start, y_start, z_start], [x_end, y_end, z_end]])

    def _pol_to_cart3d(self, point, idx):
        x, y = pol_to_cart(point.r[idx], point.phi[idx])
        return np.array([x, y, point.z[idx]])

    def _create_arc_of_circle(self, elements, idx):
        p_start = self._pol_to_cart3d(elements.start_points, idx)
        p_intermediate = self._pol_to_cart3d(elements.intermediate_points, idx)
        p_end = self._pol_to_cart3d(elements.end_points, idx)
        p_centre = self._pol_to_cart3d(elements.centres, idx)

        # Vectors from center of circle to start/end points
        v_start = p_start - p_centre
        v_end = p_end - p_centre
        radius = np.linalg.norm(v_start)
        v_start /= radius

        binormal = np.cross(v_start, p_intermediate - p_centre)
        binormal /= np.linalg.norm(binormal)

        # Tangent at start point
        tangent = np.cross(-v_start, binormal)

        # Sweep circle arc from start to end point
        angle = np.arctan2(np.dot(v_end, tangent), np.dot(v_end, v_start))
        if angle < 0:
            angle += 2 * np.pi

        t = np.linspace(0, angle, self.resolution)[:, np.newaxis]
        return p_centre + radius * (np.cos(t) * v_start + np.sin(t) * tangent)
