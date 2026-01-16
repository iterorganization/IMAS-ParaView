"""Plugin to visualize coils and loops from the pf_active, pf_passive and
coils_non_axisymmetric IDSs"""

import logging

import numpy as np
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet
from vtkmodules.vtkFiltersCore import vtkAppendPolyData, vtkTubeFilter

from imas_paraview.paraview_support.servermanager_tools import intvector, propertygroup
from imas_paraview.plugins.base_class import GGDVTKPluginBase
from imas_paraview.util import points_to_vtkpoly, pol_to_cart

logger = logging.getLogger("imas_paraview")

SUPPORTED_IDS_NAMES = ["coils_non_axisymmetric"]

# TODO: add docs / docstrings


@smproxy.source(label="Non-Axisymmetric Coils Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class CoilsNonAxisymmetricReader(GGDVTKPluginBase):
    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", SUPPORTED_IDS_NAMES)
        self.selectable_map = {}
        self.resolution = 10
        self.cs_resolution = 10

    @intvector(label="Resolution", name="resolution", default_values=10)
    def P98_SetResolution(self, val):
        """Sets the number of points for interpolating the 'arcs_of_circle' and
        'circle' geometrical element type, if they are available in the loaded IDS."""
        self._update_property("resolution", val)

    @intvector(
        label="Cross-sectional Resolution", name="cs_resolution", default_values=10
    )
    def P99_SetCrossSectionalResolution(self, val):
        """Sets the number of points for interpolating the 'annulus' cross-section, if
        it is available in the loaded IDS."""
        self._update_property("cs_resolution", val)

    @propertygroup(
        "Non-Axisymmetric Coils Reader Settings", ["resolution", "cs_resolution"]
    )
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
            # Coil names are not unique in some machine description IDSs
            coil_name = f"{coil.name} / {coil.identifier}"
            if not coil_name:
                coil_name = f"coil {i}"
                logger.warning(
                    "Non-axisymmetric coil without name found. Using %r", coil_name
                )

            if len(coil.conductor) == 0:
                logger.warning("%r has no conductors, skipping it.", coil_name)
                continue

            has_elements = False
            for conductor in coil.conductor:
                if len(conductor.elements.types) > 0:
                    has_elements = True
                    break

            if not has_elements:
                logger.warning(
                    "%r has no elements in any conductor, skipping it.", coil_name
                )
                continue

            self.selectable_map[str(coil_name)] = coil
        self._selectable = list(self.selectable_map.keys())

    def _convert_to_vtk(self, output: vtkMultiBlockDataSet):
        block_index = 0
        for coil_name in self._selected:
            coil = self.selectable_map[coil_name]
            logger.info("Loading non-axisymmetric coil %r...", coil_name)

            for conductor in coil.conductor:
                if len(conductor.elements.types) == 0:
                    logger.warning(
                        "The geometrical elements of the conductor do not have a type, "
                        "this conductor will be skipped",
                    )
                    continue

                conductor_vtk_geometry = self.create_conductor_geometry(conductor)
                output.SetBlock(block_index, conductor_vtk_geometry)
                block_index += 1

            logger.info(
                "Loaded non-axisymmetric coil %r with %d conductor(s)",
                coil_name,
                len(coil.conductor),
            )

    def create_conductor_geometry(self, conductor):
        elements = conductor.elements
        conductor_block = vtkAppendPolyData()

        for elem_idx, elem_type in enumerate(elements.types):
            if elem_type == 1:  # Line segment
                elem_points = self._create_line_segment(elements, elem_idx)
            elif elem_type == 2:  # Arc of a circle
                elem_points = self._create_circular_geometry(
                    elements, elem_idx, is_full_circle=False
                )
            elif elem_type == 3:  # Full circle
                elem_points = self._create_circular_geometry(
                    elements, elem_idx, is_full_circle=True
                )
            else:
                logger.warning(
                    "Conductor element %d has unsupported element type %d, skipping",
                    elem_idx,
                    elem_type,
                )
                continue

            vtk_conductor = points_to_vtkpoly(elem_points)

            if len(conductor.cross_section) > 0:
                vtk_conductor = self._add_cross_section(
                    vtk_conductor, conductor, elem_idx
                )

            conductor_block.AddInputData(vtk_conductor)
        conductor_block.Update()
        return conductor_block.GetOutput()

    def _pol_to_cart3d(self, point, idx):
        x, y = pol_to_cart(point.r[idx], point.phi[idx])
        return np.array([x, y, point.z[idx]])

    def _create_line_segment(self, elements, idx):
        p_start = self._pol_to_cart3d(elements.start_points, idx)
        p_end = self._pol_to_cart3d(elements.end_points, idx)
        return np.array([p_start, p_end])

    def _create_circular_geometry(self, elements, idx, is_full_circle):
        p_start = self._pol_to_cart3d(elements.start_points, idx)
        p_intermediate = self._pol_to_cart3d(elements.intermediate_points, idx)
        p_centre = self._pol_to_cart3d(elements.centres, idx)

        # Vector from center of circle to start point
        v_start = p_start - p_centre
        radius = np.linalg.norm(v_start)
        v_start /= radius

        binormal = np.cross(v_start, p_intermediate - p_centre)
        binormal /= np.linalg.norm(binormal)

        # Tangent at start point
        tangent = np.cross(-v_start, binormal)

        if is_full_circle:
            max_angle = 2 * np.pi
            t = np.linspace(0, max_angle, self.resolution + 1, endpoint=True)[:, None]
        else:
            # Sweep circle arc from start to end point
            p_end = self._pol_to_cart3d(elements.end_points, idx)
            v_end = p_end - p_centre
            max_angle = np.arctan2(np.dot(v_end, tangent), np.dot(v_end, v_start))
            if max_angle < 0:
                max_angle += 2 * np.pi
            t = np.linspace(0, max_angle, self.resolution)[:, np.newaxis]
        return p_centre + radius * (np.cos(t) * v_start + np.sin(t) * tangent)

    def _add_cross_section(self, input_poly_line, conductor, elem_idx):
        cross_section_idx = 0 if len(conductor.cross_section) == 1 else elem_idx
        cross_section = conductor.cross_section[cross_section_idx]

        if cross_section.geometry_type.index == 5:  # Annulus
            return self._add_annulus_cross_section(input_poly_line, cross_section)
        else:
            logger.warning(
                "Cross-section %d is not supported, it will be represented as a line ",
                cross_section.geometry_type.index,
            )
            return input_poly_line

    def _add_annulus_cross_section(self, poly_line, cross_section):
        outer_tube = vtkTubeFilter()
        outer_tube.SetInputData(poly_line)
        outer_tube.SetRadius(cross_section.width / 2.0)
        outer_tube.SetNumberOfSides(self.cs_resolution)
        outer_tube.Update()

        inner_tube = vtkTubeFilter()
        inner_tube.SetInputData(poly_line)
        inner_tube.SetRadius(cross_section.radius_inner)
        inner_tube.SetNumberOfSides(self.cs_resolution)
        inner_tube.Update()

        append = vtkAppendPolyData()
        append.AddInputData(outer_tube.GetOutput())
        append.AddInputData(inner_tube.GetOutput())
        append.Update()

        return append.GetOutput()
