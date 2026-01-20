"""Plugin to visualize coil conductors geometries from the coils_non_axisymmetric IDS"""

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


@smproxy.source(label="Non-Axisymmetric Coils Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class CoilsNonAxisymmetricReader(GGDVTKPluginBase):
    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", SUPPORTED_IDS_NAMES)
        self.selectable_map = {}
        self.resolution = 10
        """Resolution for arcs_of_circle and full_circle geometries"""
        self.cs_resolution = 10
        """Cross-sectional resolution"""

    @intvector(label="Resolution", name="resolution", default_values=10)
    def P98_SetResolution(self, val):
        """Sets the number of points for interpolating the 'arcs of circle' and 'full
        circle' geometrical element type, if they are available in the loaded IDS."""
        self._update_property("resolution", val)

    @intvector(
        label="Cross-sectional Resolution", name="cs_resolution", default_values=10
    )
    def P99_SetCrossSectionalResolution(self, val):
        """Sets the number of points for interpolating a circular cross-section, if
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
        """Select which coils to show in the array domain selector. Skips coils without
        any conductor elements."""
        assert self._ids is not None, "IDS cannot be empty during setup."

        self.selectable_map = {}
        for i, coil in enumerate(self._ids.coil):
            coil_name = coil.name
            if not coil_name:
                coil_name = f"coil {i}"
                logger.warning(
                    "Non-axisymmetric coil without name found. Using '%s'", coil_name
                )
            # Coil names are not unique for DD3.x
            if hasattr(coil, "identifier"):
                coil_name = f"{coil_name} / {coil.identifier}"

            if len(coil.conductor) == 0:
                logger.warning("'%s' has no conductors, skipping it.", coil_name)
                continue

            has_elements = False
            for conductor in coil.conductor:
                if len(conductor.elements.types) > 0:
                    has_elements = True
                    break

            if not has_elements:
                logger.warning(
                    "'%s' has no elements in any conductor, skipping it.", coil_name
                )
                continue

            self.selectable_map[str(coil_name)] = coil
        self._selectable = list(self.selectable_map.keys())

    def _convert_to_vtk(self, output: vtkMultiBlockDataSet):
        """Convert selected coil conductors to a VTK multi-block dataset.

        Args:
            output: vtkMultiBlockDataSet where each conductor will be added as a block.
        """
        block_index = 0
        for coil_name in self._selected:
            coil = self.selectable_map[coil_name]
            logger.info("Loading non-axisymmetric coil '%s'...", coil_name)

            for cond_idx, conductor in enumerate(coil.conductor):
                if len(conductor.elements.types) == 0:
                    logger.warning(
                        "The geometrical elements of the conductor do not have a type, "
                        "this conductor will be skipped",
                    )
                    continue

                conductor_vtk_geometry = self.create_conductor_geometry(conductor)
                if conductor_vtk_geometry is None:
                    logger.warning(
                        "Conductor %d does not have a valid geometry, skipping",
                        cond_idx,
                        coil_name,
                    )
                    continue
                output.SetBlock(block_index, conductor_vtk_geometry)
                block_index += 1

            logger.info(
                "Loaded non-axisymmetric coil '%s' with %d conductor(s)",
                coil_name,
                block_index,
            )

    def create_conductor_geometry(self, conductor):
        """Create VTK geometry representation of the conductor.

        Args:
            conductor: Conductor IDS object containing geometrical elements and
                optionally cross-sections.

        Returns:
            vtkAppendPolyData containing the complete conductor geometry, or None if
            there are no valid elements.
        """
        elements = conductor.elements
        conductor_block = vtkAppendPolyData()

        has_input = False
        has_cross_section = False

        if len(conductor.cross_section) == 1 or len(conductor.cross_section) == len(
            elements.types
        ):
            has_cross_section = True
        elif len(conductor.cross_section) == 0:
            logger.warning(
                "Conductor does not have a cross-section, only the centreline will "
                "be shown.",
            )
        else:
            logger.warning(
                "Conductor must have either 1 universal cross-section or a separate "
                "cross-section for each element. Only the centreline will be shown."
            )

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

            if elem_points is None:
                logger.warning(
                    "element %d has an invalid geometrical element, skipping", elem_idx
                )
                continue
            vtk_conductor = points_to_vtkpoly(elem_points)

            if has_cross_section:
                vtk_conductor = self._add_cross_section(
                    vtk_conductor, conductor, elem_idx
                )

            conductor_block.AddInputData(vtk_conductor)
            has_input = True

        if not has_input:
            logger.warning(
                "Conductor does not have any valid geometrical elements, skipping"
            )
            return None
        conductor_block.Update()
        return conductor_block.GetOutput()

    def _pol_to_cart3d(self, point, idx):
        """Convert cylindrical coordinates of a point to cartesian coordinates."""
        x, y = pol_to_cart(point.r[idx], point.phi[idx])
        return np.array([x, y, point.z[idx]])

    def _create_line_segment(self, elements, idx):
        """Create line segment points for a conductor element.

        Args:
            elements: Element IDS quantity with start and end points.
            idx: Index of the element.

        Returns:
            Array of points for the line segment.
        """
        p_start = self._pol_to_cart3d(elements.start_points, idx)
        p_end = self._pol_to_cart3d(elements.end_points, idx)
        return np.array([p_start, p_end])

    def _create_circular_geometry(self, elements, idx, is_full_circle):
        """Create points for circular conductor elements.

        Args:
            elements: Element container with start, intermediate, and center points.
            idx: Index of the element.
            is_full_circle: True if full circle, False if arc of circle.

        Returns:
            Array of points representing the circular geometry, or None if
            the circular element is invalid.
        """
        p_start = self._pol_to_cart3d(elements.start_points, idx)
        p_intermediate = self._pol_to_cart3d(elements.intermediate_points, idx)
        p_centre = self._pol_to_cart3d(elements.centres, idx)

        # Vector from center of circle to start point
        v_start = p_start - p_centre
        radius = np.linalg.norm(v_start)
        v_start /= radius

        if not np.isclose(radius, np.linalg.norm(p_intermediate - p_centre)):
            logger.warning(
                "Start and intermediate point of element %d are not equidistant from "
                "the centre point",
                idx,
            )
            return None

        binormal = np.cross(v_start, p_intermediate - p_centre)
        binormal /= np.linalg.norm(binormal)

        # Tangent at start point
        tangent = np.cross(-v_start, binormal)

        if is_full_circle:
            max_angle = 2 * np.pi
            resolution = self.resolution + 1  # include end point
        else:
            # Sweep circle arc from start to end point
            p_end = self._pol_to_cart3d(elements.end_points, idx)
            v_end = p_end - p_centre

            if not np.isclose(radius, np.linalg.norm(v_end)):
                logger.warning(
                    "Start and end point of element %d are not equidistant from "
                    "the centre point",
                    idx,
                )
                return None

            if not self._are_points_coplanar(p_start, p_intermediate, p_end, p_centre):
                logger.warning("Element %d points are not coplanar", idx)
                return None

            max_angle = np.arctan2(np.dot(v_end, tangent), np.dot(v_end, v_start))
            if max_angle < 0:
                max_angle += 2 * np.pi
            resolution = self.resolution
        t = np.linspace(0, max_angle, resolution)[:, None]
        return p_centre + radius * (np.cos(t) * v_start + np.sin(t) * tangent)

    def _are_points_coplanar(self, p0, p1, p2, p3):
        """Checks if 4 points are coplanar."""
        n = np.cross(p1 - p0, p2 - p0)
        if np.linalg.norm(n) == 0.0:
            return False
        dist = np.dot(p3 - p0, n)
        return np.isclose(dist, 0.0)

    def _add_cross_section(self, input_poly_line, conductor, elem_idx):
        """Add cross-sectional representation to a conductor polyline.

        Args:
            input_poly_line: VTK polyline of the conductor element.
            conductor: Conductor object containing cross-sections.
            elem_idx: Index of the conductor element.

        Returns:
            VTK polydata with cross-section geometry.
        """
        cross_section_idx = 0 if len(conductor.cross_section) == 1 else elem_idx
        cross_section = conductor.cross_section[cross_section_idx]

        if cross_section.geometry_type.index == 2:  # Circle
            return self._add_circular_cross_section(
                input_poly_line, cross_section, "circle"
            )
        elif cross_section.geometry_type.index == 5:  # Annulus
            return self._add_circular_cross_section(
                input_poly_line, cross_section, "annulus"
            )
        else:
            logger.warning(
                "Cross-section %d is not supported, it will be represented as a line ",
                cross_section.geometry_type.index,
            )
            return input_poly_line

    def _add_circular_cross_section(self, poly_line, cross_section, cs_type):
        """Create a circular cross-section around a polyline using VTK tube filters.

        Args:
            poly_line: VTK polyline representing conductor element.
            cross_section: Cross-section object containing inner radius and width.
            cs_type: Type of cross section to add, either "circle" or "annulus"

        Returns:
            VTK polydata representing the annulus geometry.
        """
        assert cs_type in ["circle", "annulus"], "invalid cross-section type"
        append = vtkAppendPolyData()

        outer_tube = vtkTubeFilter()
        outer_tube.SetInputData(poly_line)
        outer_tube.SetRadius(cross_section.width / 2.0)
        outer_tube.SetNumberOfSides(self.cs_resolution)
        outer_tube.Update()
        append.AddInputData(outer_tube.GetOutput())

        if cs_type == "annulus":
            inner_tube = vtkTubeFilter()
            inner_tube.SetInputData(poly_line)
            inner_tube.SetRadius(cross_section.radius_inner)
            inner_tube.SetNumberOfSides(self.cs_resolution)
            inner_tube.Update()
            append.AddInputData(inner_tube.GetOutput())

        append.Update()
        return append.GetOutput()
