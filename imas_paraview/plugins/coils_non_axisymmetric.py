"""Plugin to visualize coil conductor geometries from the coils_non_axisymmetric
and tf IDSs"""

import logging

import imas
import numpy as np
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkMultiBlockDataSet,
    vtkPolyData,
)
from vtkmodules.vtkFiltersCore import vtkAppendPolyData, vtkTubeFilter

from imas_paraview.paraview_support.servermanager_tools import intvector, propertygroup
from imas_paraview.plugins.base_class import GGDVTKPluginBase
from imas_paraview.util import points_to_vtkpoly, pol_to_cart

logger = logging.getLogger("imas_paraview")

SUPPORTED_IDS_NAMES = ["coils_non_axisymmetric", "tf"]


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

            for conductor in coil.conductor:
                if len(conductor.elements.types) == 0:
                    logger.warning(
                        "The geometrical elements of the conductor do not have a type, "
                        "this conductor will be skipped",
                    )
                    continue

                conductor_vtk_geometry = self.create_conductor_geometry(conductor)
                if conductor_vtk_geometry is None:
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
        else:
            logger.warning(
                "Conductor '%s' must have either a single cross-section or a separate "
                "cross-section for each element. Only the centreline will be shown.",
                imas.util.get_full_path(conductor),
            )

        # (normal, binormal, tangent) frame at the end of the last successfully
        # processed element. Used to parallel-transport the cross-section
        # orientation across line segments that do not provide an intermediate
        # point (which is only needed to fix their orientation).
        prev_frame = None

        for elem_idx, elem_type in enumerate(elements.types):
            is_closed = elem_type == 3
            line_frame = None
            circ_frame = None
            if elem_type == 1:  # Line segment
                elem_points = self._create_line_segment(elements, elem_idx)
            elif elem_type in (2, 3):  # Arc of a circle or full circle
                circular_result = self._circular_points_and_frames(
                    elements, elem_idx, is_full_circle=is_closed
                )
                if circular_result is None:
                    elem_points = None
                else:
                    elem_points, *circ_frame = circular_result
            else:
                logger.warning(
                    "Element %d of '%s' has unsupported element type %d, skipping",
                    elem_idx,
                    imas.util.get_full_path(elements),
                    elem_type,
                )
                continue

            if elem_points is None:
                continue

            if elem_type == 1:
                line_frame = self._create_line_frame(elements, elem_idx, prev_frame)
                frame = None if line_frame is None else line_frame[:2]
            else:
                frame = None if circ_frame is None else circ_frame[:2]

            vtk_conductor = points_to_vtkpoly(elem_points)

            if has_cross_section:
                vtk_conductor = self._add_cross_section(
                    vtk_conductor, conductor, elem_idx, elem_points, frame, is_closed
                )

            conductor_block.AddInputData(vtk_conductor)
            has_input = True

            if line_frame is not None:
                prev_frame = line_frame
            elif circ_frame is not None:
                prev_frame = (circ_frame[0][-1], circ_frame[1][-1], circ_frame[2][-1])

        if not has_input:
            logger.warning(
                "Conductor '%s' does not have any valid geometrical elements, skipping",
                imas.util.get_full_path(conductor),
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

    def _circular_arc_params(self, elements, idx, is_full_circle):
        """Compute the shared geometric parameters of an arc/circle conductor element.

        Args:
            elements: Element container with start, intermediate, and center points.
            idx: Index of the element.
            is_full_circle: True if full circle, False if arc of circle.

        Returns:
            Tuple of (p_centre, radius, v_start, tangent, binormal, t), where v_start
            and tangent are the unit vectors pointing from the centre to the start
            point, and along the direction of travel at the start point,
            respectively, binormal is the (constant) unit vector normal to the plane
            of the circle, and t is the array of angles (in radians, relative to
            v_start) sampled along the arc/circle. Returns None if the points
            defining the element are invalid.
        """
        p_start = self._pol_to_cart3d(elements.start_points, idx)
        p_intermediate = self._pol_to_cart3d(elements.intermediate_points, idx)
        p_centre = self._pol_to_cart3d(elements.centres, idx)
        p_end = (  # End point is only defined for arc of circle
            None if is_full_circle else self._pol_to_cart3d(elements.end_points, idx)
        )

        if not self._are_circular_points_valid(
            p_start, p_intermediate, p_end, p_centre, is_full_circle
        ):
            return None

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
            resolution = self.resolution + 1  # include end point
        else:
            # Sweep circle arc from start to end point
            v_end = p_end - p_centre
            max_angle = np.arctan2(np.dot(v_end, tangent), np.dot(v_end, v_start))
            if max_angle < 0:
                max_angle += 2 * np.pi
            resolution = self.resolution
        t = np.linspace(0, max_angle, resolution)[:, None]
        return p_centre, radius, v_start, tangent, binormal, t

    def _circular_points_and_frames(self, elements, idx, is_full_circle):
        """Compute path points and (normal, binormal, tangent) orientation frames
        for a circular conductor element, computing the shared
        :meth:`_circular_arc_params` only once.

        Per the Data Dictionary convention, the binormal is constant along the arc
        (the rotation axis of the circle), while the normal and tangent rotate
        together with the point position (normal = centre - point on curve).

        Args:
            elements: Element container with start, intermediate, and center points.
            idx: Index of the element.
            is_full_circle: True if full circle, False if arc of circle.

        Returns:
            Tuple of (points, normals, binormals, tangents), where points is an
            array of 3D points representing the circular geometry, and normals,
            binormals, tangents are arrays with one unit 3-vector per path point.
            Returns None if the circular element is invalid.
        """
        params = self._circular_arc_params(elements, idx, is_full_circle)
        if params is None:
            return None
        p_centre, radius, v_start, tangent_start, binormal, t = params
        radial_dir = np.cos(t) * v_start + np.sin(t) * tangent_start
        points = p_centre + radius * radial_dir
        normals = -radial_dir
        tangents = -np.sin(t) * v_start + np.cos(t) * tangent_start
        binormals = np.broadcast_to(binormal, normals.shape)
        return points, normals, binormals, tangents

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
        result = self._circular_points_and_frames(elements, idx, is_full_circle)
        return None if result is None else result[0]

    def _create_circular_frames(self, elements, idx, is_full_circle):
        """Compute the (normal, binormal, tangent) cross-section orientation frame
        at every point returned by :meth:`_create_circular_geometry` for the same
        element.

        Args:
            elements: Element container with start, intermediate, and center points.
            idx: Index of the element.
            is_full_circle: True if full circle, False if arc of circle.

        Returns:
            Tuple of (normals, binormals, tangents) arrays with one unit 3-vector
            per path point, or None if the element is invalid.
        """
        result = self._circular_points_and_frames(elements, idx, is_full_circle)
        return None if result is None else result[1:]

    def _create_line_frame(self, elements, idx, prev_frame):
        """Compute the (normal, binormal, tangent) cross-section orientation frame
        of a line segment conductor element.

        Per the Data Dictionary convention: tangent = end point - start point;
        normal = intermediate point - start point; binormal = tangent x normal.
        The intermediate point is only used to fix the orientation of the
        cross-section and may not always be filled (e.g. some machine description
        datasets omit it for line elements). In that case, the frame is instead
        parallel-transported from the previous conductor element (or picked
        arbitrarily if there is none), to keep the cross-section orientation
        continuous instead of discarding it entirely.

        Args:
            elements: Element container with start, end, and intermediate points.
            idx: Index of the element.
            prev_frame: (normal, binormal, tangent) frame at the end of the
                previously processed conductor element, or None if there is none.

        Returns:
            Tuple of (normal, binormal, tangent) unit 3-vectors, or None if the
            line segment itself is degenerate (zero length).
        """
        p_start = self._pol_to_cart3d(elements.start_points, idx)
        p_end = self._pol_to_cart3d(elements.end_points, idx)

        tangent = p_end - p_start
        tangent_norm = np.linalg.norm(tangent)
        if np.isclose(tangent_norm, 0):
            logger.warning(
                "Cannot determine cross-section orientation of line segment: "
                "start and end point coincide"
            )
            return None
        tangent = tangent / tangent_norm

        normal = binormal = None
        if idx < len(elements.intermediate_points.r):
            p_intermediate = self._pol_to_cart3d(elements.intermediate_points, idx)
            candidate_normal = p_intermediate - p_start
            candidate_normal_norm = np.linalg.norm(candidate_normal)
            if not np.isclose(candidate_normal_norm, 0):
                candidate_normal = candidate_normal / candidate_normal_norm
                candidate_binormal = np.cross(tangent, candidate_normal)
                candidate_binormal_norm = np.linalg.norm(candidate_binormal)
                if not np.isclose(candidate_binormal_norm, 0):
                    binormal = candidate_binormal / candidate_binormal_norm
                    # Re-orthogonalize the normal, in case the intermediate point
                    # was not exactly perpendicular to the tangent
                    normal = np.cross(binormal, tangent)

        if normal is None:
            # No (valid) intermediate point to derive the orientation from: keep
            # the cross-section orientation continuous with the previous element.
            if prev_frame is None:
                normal = self._arbitrary_perpendicular(tangent)
                binormal = np.cross(tangent, normal)
            else:
                prev_normal, prev_binormal, prev_tangent = prev_frame
                normal, binormal = self._parallel_transport_frame(
                    prev_normal, prev_binormal, prev_tangent, tangent
                )

        return normal, binormal, tangent

    def _arbitrary_perpendicular(self, tangent):
        """Return an arbitrary unit vector perpendicular to tangent.

        Args:
            tangent: Unit 3-vector to find a perpendicular vector for.

        Returns:
            Unit 3-vector perpendicular to tangent.
        """
        candidate = np.array([0.0, 0.0, 1.0])
        if abs(np.dot(candidate, tangent)) > 0.9:
            candidate = np.array([1.0, 0.0, 0.0])
        normal = candidate - np.dot(candidate, tangent) * tangent
        return normal / np.linalg.norm(normal)

    def _parallel_transport_frame(
        self, prev_normal, prev_binormal, prev_tangent, new_tangent
    ):
        """Parallel-transport a (normal, binormal) frame from prev_tangent to
        new_tangent, using the minimal rotation that maps one tangent onto the
        other (Rodrigues' rotation formula). This avoids introducing unnecessary
        twist in the cross-section when explicit orientation data is unavailable.

        Args:
            prev_normal: Normal unit vector of the frame to transport.
            prev_binormal: Binormal unit vector of the frame to transport.
            prev_tangent: Tangent unit vector the frame is currently aligned with.
            new_tangent: Tangent unit vector to align the frame with.

        Returns:
            Tuple of (normal, binormal) unit 3-vectors, transported to new_tangent.
        """
        cos_angle = np.clip(np.dot(prev_tangent, new_tangent), -1.0, 1.0)
        axis = np.cross(prev_tangent, new_tangent)
        axis_norm = np.linalg.norm(axis)
        if np.isclose(axis_norm, 0):
            if cos_angle > 0:
                return prev_normal, prev_binormal  # tangents are (nearly) identical
            axis = prev_normal  # tangents are opposite: any perpendicular axis works
            axis_norm = np.linalg.norm(axis)
        axis = axis / axis_norm
        angle = np.arccos(cos_angle)
        sin_a, cos_a = np.sin(angle), np.cos(angle)

        def rotate(v):
            return (
                v * cos_a
                + np.cross(axis, v) * sin_a
                + axis * np.dot(axis, v) * (1 - cos_a)
            )

        return rotate(prev_normal), rotate(prev_binormal)

    def _are_circular_points_valid(
        self, p_start, p_intermediate, p_end, p_centre, is_full_circle
    ):
        """Check geometric validity of a set of points of a circular element.

        Args:
            p_start: Start point of the circular element.
            p_intermediate: Intermediate point of the circular element.
            p_end: End point of the circular element.
            p_centre: Centre point of the circular element.
            is_full_circle: True if full circle, False if arc of circle.

        Returns:
            True if points make a valid circular element, or False if they don't
        """
        v_start = p_start - p_centre
        radius = np.linalg.norm(v_start)
        r_intermediate = np.linalg.norm(p_intermediate - p_centre)

        if np.isclose(radius, 0):
            logger.warning("Start point coincides with centre: zero radius")
            return False
        if not np.isclose(radius, r_intermediate):
            logger.warning(
                "Start and intermediate points are not equidistant from centre"
            )
            return False
        if not is_full_circle:
            r_end = np.linalg.norm(p_end - p_centre)
            if not np.isclose(radius, r_end):
                logger.warning("Start and end points are not equidistant from centre")
                return False

        cross = np.cross(v_start, p_intermediate - p_centre)
        if np.linalg.norm(cross) == 0:
            logger.warning(
                "The plane defined by start, intermediate, and centre point is "
                "degenerate"
            )
            return False

        if not is_full_circle:
            triple_product = np.dot(
                v_start, np.cross(p_intermediate - p_centre, p_end - p_centre)
            )
            if not np.isclose(triple_product, 0.0):
                logger.warning(
                    "Start, intermediate, end, and centre points are not coplanar"
                )
                return False
        return True

    def _add_cross_section(
        self, input_poly_line, conductor, elem_idx, elem_points, frame, is_closed
    ):
        """Add cross-sectional representation to a conductor polyline.

        Args:
            input_poly_line: VTK polyline of the conductor element.
            conductor: Conductor object containing cross-sections.
            elem_idx: Index of the conductor element.
            elem_points: Array of 3D points forming the conductor element centreline,
                as returned by e.g. :meth:`_create_line_segment`. Used to sweep the
                polygonal cross-section outline.
            frame: (normal, binormal) orientation frame of the cross-section along
                elem_points, as returned by e.g. :meth:`_create_line_frame`, or None.
            is_closed: True if elem_points forms a closed loop (full circle).

        Returns:
            VTK polydata with cross-section geometry.
        """
        cross_section_idx = 0 if len(conductor.cross_section) == 1 else elem_idx
        cross_section = conductor.cross_section[cross_section_idx]

        geometry_type = cross_section.geometry_type.index
        if geometry_type in (2, 5):  # circle (2) or annulus (5)
            return self._add_circular_cross_section(input_poly_line, cross_section)
        elif geometry_type == 1:  # polygon
            polygon = self._add_polygon_cross_section(
                elem_points, frame, cross_section, is_closed
            )
            return input_poly_line if polygon is None else polygon
        else:
            logger.warning(
                "Cross-section identifier %d of '%s' is not supported, it will be "
                "represented as a line instead",
                cross_section.geometry_type.index,
                imas.util.get_full_path(cross_section),
            )
            return input_poly_line

    def _add_polygon_cross_section(self, elem_points, frame, cross_section, is_closed):
        """Sweep a polygonal cross-section outline along a conductor element
        centreline.

        Args:
            elem_points: Array of 3D points forming the conductor element
                centreline.
            frame: Tuple of (normal, binormal) unit vectors describing the
                orientation of the cross-section along elem_points. Each of normal
                and binormal is either a single 3-vector (constant orientation, as
                used for line segments) or an array with one vector per point in
                elem_points (rotating orientation, as used for arcs/circles).
            cross_section: Cross-section IDS object containing the polygonal
                outline, given in local (normal, binormal) coordinates.
            is_closed: True if elem_points forms a closed loop (full circle), in
                which case the swept surface is not capped at the ends.

        Returns:
            VTK polydata representing the swept polygonal cross-section, or None if
            the outline or orientation frame is invalid.
        """
        if frame is None:
            return None
        normal, binormal = frame

        outline = cross_section.outline
        n_coords = np.asarray(outline.normal)
        b_coords = np.asarray(outline.binormal)
        n_ring = len(n_coords)
        if n_ring < 3 or len(b_coords) != n_ring:
            logger.warning(
                "Polygon cross-section of '%s' must have at least 3 outline points, "
                "with matching outline.normal/outline.binormal array lengths, it "
                "will be represented as a line instead",
                imas.util.get_full_path(cross_section),
            )
            return None

        elem_points = np.asarray(elem_points)
        n_path = len(elem_points)

        normals = np.broadcast_to(normal, (n_path, 3))
        binormals = np.broadcast_to(binormal, (n_path, 3))

        # Sweep the outline along the path: one ring of n_ring points per path point
        vertices = (
            elem_points[:, None, :]
            + n_coords[None, :, None] * normals[:, None, :]
            + b_coords[None, :, None] * binormals[:, None, :]
        )
        vertices = vertices.reshape(-1, 3)

        vtk_points = vtkPoints()
        vtk_points.SetData(numpy_to_vtk(vertices))

        polys = vtkCellArray()
        for i in range(n_path - 1):
            for j in range(n_ring):
                j_next = (j + 1) % n_ring
                polys.InsertNextCell(4)
                polys.InsertCellPoint(i * n_ring + j)
                polys.InsertCellPoint(i * n_ring + j_next)
                polys.InsertCellPoint((i + 1) * n_ring + j_next)
                polys.InsertCellPoint((i + 1) * n_ring + j)

        if not is_closed:
            # Cap both open ends of the swept surface with the outline polygon
            polys.InsertNextCell(n_ring)
            for j in range(n_ring):
                polys.InsertCellPoint(j)

            last_ring_start = (n_path - 1) * n_ring
            polys.InsertNextCell(n_ring)
            for j in reversed(range(n_ring)):
                polys.InsertCellPoint(last_ring_start + j)

        poly = vtkPolyData()
        poly.SetPoints(vtk_points)
        poly.SetPolys(polys)
        return poly

    def _add_circular_cross_section(self, poly_line, cross_section):
        """Create a circular cross-section around a polyline using VTK tube filters.

        Args:
            poly_line: VTK polyline representing conductor element.
            cross_section: Cross-section object containing inner radius and width.

        Returns:
            VTK polydata representing the geometry with a cross-section.
        """
        append = vtkAppendPolyData()

        outer_tube = vtkTubeFilter()
        outer_tube.SetInputData(poly_line)
        outer_tube.SetRadius(cross_section.width / 2.0)
        outer_tube.SetNumberOfSides(self.cs_resolution)
        outer_tube.Update()
        append.AddInputData(outer_tube.GetOutput())

        if cross_section.geometry_type.index == 5:  # annulus
            inner_tube = vtkTubeFilter()
            inner_tube.SetInputData(poly_line)
            inner_tube.SetRadius(cross_section.radius_inner)
            inner_tube.SetNumberOfSides(self.cs_resolution)
            inner_tube.Update()
            append.AddInputData(inner_tube.GetOutput())

        append.Update()
        return append.GetOutput()
