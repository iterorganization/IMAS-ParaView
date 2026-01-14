"""Plugin to visualize coils and loops from the pf_active, pf_passive and
coils_non_axisymmetric IDSs"""

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

SUPPORTED_IDS_NAMES = ["coils_non_axisymmetric"]

# TODO: add tests
# TODO: add docs


@smproxy.source(label="Non-Axisymmetric Coils Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class CoilsNonAxisymmetricReader(GGDVTKPluginBase):
    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", SUPPORTED_IDS_NAMES)
        self.selectable_map = {}

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

    def _load_coils(self, ids_quantityuantity):
        for i, coil in enumerate(ids_quantityuantity):
            coil_name = coil.name
            if not coil_name:
                coil_name = f"coil {i}"
                logger.warning(
                    f"Non-axisymmetric coil without name found. Using {coil_name}"
                )

            if len(coil.conductor) == 0:
                logger.warning(f"{coil_name} has no conductors, skipping it.")
                continue

            has_elements = False
            for conductor in coil.conductor:
                if len(conductor.elements.types) > 0:
                    has_elements = True
                    break

            if not has_elements:
                logger.warning(
                    f"{coil_name} has no elements in any conductor, skipping it."
                )
                continue

            self.selectable_map[str(coil_name)] = coil
        self._selectable = list(self.selectable_map.keys())

    def _convert_to_vtk(self, output: vtkMultiBlockDataSet):
        """Create VTK objects for non-axisymmetric coils with 3D conductor elements"""
        block_index = 0

        for coil_name in self._selected:
            coil = self.selectable_map[coil_name]

            for conductor_idx, conductor in enumerate(coil.conductor):
                elements = conductor.elements

                if len(elements.types) == 0:
                    continue

                pts = vtkPoints()
                cells = vtkCellArray()

                num_elements_added = 0
                for elem_idx in range(len(elements.types)):
                    elem_type = elements.types[elem_idx]

                    if elem_type == 1:
                        if self._add_line_segment(pts, cells, elements, elem_idx):
                            num_elements_added += 1
                    # elif elem_type == 2:  # Arc of circle
                    #     if self._add_arc_3d(pts, cells, elements, elem_idx):
                    #         num_elements_added += 1
                    else:
                        logger.warning(
                            f"{coil_name} conductor {conductor_idx} element {elem_idx} has "
                            f"unsupported element type {elem_type}: skipping"
                        )

                if pts.GetNumberOfPoints() > 0:
                    poly = vtkPolyData()
                    poly.SetPoints(pts)
                    poly.SetLines(cells)

                    output.SetBlock(block_index, poly)
                    conductor_name = f"{coil_name}_conductor_{conductor_idx}"
                    output.GetMetaData(block_index).Set(output.NAME(), conductor_name)
                    block_index += 1

                    logger.info(
                        f"Loaded {coil_name} conductor {conductor_idx} with {num_elements_added} element(s)"
                    )

            logger.info(
                f"Loaded non-axisymmetric coil {coil_name} with {len(coil.conductor)} conductor(s)"
            )

    def _add_line_segment(self, pts, cells, elements, idx, resolution=10):
        r_start = elements.start_points.r[idx]
        phi_start = elements.start_points.phi[idx]
        z_start = elements.start_points.z[idx]

        r_end = elements.end_points.r[idx]
        phi_end = elements.end_points.phi[idx]
        z_end = elements.end_points.z[idx]

        cells.InsertNextCell(resolution + 1)
        for i in range(resolution + 1):
            t = i / resolution
            r = r_start + t * (r_end - r_start)
            phi = phi_start + t * (phi_end - phi_start)
            z = z_start + t * (z_end - z_start)

            x = r * np.cos(phi)
            y = r * np.sin(phi)

            pts.InsertNextPoint(x, y, z)
            cells.InsertCellPoint(pts.GetNumberOfPoints() - 1)

        return True

    def _add_arc_3d(self, pts, cells, elements, idx, resolution=20):
        """Add an arc of circle in 3D (r, phi, z) coordinates to the points and cells"""
        r_start = elements.start_points.r[idx]
        phi_start = elements.start_points.phi[idx]
        z_start = elements.start_points.z[idx]

        r_mid = elements.intermediate_points.r[idx]
        phi_mid = elements.intermediate_points.phi[idx]
        z_mid = elements.intermediate_points.z[idx]

        r_end = elements.end_points.r[idx]
        phi_end = elements.end_points.phi[idx]
        z_end = elements.end_points.z[idx]

        r_center = elements.centres.r[idx]
        phi_center = elements.centres.phi[idx]
        z_center = elements.centres.z[idx]

        # Check for invalid data (placeholder value is -9e40)
        if (
            abs(r_start) > 1e30
            or abs(r_mid) > 1e30
            or abs(r_end) > 1e30
            or abs(r_center) > 1e30
        ):
            logger.warning(f"Arc element {idx} has invalid data, skipping")
            return False

        # Convert cylindrical to Cartesian
        x_start = r_start * np.cos(phi_start)
        y_start = r_start * np.sin(phi_start)

        x_mid = r_mid * np.cos(phi_mid)
        y_mid = r_mid * np.sin(phi_mid)

        x_end = r_end * np.cos(phi_end)
        y_end = r_end * np.sin(phi_end)

        x_center = r_center * np.cos(phi_center)
        y_center = r_center * np.sin(phi_center)

        # Calculate vectors from center to start and mid points
        vec_start = np.array(
            [x_start - x_center, y_start - y_center, z_start - z_center]
        )
        vec_mid = np.array([x_mid - x_center, y_mid - y_center, z_mid - z_center])
        vec_end = np.array([x_end - x_center, y_end - y_center, z_end - z_center])

        radius = np.linalg.norm(vec_start)

        if radius < 1e-10:
            logger.warning(f"Arc element {idx} has zero radius, skipping")
            return False

        # Calculate the normal to the plane of the arc (binormal)
        binormal = np.cross(vec_start, vec_mid)
        binormal_norm = np.linalg.norm(binormal)

        if binormal_norm < 1e-10:
            logger.warning(f"Arc element {idx} has collinear points, skipping")
            return False

        binormal = binormal / binormal_norm

        # Calculate angles
        # Normalize vectors
        vec_start_norm = vec_start / np.linalg.norm(vec_start)
        vec_mid_norm = vec_mid / np.linalg.norm(vec_mid)
        vec_end_norm = vec_end / np.linalg.norm(vec_end)

        # Calculate angle from start to end going through mid
        # Use atan2 for proper quadrant handling
        def angle_in_plane(vec):
            """Calculate angle of vector in the plane defined by binormal"""
            # Project vector onto plane and get angle
            local_x = vec_start_norm
            local_y = np.cross(binormal, vec_start_norm)
            x_comp = np.dot(vec, local_x)
            y_comp = np.dot(vec, local_y)
            return np.arctan2(y_comp, x_comp)

        angle_start = 0.0
        angle_mid = angle_in_plane(vec_mid_norm)
        angle_end = angle_in_plane(vec_end_norm)

        # Ensure we go through the intermediate point (aperture < pi)
        if angle_mid < 0:
            angle_mid += 2 * np.pi
        if angle_end < 0:
            angle_end += 2 * np.pi

        # If end angle is less than mid, we need to wrap around
        if angle_end < angle_mid:
            angle_end += 2 * np.pi

        # Create arc points
        cells.InsertNextCell(resolution + 1)
        for i in range(resolution + 1):
            t = i / resolution
            angle = angle_start + t * (angle_end - angle_start)

            # Rotate vec_start around binormal by angle
            cos_a = np.cos(angle)
            sin_a = np.sin(angle)

            # Rodrigues' rotation formula
            vec_rotated = (
                vec_start * cos_a
                + np.cross(binormal, vec_start) * sin_a
                + binormal * np.dot(binormal, vec_start) * (1 - cos_a)
            )

            # Add center offset
            x = x_center + vec_rotated[0]
            y = y_center + vec_rotated[1]
            z = z_center + vec_rotated[2]

            pts.InsertNextPoint(x, y, z)
            cells.InsertCellPoint(pts.GetNumberOfPoints() - 1)

        return True
