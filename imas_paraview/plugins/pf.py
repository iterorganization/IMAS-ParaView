"""Plugin to visualize the axisymmetric active poloidal field coils from the pf_active
IDS, and the axisymmetric passive conductors from the pf_active IDS."""

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

SUPPORTED_IDS_NAMES = ["pf_active", "pf_passive"]

# TODO: add documentation


@smproxy.source(label="PF Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class PFReader(GGDVTKPluginBase):
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
        """Select which coils or loops names to show in the array domain selector."""
        assert self._ids is not None, "IDS cannot be empty during setup."

        self.selectable_map = {}

        if self._ids.metadata.name == "pf_active":
            self._load_ids_quantities(self._ids.coil)
        elif self._ids.metadata.name == "pf_passive":
            self._load_ids_quantities(self._ids.loop, default_name="Loop")
        else:
            raise NotImplementedError(f"Unable to load {self._ids.metadata.name}.")

    def _load_ids_quantities(self, ids_quantity, default_name="Coil"):
        """Populate selectable elements from IDS content.

        Args:
            ids_quantity: Coil or loop array of structure of an IDS.
            default_name: Default name of the quantity.
        """
        for i, quantity in enumerate(ids_quantity):
            quantity_name = quantity.name
            if not quantity_name:
                quantity_name = f"{default_name} {i}"
                logger.warning(
                    f"{default_name} without name found. "
                    f"Renaming it to {quantity_name!r}"
                )
            if len(quantity.element) == 0:
                logger.warning(f"{quantity_name!r} has no elements, skipping it.")
                continue

            self.selectable_map[str(quantity_name)] = quantity
        self._selectable = list(self.selectable_map.keys())

    def _convert_to_vtk(self, output: vtkMultiBlockDataSet):
        """Convert each selected IDS quantity into a vtk object based on its geometry
        type, and store in as a separate block in the output vtkMultiBlockDataSet.

        Args:
            output: The vtkMultiBlockDataSet containing the converted vtk objects.
        """
        block_id = 0
        for quantity_name in self._selected:
            quantity = self.selectable_map[quantity_name]
            for element in quantity.element:
                geom_type = element.geometry.geometry_type

                if geom_type == 2:
                    vtk_coil = self._create_rectangle(element.geometry.rectangle)
                elif geom_type == 3:
                    vtk_coil = self._create_oblique(element.geometry.oblique)
                elif geom_type == 5:
                    vtk_coil = self._create_annulus(element.geometry.annulus)
                elif geom_type == 6:
                    vtk_coil = self._create_thick_line(element.geometry.thick_line)
                else:
                    logger.warning(
                        f"{quantity_name!r} has unsupported geometry type: {geom_type},"
                        " it will be skipped."
                    )
                    continue

                output.SetBlock(block_id, vtk_coil)
                block_id += 1

            logger.info(
                f"Loaded {quantity_name!r} with {len(quantity.element)} element(s)."
            )

    def _polyline_from_points(self, points):
        """Create a vtk polyline by connecting a list of points.

        Args:
            points: List of tuples containing x,y,z-coordinates of the points.
        """
        pts = vtkPoints()
        cells = vtkCellArray()
        cells.InsertNextCell(len(points))
        for x_i, y_i, z_i in points:
            pts.InsertNextPoint(x_i, y_i, z_i)
            cells.InsertCellPoint(pts.GetNumberOfPoints() - 1)
        poly = vtkPolyData()
        poly.SetPoints(pts)
        poly.SetLines(cells)
        return poly

    def _create_rectangle(self, rectangle):
        """Create vtkPolyData object from rectangle geometry"""
        r, z = rectangle.r, rectangle.z
        width, height = rectangle.width, rectangle.height

        r0 = r - width / 2.0
        r1 = r + width / 2.0
        z0 = z - height / 2.0
        z1 = z + height / 2.0

        points = [
            (r0, 0.0, z0),
            (r1, 0.0, z0),
            (r1, 0.0, z1),
            (r0, 0.0, z1),
            (r0, 0.0, z0),  # close the loop
        ]
        return self._polyline_from_points(points)

    def _create_oblique(self, oblique):
        """Create vtkPolyData object from oblique geometry"""
        r0, z0 = oblique.r, oblique.z
        la, lb = oblique.length_alpha, oblique.length_beta
        alpha, beta = oblique.alpha, oblique.beta

        dr_alpha = la * np.cos(alpha)
        dz_alpha = la * np.sin(alpha)

        dr_beta = -lb * np.sin(beta)
        dz_beta = lb * np.cos(beta)

        points = [
            (r0, 0.0, z0),
            (r0 + dr_alpha, 0.0, z0 + dz_alpha),
            (r0 + dr_alpha + dr_beta, 0.0, z0 + dz_alpha + dz_beta),
            (r0 + dr_beta, 0.0, z0 + dz_beta),
            (r0, 0.0, z0),  # close the loop
        ]
        return self._polyline_from_points(points)

    def _create_annulus(self, annulus, resolution=10):
        """Create vtkPolyData object from annulus geometry"""
        r0, z0 = annulus.r, annulus.z
        r_in, r_out = annulus.radius_inner, annulus.radius_outer

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

        polydata = vtkPolyData()
        polydata.SetPoints(pts)
        polydata.SetLines(cells)
        return polydata

    def _create_thick_line(self, thick_line):
        """Create vtkPolyData object from thick_line geometry"""
        p1 = thick_line.first_point
        p2 = thick_line.second_point
        thickness = thick_line.thickness

        dr = p2.r - p1.r
        dz = p2.z - p1.z
        length = np.hypot(dr, dz)
        if length == 0.0:
            logger.warning("Thick line with zero length, skipping")
            return None

        perp_r = -dz / length
        perp_z = dr / length

        offset_r = perp_r * thickness / 2.0
        offset_z = perp_z * thickness / 2.0

        corners = [
            (p1.r + offset_r, 0.0, p1.z + offset_z),
            (p2.r + offset_r, 0.0, p2.z + offset_z),
            (p2.r - offset_r, 0.0, p2.z - offset_z),
            (p1.r - offset_r, 0.0, p1.z - offset_z),
            (p1.r + offset_r, 0.0, p1.z + offset_z),  # close loop
        ]
        return self._polyline_from_points(corners)
