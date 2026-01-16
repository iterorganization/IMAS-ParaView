"""Plugin to visualize the axisymmetric active poloidal field coils in the pf_active
IDS, as well as the axisymmetric passive conductors in the pf_passive IDS."""

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
from imas_paraview.util import points_to_vtkpoly

logger = logging.getLogger("imas_paraview")

SUPPORTED_IDS_NAMES = ["pf_active", "pf_passive"]


@smproxy.source(label="PF Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class PFReader(GGDVTKPluginBase):
    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", SUPPORTED_IDS_NAMES)
        self.selectable_map = {}
        self.resolution = 10

    @intvector(label="Resolution", name="resolution", default_values=10)
    def P99_SetResolution(self, val):
        """Sets the number of points for the 'arcs_of_circle' and 'annulus' geometry
        types, if they are available in the loaded IDS."""
        self._update_property("resolution", val)

    @propertygroup("PF Reader Settings", ["resolution"])
    def PG3_PFReaderGroup(self):
        """Dummy function to define a PropertyGroup."""

    def setup_ids(self):
        """Select which coil or loop names to show in the array domain selector."""
        assert self._ids is not None, "IDS cannot be empty during setup."

        self.selectable_map = {}

        if self._ids.metadata.name == "pf_active":
            self._load_ids_quantities(self._ids.coil, "Coil")
        elif self._ids.metadata.name == "pf_passive":
            self._load_ids_quantities(self._ids.loop, "Loop")
        else:
            raise NotImplementedError(f"Unable to load {self._ids.metadata.name}.")

    def _load_ids_quantities(self, ids_quantity, default_name):
        """Populate selectable elements from the IDS content.

        Args:
            ids_quantity: Coil or loop AoS of an IDS.
            default_name: Default name of the quantity.
        """
        for i, quantity in enumerate(ids_quantity):
            quantity_name = quantity.name
            if not quantity_name:
                quantity_name = f"{default_name} {i}"
                logger.warning(
                    "%s without name found. Renaming it to %r",
                    default_name,
                    quantity_name,
                )
            if len(quantity.element) == 0:
                logger.warning("%r has no elements, skipping it.", quantity_name)
                continue

            self.selectable_map[str(quantity_name)] = quantity
        self._selectable = list(self.selectable_map.keys())

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        if len(self._selected) > 0:
            output = vtkMultiBlockDataSet.GetData(outInfo)
            self._convert_to_vtk(output)
        return 1

    def _convert_to_vtk(self, output: vtkMultiBlockDataSet):
        """Convert each selected IDS quantity into a vtk object based on its geometry
        type, and store it as a separate block in the output vtkMultiBlockDataSet.

        Args:
            output: The vtkMultiBlockDataSet containing the converted vtk objects.
        """
        block_id = 0
        for quantity_name in self._selected:
            quantity = self.selectable_map[quantity_name]
            for element in quantity.element:
                geom_type = element.geometry.geometry_type

                if geom_type == 1:
                    vtk_geom = self._create_outline(element.geometry.outline)
                elif geom_type == 2:
                    vtk_geom = self._create_rectangle(element.geometry.rectangle)
                elif geom_type == 3:
                    vtk_geom = self._create_oblique(element.geometry.oblique)
                elif geom_type == 4:
                    vtk_geom = self._create_arcs_of_circle(
                        element.geometry.arcs_of_circle
                    )
                elif geom_type == 5:
                    vtk_geom = self._create_annulus(element.geometry.annulus)
                elif geom_type == 6:
                    vtk_geom = self._create_thick_line(element.geometry.thick_line)
                else:
                    logger.warning(
                        "%r has unsupported geometry type: %s, it will be skipped.",
                        quantity_name,
                        geom_type,
                    )
                    continue

                output.SetBlock(block_id, vtk_geom)
                block_id += 1

            logger.info(
                "Loaded %r with %s element(s).", quantity_name, len(quantity.element)
            )

    def _create_outline(self, outline):
        """Create vtkPolyData object from outline geometry"""
        r, z = outline.r, outline.z
        points = [(r_i, 0.0, z_i) for r_i, z_i in zip(r, z, strict=True)]
        return points_to_vtkpoly(points, is_closed=True)

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
        ]
        return points_to_vtkpoly(points, is_closed=True)

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
        ]
        return points_to_vtkpoly(points, is_closed=True)

    def _create_arcs_of_circle(self, arcs_of_circle):
        """Create vtkPolyData object from arcs_of_circle geometry."""
        r = arcs_of_circle.r
        z = arcs_of_circle.z
        radii = arcs_of_circle.curvature_radii

        r1, z1 = r, z
        r2, z2 = np.roll(r, -1), np.roll(z, -1)

        radius = np.abs(radii)
        curvature_sign = np.sign(radii)

        dr = r2 - r1
        dz = z2 - z1

        # Arc chord length
        d = np.hypot(dr, dz)

        if np.any(d > 2.0 * radius):
            logger.warning("Arc chord is longer than diameter, skipping element")
            return None

        # Midpoint of chords
        mr = (r1 + r2) / 2
        mz = (z1 + z2) / 2

        # Distance from midpoint to circle center
        h = np.sqrt(radius**2 - (d / 2.0) ** 2)

        # Unit perpendicular vectors
        nr = -dz / d
        nz = dr / d

        # Centers of circles
        cr = mr + curvature_sign * h * nr
        cz = mz + curvature_sign * h * nz

        # Angles of endpoints relative to centers
        theta1 = np.arctan2(z1 - cz, r1 - cr)
        theta2 = np.arctan2(z2 - cz, r2 - cr)

        mask_pos = (curvature_sign > 0) & (theta2 < theta1)
        mask_neg = (curvature_sign < 0) & (theta2 > theta1)
        theta2[mask_pos] += 2 * np.pi
        theta2[mask_neg] -= 2 * np.pi

        theta = np.linspace(theta1, theta2, self.resolution, endpoint=False).T

        points_r = cr[:, None] + radius[:, None] * np.cos(theta)
        points_z = cz[:, None] + radius[:, None] * np.sin(theta)

        points = self._stack_r_z(points_r.ravel(), points_z.ravel())
        return points_to_vtkpoly(points, is_closed=True)

    def _create_annulus(self, annulus):
        """Create vtkPolyData object from annulus geometry"""
        r0, z0 = annulus.r, annulus.z
        r_in, r_out = annulus.radius_inner, annulus.radius_outer

        theta = np.linspace(0, 2.0 * np.pi, self.resolution, endpoint=False)

        inner_points = self._stack_r_z(
            r0 + r_in * np.cos(theta), z0 + r_in * np.sin(theta)
        )
        outer_points = self._stack_r_z(
            r0 + r_out * np.cos(theta), z0 + r_out * np.sin(theta)
        )

        points = vtkPoints()
        points.SetData(numpy_to_vtk(np.vstack([inner_points, outer_points])))

        inner_ids = np.arange(self.resolution)
        outer_ids = np.arange(self.resolution, 2 * self.resolution)

        inner_loop = np.append(inner_ids, inner_ids[0])
        outer_loop = np.append(outer_ids, outer_ids[0])

        cells = vtkCellArray()
        cells.InsertNextCell(self.resolution + 1, inner_loop)
        cells.InsertNextCell(self.resolution + 1, outer_loop)

        polydata = vtkPolyData()
        polydata.SetPoints(points)
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

        points = [
            (p1.r + offset_r, 0.0, p1.z + offset_z),
            (p2.r + offset_r, 0.0, p2.z + offset_z),
            (p2.r - offset_r, 0.0, p2.z - offset_z),
            (p1.r - offset_r, 0.0, p1.z - offset_z),
        ]
        return points_to_vtkpoly(points, is_closed=True)

    def _stack_r_z(self, r, z):
        return np.column_stack((r, np.zeros_like(r), z))
