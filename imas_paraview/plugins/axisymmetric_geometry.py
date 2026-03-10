"""Plugin to visualize axisymmetric geometry structures."""

import logging

import numpy as np
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkMultiBlockDataSet,
    vtkPolyData,
)

from imas_paraview.ids_util import create_name_recursive
from imas_paraview.paraview_support.servermanager_tools import intvector, propertygroup
from imas_paraview.plugins.base_class import GGDVTKPluginBase
from imas_paraview.util import points_to_vtkpoly

logger = logging.getLogger("imas_paraview")

SUPPORTED_IDS_NAMES = [
    "pf_active",  # coil(i1)/element(i2)/geometry
    "pf_passive",  # loop(i1)/element(i2)/geometry
    "ferritic",  # object(i1)/axisymmetric(i2)
    "ic_antennas",  # antenna(i1)/module(i2)/strap(i3)/geometry
    "iron_core",  # segment(i1)/geometry
]


@smproxy.source(label="Geometry Reader (Axisymmetric)")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class AxisymmetricGeometryReader(GGDVTKPluginBase):
    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", SUPPORTED_IDS_NAMES)
        self.selectable_map = {}
        self.resolution = 10

    @intvector(label="Resolution", name="resolution", default_values=10)
    def P99_SetResolution(self, val):
        """Sets the number of points for the 'arcs_of_circle' and 'annulus' geometry
        types, if they are available in the loaded IDS."""
        self._update_property("resolution", val)

    @propertygroup("Geometry Reader (Axisymmetric) Settings", ["resolution"])
    def PG3_AxisymmetricGeometryReaderGroup(self):
        """Dummy function to define a PropertyGroup."""

    def setup_ids(self):
        """Select which coil or loop names to show in the array domain selector."""
        assert self._ids is not None, "IDS cannot be empty during setup."

        self.selectable_map = {}

        if self._ids.metadata.name not in SUPPORTED_IDS_NAMES:
            raise NotImplementedError(f"{self._ids.metadata.name} is not supported")

        if self._ids.metadata.name == "pf_active":
            for coil in self._ids.coil:
                for element in coil.element:
                    self._add_geometry_to_map(element.geometry)
        elif self._ids.metadata.name == "pf_passive":
            for loop in self._ids.loop:
                for element in loop.element:
                    self._add_geometry_to_map(element.geometry)
        elif self._ids.metadata.name == "ferritic":
            for obj in self._ids.object:
                for axisymmetric in obj.axisymmetric:
                    self._add_geometry_to_map(axisymmetric)
        elif self._ids.metadata.name == "ic_antennas":
            for antenna in self._ids.antenna:
                for module in antenna.module:
                    for strap in module.strap:
                        self._add_geometry_to_map(strap.geometry)
        elif self._ids.metadata.name == "iron_core":
            for segment in self._ids.segment:
                self._add_geometry_to_map(segment.geometry)

        self._selectable = list(self.selectable_map.keys())

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        if len(self._selected) > 0:
            output = vtkMultiBlockDataSet.GetData(outInfo)
            self._convert_to_vtk(output)
        return 1

    def _add_geometry_to_map(self, geometry):
        quantity_name = create_name_recursive(geometry)
        original_name = str(quantity_name)

        # Add suffix if given name is not unique
        if original_name in self.selectable_map:
            counter = 1
            while f"{original_name} #{counter}" in self.selectable_map:
                counter += 1
            quantity_name = f"{original_name} #{counter}"
        else:
            quantity_name = original_name

        self.selectable_map[str(quantity_name)] = geometry

    def _convert_to_vtk(self, output: vtkMultiBlockDataSet):
        """Convert each selected IDS quantity into a vtk object based on its geometry
        type, and store it as a separate block in the output vtkMultiBlockDataSet.

        Args:
            output: The vtkMultiBlockDataSet containing the converted vtk objects.
        """
        for block_id, quantity_name in enumerate(self._selected):
            quantity = self.selectable_map[quantity_name]
            geom_type = quantity.geometry_type

            if geom_type == 1:
                vtk_geom = self._create_outline(quantity.outline)
            elif geom_type == 2:
                vtk_geom = self._create_rectangle(quantity.rectangle)
            elif geom_type == 3:
                vtk_geom = self._create_oblique(quantity.oblique)
            elif geom_type == 4:
                vtk_geom = self._create_arcs_of_circle(quantity.arcs_of_circle)
            elif geom_type == 5:
                vtk_geom = self._create_annulus(quantity.annulus)
            elif geom_type == 6:
                vtk_geom = self._create_thick_line(quantity.thick_line)
            else:
                logger.warning(
                    "%r has unsupported geometry type: %s, it will be skipped.",
                    quantity_name,
                    geom_type,
                )
                continue

            output.SetBlock(block_id, vtk_geom)
            logger.info("Loaded '%s'.", quantity_name)

    def _create_outline(self, outline):
        """Create vtkPolyData object from outline geometry"""
        r, z = outline.r, outline.z
        points = [(r_i, 0.0, z_i) for r_i, z_i in zip(r, z, strict=True)]
        return points_to_vtkpoly(points, is_closed=True, is_filled=True)

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
        return points_to_vtkpoly(points, is_closed=True, is_filled=True)

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
        return points_to_vtkpoly(points, is_closed=True, is_filled=True)

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
        return points_to_vtkpoly(points, is_closed=True, is_filled=True)

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

        # Create a quad face for each segment
        idx = np.arange(self.resolution)
        idx_next = (idx + 1) % self.resolution

        cells = np.column_stack(
            [
                np.full(self.resolution, 4, dtype=np.int64),
                idx,
                idx + self.resolution,
                idx_next + self.resolution,
                idx_next,
            ]
        ).ravel()

        polys = vtkCellArray()
        polys.SetCells(self.resolution, numpy_to_vtkIdTypeArray(cells))

        polydata = vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetPolys(polys)
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
        return points_to_vtkpoly(points, is_closed=True, is_filled=True)

    def _stack_r_z(self, r, z):
        return np.column_stack((r, np.zeros_like(r), z))
