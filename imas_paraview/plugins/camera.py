import logging
from dataclasses import dataclass

import numpy as np
from packaging.version import Version
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
from vtkmodules.vtkCommonCore import vtkPoints, vtkStringArray
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkCompositeDataSet,
    vtkMultiBlockDataSet,
    vtkPolyData,
)

from imas_paraview.paraview_support.servermanager_tools import (
    command_button_property,
    doublevector,
    propertygroup,
    stringlistdomain,
    stringvector,
)
from imas_paraview.plugins.base_class import GGDVTKPluginBase
from imas_paraview.util import ensure_unique_name, pol_to_cart

logger = logging.getLogger("imas_paraview")

# TODO: add support for camera_x_rays IDS
SUPPORTED_IDS_NAMES = ["camera_ir", "camera_visible"]


@dataclass
class CameraGeometry:
    origin: np.ndarray
    forward: np.ndarray
    up: np.ndarray
    hfov: float
    vfov: float
    target: np.ndarray


@smproxy.source(label="Camera Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class CameraReader(GGDVTKPluginBase):
    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", SUPPORTED_IDS_NAMES)
        self.frustum_length = 1.0
        self.selectable_map: dict[str, CameraGeometry] = {}
        self._snap_camera_name = ""

    @stringvector(
        name="SnapCameraList", information_only=1, si_class="vtkSIDataArrayProperty"
    )
    def P96_GetSnapCameraList(self):
        """Return the list of loaded cameras for the snap dropdown."""

        array = vtkStringArray()
        for name in self.selectable_map:
            array.InsertNextValue(name)
        return array

    @stringvector(name="SnapCameraName", label="Select Camera")
    @stringlistdomain("SnapCameraList", name="snap_camera_list")
    def P97_SetSnapCameraName(self, value):
        """Select which loaded camera to snap the ParaView view to."""
        self._snap_camera_name = str(value).strip()

    @command_button_property("SnapToCamera", "Snap View to Camera", "P98_SnapToCamera")
    def P98_SnapToCamera(self):
        """Snap the ParaView camera to the selected camera."""
        if (
            not self._snap_camera_name
            or self._snap_camera_name not in self.selectable_map
        ):
            logger.error("No valid camera selected to snap to.")
            return

        geometry = self.selectable_map[self._snap_camera_name]
        self._snap_view_to_geometry(geometry)

    @doublevector(
        label="Frustum Length (m)",
        name="frustum_length",
        default_values=1.0,
    )
    def P99_SetFrustumLength(self, val):
        """Sets the distance from the camera's origin to the base of the frustum."""
        self._update_property("frustum_length", val)

    @propertygroup("Camera Reader Settings", ["frustum_length"])
    def PG3_CameraFrustumGroup(self):
        """Dummy function to define a PropertyGroup."""

    @propertygroup("Snap View to Camera", ["SnapCameraName", "SnapToCamera"])
    def PG4_SnapCameraGroup(self):
        """Dummy function to define a PropertyGroup."""

    def setup_ids(self):
        """Populate the selectable list with one entry per camera/channel."""
        assert self._ids is not None, "IDS cannot be empty during setup."
        self.selectable_map = {}

        if self._ids.metadata.name == "camera_ir":
            self._extract_camera_ir()
        else:  # camera_visible
            self._extract_camera_visible()

        self._selectable = list(self.selectable_map.keys())

    def _cart_vector_has_value(self, vec):
        """Helper to check if a lazy-loaded cartesian IMAS vector quantity has data."""
        try:
            return vec.x.has_value and vec.y.has_value and vec.z.has_value
        except AttributeError:
            return False

    def _extract_camera_ir(self):
        """Extract camera geometries from camera_ir IDS."""
        imas_version = str(self._ids.ids_properties.version_put.data_dictionary)
        if imas_version and Version(imas_version) < Version("4.1.0"):
            logger.error(
                "The DD version of the IDS ('%s') is too old, it should be at "
                "least '4.1.0'",
                imas_version,
            )
            return

        for ch_idx, channel in enumerate(self._ids.channel):
            channel_name = str(channel.name) or f"channel {ch_idx}"

            for cam_idx, camera in enumerate(channel.camera):
                camera_name = str(camera.name) or f"camera {cam_idx}"
                geom_name = ensure_unique_name(
                    f"{channel_name} / {camera_name}", list(self.selectable_map.keys())
                )
                if not (
                    self._cart_vector_has_value(camera.pinhole)
                    and self._cart_vector_has_value(camera.direction)
                    and self._cart_vector_has_value(camera.up)
                    and self._cart_vector_has_value(channel.target_surface_center)
                    and camera.field_of_view_horizontal.has_value
                    and camera.field_of_view_vertical.has_value
                ):
                    logger.warning(
                        "'%s' is missing required geometry data. Skipping.", geom_name
                    )
                    continue

                origin = np.array(
                    [camera.pinhole.x, camera.pinhole.y, camera.pinhole.z]
                )
                forward = np.array(
                    [camera.direction.x, camera.direction.y, camera.direction.z]
                )
                up = np.array([camera.up.x, camera.up.y, camera.up.z])
                target = np.array(
                    [
                        channel.target_surface_center.x,
                        channel.target_surface_center.y,
                        channel.target_surface_center.z,
                    ]
                )

                geometry = CameraGeometry(
                    origin=origin,
                    forward=forward,
                    up=up,
                    hfov=camera.field_of_view_horizontal,
                    vfov=camera.field_of_view_vertical,
                    target=target,
                )

                self.selectable_map[geom_name] = geometry

    def _extract_camera_visible(self):
        """Extract camera geometries from camera_visible IDS."""
        for channel in self._ids.channel:
            name = ensure_unique_name(
                str(channel.name), list(self.selectable_map.keys())
            )
            if len(channel.aperture) == 0:
                logger.warning("'%s' has no aperture defined, skipping.", name)
                continue

            ap = channel.aperture[0]
            centre = ap.centre
            alpha = channel.viewing_angle_alpha_bounds
            beta = channel.viewing_angle_beta_bounds

            if not (
                centre.r.has_value
                and centre.phi.has_value
                and centre.z.has_value
                and self._cart_vector_has_value(ap.x3_unit_vector)
                and self._cart_vector_has_value(ap.x2_unit_vector)
                and alpha.has_value
                and beta.has_value
            ):
                logger.warning(
                    "'%s' is missing required geometry data. Skipping.", name
                )
                continue

            x, y = pol_to_cart(centre.r, centre.phi)

            origin = np.array([x, y, centre.z])
            forward = np.array(
                [ap.x3_unit_vector.x, ap.x3_unit_vector.y, ap.x3_unit_vector.z]
            )
            up = np.array(
                [ap.x2_unit_vector.x, ap.x2_unit_vector.y, ap.x2_unit_vector.z]
            )

            self.selectable_map[name] = CameraGeometry(
                origin=origin,
                forward=forward,
                up=up,
                hfov=alpha[1] - alpha[0],
                vfov=beta[1] - beta[0],
                target=origin + forward,
            )

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        if self._selected:
            output = vtkMultiBlockDataSet.GetData(outInfo)
            self._convert_to_vtk(output)

        return 1

    def _convert_to_vtk(self, output: vtkMultiBlockDataSet):
        """Convert selected cameras into a vtkPolyData object and store it in a block of
        a vtkMultiBlockDataSet.

        Args:
            output: The vtkMultiBlockDataSet to store the camera into.
        """
        for block_id, name in enumerate(self._selected):
            geometry = self.selectable_map[name]
            vtk_geom = self._build_camera_polydata(geometry)
            if vtk_geom is not None:
                output.SetBlock(block_id, vtk_geom)
                meta = output.GetMetaData(block_id)
                meta.Set(vtkCompositeDataSet.NAME(), name)
                logger.info("Loaded frustum for '%s'.", name)

    def _build_camera_polydata(self, geometry: CameraGeometry):
        """Convert a CameraGeometry into a vtkPolyData object.

        Args:
            geometry: CameraGeometry dataclass of the selected camera.

        Returns:
            vtkPolyData of the view pyramid.
        """
        forward = geometry.forward / np.linalg.norm(geometry.forward)

        up = geometry.up / np.linalg.norm(geometry.up)
        up_ortho = up - np.dot(up, forward) * forward
        up_ortho /= np.linalg.norm(up_ortho)

        right = np.cross(forward, up_ortho)
        right /= np.linalg.norm(right)

        half_h = self.frustum_length * np.tan(geometry.hfov / 2.0)
        half_v = self.frustum_length * np.tan(geometry.vfov / 2.0)

        base_center = geometry.origin + self.frustum_length * forward

        c0 = base_center - half_h * right - half_v * up_ortho
        c1 = base_center + half_h * right - half_v * up_ortho
        c2 = base_center + half_h * right + half_v * up_ortho
        c3 = base_center - half_h * right + half_v * up_ortho

        points = np.array([geometry.origin, c0, c1, c2, c3])

        vtk_pts = vtkPoints()
        vtk_pts.SetData(numpy_to_vtk(points))

        edges = [
            [2, 0, 1],  # origin to c0
            [2, 0, 2],  # origin to c1
            [2, 0, 3],  # origin to c2
            [2, 0, 4],  # origin to c3
            [2, 1, 2],  # c0  to c1
            [2, 2, 3],  # c1  to c2
            [2, 3, 4],  # c2  to c3
            [2, 4, 1],  # c3  to c0
        ]

        flat = np.array(edges, dtype=np.int64).flatten()
        lines = vtkCellArray()
        lines.SetCells(len(edges), numpy_to_vtkIdTypeArray(flat))

        polydata = vtkPolyData()
        polydata.SetPoints(vtk_pts)
        polydata.SetLines(lines)
        return polydata

    def _snap_view_to_geometry(self, geometry: CameraGeometry):
        """Snap the currently active RenderView to the selected camera, setting
        ParaView's camera to match with the selected camera's position, focal point,
        and field of view.

        Args:
            geometry: The CameraGeometry to snap to.
        """

        from paraview.simple import GetActiveView

        view = GetActiveView()
        if view is None:
            logger.error(
                "Cannot find your active view. Snapping the camera view is only "
                "available when running ParaView in standalone mode."
            )
            return
        if view.GetXMLName() != "RenderView":
            logger.error(
                "Cannot snap camera in the currently active viewport. Please select a "
                "RenderView as your active view."
            )
            return

        vtk_cam = view.GetActiveCamera()
        vtk_cam.SetPosition(*geometry.origin)
        vtk_cam.SetFocalPoint(*geometry.target)

        # Ensure camera is oriented upright
        if geometry.up[2] < 0:
            vtk_cam.SetViewUp(*(-geometry.up))
        else:
            vtk_cam.SetViewUp(*geometry.up)

        # Ensure full camera view fits within the RenderView
        view_width, view_height = view.ViewSize
        viewport_aspect = view_width / view_height
        camera_aspect = np.tan(geometry.hfov / 2.0) / np.tan(geometry.vfov / 2.0)
        if viewport_aspect < camera_aspect:
            vfov_to_set = 2.0 * np.arctan(np.tan(geometry.hfov / 2.0) / viewport_aspect)
        else:
            vfov_to_set = geometry.vfov
        vtk_cam.SetViewAngle(np.degrees(vfov_to_set))

        view.StillRender()
        logger.info("Snapped camera to '%s'", self._snap_camera_name)
