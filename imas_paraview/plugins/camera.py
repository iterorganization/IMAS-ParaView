import logging
from dataclasses import dataclass

import numpy as np
from paraview.simple import (
    GetActiveView,
)
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
from vtkmodules.vtkCommonCore import vtkPoints
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
    target: np.ndarray | None


@smproxy.source(label="Camera Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class CameraReader(GGDVTKPluginBase):
    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", SUPPORTED_IDS_NAMES)
        self.frustum_length = 1.0
        self.selectable_map: dict[str, CameraGeometry] = {}
        self._snap_camera_name = ""
        self._snap_requested = False

    @stringvector(
        name="SnapCameraList", information_only=1, si_class="vtkSIDataArrayProperty"
    )
    def P96_GetSnapCameraList(self):
        """Return the list of loaded cameras for the snap dropdown."""
        from vtkmodules.vtkCommonCore import vtkStringArray

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
        self._snap_requested = True
        self.Modified()

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
            self._setup_camera_ir()
        else:  # camera_visible
            self._setup_camera_visible()

        self._selectable = list(self.selectable_map.keys())

    def _setup_camera_ir(self):
        # NOTE: This requires DD version >= 4.1.0
        for ch_idx, channel in enumerate(self._ids.channel):
            channel_name = str(channel.name) or f"channel {ch_idx}"

            for cam_idx, camera in enumerate(channel.camera):
                camera_name = str(camera.name) or f"camera {cam_idx}"
                geometry = CameraGeometry(
                    origin=np.array(
                        [camera.pinhole.x, camera.pinhole.y, camera.pinhole.z]
                    ),
                    forward=np.array(
                        [camera.direction.x, camera.direction.y, camera.direction.z]
                    ),
                    up=np.array([camera.up.x, camera.up.y, camera.up.z]),
                    hfov=camera.field_of_view_horizontal,
                    vfov=camera.field_of_view_vertical,
                    target=np.array(
                        [
                            channel.target_surface_center.x,
                            channel.target_surface_center.y,
                            channel.target_surface_center.z,
                        ]
                    ),
                )

                camera_name = ensure_unique_name(
                    f"{channel_name} / {camera_name}", list(self.selectable_map.keys())
                )
                self.selectable_map[camera_name] = geometry

    def _setup_camera_visible(self):
        for channel in self._ids.channel:
            name = ensure_unique_name(
                str(channel.name), list(self.selectable_map.keys())
            )
            if len(channel.aperture) == 0:
                logger.warning("'%s' has no aperture defined, skipping.", name)
                continue

            aperture = channel.aperture[0]
            centre = aperture.centre

            x, y = pol_to_cart(centre.r, centre.phi)

            alpha = channel.viewing_angle_alpha_bounds
            beta = channel.viewing_angle_beta_bounds

            geometry = CameraGeometry(
                origin=np.array([x, y, centre.z]),
                forward=np.array(
                    [
                        aperture.x3_unit_vector.x,
                        aperture.x3_unit_vector.y,
                        aperture.x3_unit_vector.z,
                    ]
                ),
                up=np.array(
                    [
                        aperture.x2_unit_vector.x,
                        aperture.x2_unit_vector.y,
                        aperture.x2_unit_vector.z,
                    ]
                ),
                hfov=alpha[1] - alpha[0],
                vfov=beta[1] - beta[0],
                target=None,
            )

            self.selectable_map[name] = geometry

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        if self._selected:
            output = vtkMultiBlockDataSet.GetData(outInfo)
            self._convert_to_vtk(output)

        if self._snap_requested:
            self._snap_requested = False
            logger.info("Snapping to camera '%s'", self._snap_camera_name)
            geometry = self.selectable_map[self._snap_camera_name]
            self._snap_view_to_geometry(geometry)

        return 1

    def _convert_to_vtk(self, output: vtkMultiBlockDataSet):
        for block_id, name in enumerate(self._selected):
            geometry = self.selectable_map[name]
            vtk_geom = self._build_frustum_polydata(geometry)
            if vtk_geom is not None:
                output.SetBlock(block_id, vtk_geom)
                meta = output.GetMetaData(block_id)
                meta.Set(vtkCompositeDataSet.NAME(), name)
                logger.info("Loaded frustum for '%s'.", name)

    def _build_frustum_polydata(self, geometry: CameraGeometry):
        forward = geometry.forward / np.linalg.norm(geometry.forward)

        up_raw = geometry.up / np.linalg.norm(geometry.up)
        up_ortho = up_raw - np.dot(up_raw, forward) * forward
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

        pts_np = np.array([geometry.origin, c0, c1, c2, c3], dtype=np.float64)

        vtk_pts = vtkPoints()
        vtk_pts.SetData(numpy_to_vtk(pts_np))

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

        flat = np.array([v for edge in edges for v in edge], dtype=np.int64)
        lines = vtkCellArray()
        lines.SetCells(len(edges), numpy_to_vtkIdTypeArray(flat))

        polydata = vtkPolyData()
        polydata.SetPoints(vtk_pts)
        polydata.SetLines(lines)
        return polydata

    def _snap_view_to_geometry(self, geometry: CameraGeometry):
        view = GetActiveView()
        if view.GetXMLName() != "RenderView":
            logger.error(
                "Cannot snap in current active viewport. Please select a RenderView as"
                "your active view."
            )
            return

        target = (
            geometry.target
            if geometry.target is not None
            else geometry.origin + geometry.forward
        )

        vtk_cam = view.GetActiveCamera()
        vtk_cam.SetPosition(*geometry.origin)
        vtk_cam.SetFocalPoint(*target)
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
