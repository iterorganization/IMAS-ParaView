import logging
import math
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
from imas_paraview.util import pol_to_cart

logger = logging.getLogger("imas_paraview")

# TODO: add support for camera_x_rays and spectrometer_x_ray_crystal IDSs
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

        arr = vtkStringArray()
        for name in self.selectable_map:
            arr.InsertNextValue(name)
        return arr

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
        """Distance from the origin to the base of the frustum, in meters."""
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

        ids_name = self._ids.metadata.name
        if ids_name not in SUPPORTED_IDS_NAMES:
            raise NotImplementedError(f"{ids_name} is not supported")

        self.selectable_map = {}

        if ids_name == "camera_ir":
            self._setup_camera_ir()
        else:  # camera_visible
            self._setup_camera_visible()

        self._selectable = list(self.selectable_map.keys())

    def _setup_camera_ir(self):
        for ch_idx, channel in enumerate(self._ids.channel):
            channel_name = str(channel.name) or f"channel {ch_idx}"

            for cam_idx, camera in enumerate(channel.camera):
                camera_name = str(camera.name) or f"camera {cam_idx}"
                geom = CameraGeometry(
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

                entry_name = self._unique_name(f"{channel_name} / {camera_name}")
                self.selectable_map[entry_name] = geom

    def _setup_camera_visible(self):
        for ch_idx, channel in enumerate(self._ids.channel):
            channel_name = str(channel.name) or f"channel {ch_idx}"
            geom = self._extract_camera_visible_geometry(channel)
            if geom is None:
                continue

            self.selectable_map[self._unique_name(channel_name)] = geom

    def _extract_camera_visible_geometry(self, channel) -> CameraGeometry | None:
        if len(channel.aperture) == 0:
            logger.warning(
                "channel_visible channel '%s' has no aperture defined — skipping.",
                channel.name,
            )
            return None

        ap = channel.aperture[0]
        centre = ap.centre

        # general util for pol_to_cart3d?
        x, y = pol_to_cart(centre.r, centre.phi)

        # TODO: add geometry type? outline/circular/rectangle detector surface
        alpha = channel.viewing_angle_alpha_bounds
        beta = channel.viewing_angle_beta_bounds

        return CameraGeometry(
            origin=np.array([x, y, centre.z]),
            forward=np.array(
                [ap.x3_unit_vector.x, ap.x3_unit_vector.y, ap.x3_unit_vector.z]
            ),
            up=np.array(
                [ap.x2_unit_vector.x, ap.x2_unit_vector.y, ap.x2_unit_vector.z]
            ),
            hfov=alpha[1] - alpha[0],
            vfov=beta[1] - beta[0],
            target=None,
        )

    def _unique_name(self, name):
        # TODO: move to general util function?
        if name not in self.selectable_map:
            return name
        counter = 1
        while f"{name} #{counter}" in self.selectable_map:
            counter += 1
        return f"{name} #{counter}"

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        if self._selected:
            output = vtkMultiBlockDataSet.GetData(outInfo)
            self._convert_to_vtk(output)

        if self._snap_requested:
            self._snap_requested = False
            name = str(self._snap_camera_name).strip()
            logger.info("Snap requested for camera %r", name)
            logger.info("Available cameras: %s", list(self.selectable_map.keys()))
            geom = self.selectable_map.get(name)
            if geom is None:
                logger.warning("Snap: camera %r not found in selectable_map.", name)
            else:
                self._snap_view_to_geometry(geom)

        return 1

    def _convert_to_vtk(self, output: vtkMultiBlockDataSet):
        for block_id, entry_name in enumerate(self._selected):
            geom = self.selectable_map.get(entry_name)
            if geom is None:
                logger.warning(
                    "Selected camera '%s' not found in selectable_map — skipping.",
                    entry_name,
                )
                continue

            vtk_geom = self._build_frustum_polydata(
                geom, self.frustum_length, entry_name
            )
            if vtk_geom is not None:
                output.SetBlock(block_id, vtk_geom)
                meta = output.GetMetaData(block_id)
                meta.Set(vtkCompositeDataSet.NAME(), entry_name)
                logger.info("Loaded frustum for '%s'.", entry_name)

    def _build_frustum_polydata(
        self,
        geom: CameraGeometry,
        frustum_length: float,
        entry_name: str,
    ) -> vtkPolyData | None:
        forward = geom.forward / np.linalg.norm(geom.forward)

        up_raw = geom.up / np.linalg.norm(geom.up)
        up_ortho = up_raw - np.dot(up_raw, forward) * forward
        if np.linalg.norm(up_ortho) < 1e-12:
            logger.warning(
                "'%s' up vector is parallel to direction — skipping.", entry_name
            )
            return None
        up_ortho /= np.linalg.norm(up_ortho)

        right = np.cross(forward, up_ortho)
        right /= np.linalg.norm(right)

        d = frustum_length
        half_h = d * math.tan(geom.hfov / 2.0)
        half_v = d * math.tan(geom.vfov / 2.0)

        base_center = geom.origin + d * forward

        c0 = base_center - half_h * right - half_v * up_ortho
        c1 = base_center + half_h * right - half_v * up_ortho
        c2 = base_center + half_h * right + half_v * up_ortho
        c3 = base_center - half_h * right + half_v * up_ortho

        pts_np = np.array([geom.origin, c0, c1, c2, c3], dtype=np.float64)

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

    def _snap_view_to_geometry(self, geom: CameraGeometry):
        view = GetActiveView()
        if view is None:
            logger.warning("Snap: no active view found.")
            return

        forward = geom.forward / np.linalg.norm(geom.forward)
        up_vector = geom.up / np.linalg.norm(geom.up)

        focal = (
            np.array(geom.target, dtype=float)
            if geom.target is not None
            else geom.origin + forward
        )

        # NOTE: Would be nice to also change the aspect ratio of the renderview to
        # match the camera. Trying to set the view's ViewSize doesn't seem to work,
        # it just stretches the resolution to fit into viewport, distorting the UI.

        vtk_cam = view.GetActiveCamera()
        vtk_cam.SetPosition(*geom.origin)
        vtk_cam.SetFocalPoint(*focal)
        vtk_cam.SetViewUp(*up_vector)
        vtk_cam.SetViewAngle(math.degrees(geom.vfov))
        view.StillRender()
        logger.info(
            "Snapped to %s: Position=%s, ViewUp=%s, VFOV=%.2f deg",
            self._snap_camera_name,
            geom.origin,
            up_vector,
            math.degrees(geom.vfov),
        )
