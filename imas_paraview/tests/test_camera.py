import imas
from vtk.util.numpy_support import vtk_to_numpy

from imas_paraview.plugins.camera import CameraReader


def test_camera_visible():
    ids = imas.IDSFactory(version="4.1.0").new("camera_visible")
    name = "test_camera"

    ids.channel.resize(1)
    channel = ids.channel[0]
    channel.name = name

    channel.aperture.resize(1)
    aperture = channel.aperture[0]

    aperture.centre.r = 1.0
    aperture.centre.phi = 0.0
    aperture.centre.z = 0.0

    aperture.x3_unit_vector.x = 1.0
    aperture.x3_unit_vector.y = 0.0
    aperture.x3_unit_vector.z = 0.0

    aperture.x2_unit_vector.x = 0.0
    aperture.x2_unit_vector.y = 0.0
    aperture.x2_unit_vector.z = 1.0

    channel.viewing_angle_alpha_bounds = [-0.1, 0.1]
    channel.viewing_angle_beta_bounds = [-0.1, 0.1]

    reader = CameraReader()
    reader._ids = ids
    reader.setup_ids()

    assert len(reader.selectable_map) > 0

    geometry = reader.selectable_map[name]
    vtk_obj = reader._build_frustum_polydata(geometry)

    pts = vtk_to_numpy(vtk_obj.GetPoints().GetData())
    assert len(pts) > 0


def test_camera_ir():
    ids = imas.IDSFactory(version="4.1.0").new("camera_ir")

    ids.channel.resize(1)
    channel = ids.channel[0]
    channel.name = "test_channel"

    channel.camera.resize(1)
    camera = channel.camera[0]
    camera.name = "test_camera"

    camera.pinhole.x = 0.0
    camera.pinhole.y = 0.0
    camera.pinhole.z = 0.0

    camera.direction.x = 1.0
    camera.direction.y = 0.0
    camera.direction.z = 0.0

    camera.up.x = 0.0
    camera.up.y = 0.0
    camera.up.z = 1.0

    camera.field_of_view_horizontal = 0.2
    camera.field_of_view_vertical = 0.2

    channel.target_surface_center.x = 1.0
    channel.target_surface_center.y = 0.0
    channel.target_surface_center.z = 0.0

    reader = CameraReader()
    reader._ids = ids
    reader.setup_ids()

    assert len(reader.selectable_map) > 0

    geometry = reader.selectable_map["test_channel / test_camera"]
    vtk_obj = reader._build_frustum_polydata(geometry)

    pts = vtk_to_numpy(vtk_obj.GetPoints().GetData())
    assert len(pts) > 0
