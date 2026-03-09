import imas
import numpy as np
from imas.ids_defs import IDS_TIME_MODE_HOMOGENEOUS
from vtk.util.numpy_support import vtk_to_numpy

from imas_paraview.plugins.pellets import PelletReader


def test_spi_fragments():
    ids = imas.IDSFactory(version="4.1.0").new("spi")
    ids.ids_properties.homogeneous_time = IDS_TIME_MODE_HOMOGENEOUS

    ids.time = np.array([0.0, 1.0, 2.0])

    ids.injector.resize(1)
    injector = ids.injector[0]
    injector.name = "test_injector"
    injector.fragment.resize(1)
    injector.fragment[0].position.r = [1.0, 2.0, 3.0]
    injector.fragment[0].position.phi = [4.0, 5.0, 6.0]
    injector.fragment[0].position.z = [7.0, 8.0, 9.0]
    injector.fragment[0].volume = [1.11, 2.22, 3.33]
    injector.fragment[0].velocity_r = [1.1, 2.2, 3.3]
    injector.fragment[0].velocity_phi = [4.4, 5.5, 6.6]
    injector.fragment[0].velocity_z = [7.7, 8.8, 9.9]

    # Fill shatter cone origin and velocity centre of mass
    injector.shatter_cone.origin.r = 0
    injector.shatter_cone.origin.phi = 0.0
    injector.shatter_cone.origin.z = 0.0
    injector.velocity_mass_centre_fragments_r = 1.0
    injector.velocity_mass_centre_fragments_phi = np.pi
    injector.velocity_mass_centre_fragments_z = 3.0

    reader = PelletReader()
    reader._ids = ids
    reader.setup_ids()

    fragments_pos = reader._create_fragment_positions(
        reader.selectable_map["Shattered Fragment Positions (test_injector)"], 0
    )
    assert len(vtk_to_numpy(fragments_pos.GetPoints().GetData())) > 0

    fragments_vel = reader._create_fragment_velocities(
        reader.selectable_map["Shattered Fragment Velocities (test_injector)"], 0
    )
    assert len(vtk_to_numpy(fragments_vel.GetPoints().GetData())) > 0

    fragments_com_vel = reader._create_centre_mass_velocity(
        reader.selectable_map["Fragment centre of mass velocity (test_injector)"],
    )
    assert len(vtk_to_numpy(fragments_com_vel.GetPoints().GetData())) > 0
