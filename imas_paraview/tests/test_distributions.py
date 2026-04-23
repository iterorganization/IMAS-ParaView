import imas
import numpy as np
from imas import identifiers
from imas.ids_defs import IDS_TIME_MODE_HOMOGENEOUS
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from imas_paraview.plugins.distributions import DistributionsReader


def test_distribution_name():
    ids = imas.IDSFactory(version="4.1.0").new("distributions")
    reader = DistributionsReader()
    ids.ids_properties.homogeneous_time = IDS_TIME_MODE_HOMOGENEOUS
    ids.distribution.resize(5)

    d = ids.distribution[0]
    d.species.type = identifiers.species_reference_identifier.electron
    assert reader._create_dist_name(d) == "Electron"

    d = ids.distribution[1]
    d.species.type = identifiers.species_reference_identifier.ion
    d.species.ion.name = "D+"
    assert reader._create_dist_name(d) == "Ion (D+)"

    d = ids.distribution[2]
    d.species.type = identifiers.species_reference_identifier.ion_state
    d.species.ion.name = "C+"
    d.species.ion.state.name = "C+2"
    assert reader._create_dist_name(d) == "Ion (C+) State (C+2)"

    d = ids.distribution[3]
    d.species.type = identifiers.species_reference_identifier.neutral
    d.species.neutral.name = "D"
    assert reader._create_dist_name(d) == "Neutral (D)"

    d = ids.distribution[4]
    d.species.type = identifiers.species_reference_identifier.neutral_state
    d.species.neutral.name = "H"
    d.species.neutral.state.name = "excited"
    assert reader._create_dist_name(d) == "Neutral (H) State (excited)"


def test_load_markers():
    ids = imas.IDSFactory(version="4.1.0").new("distributions")
    ids.ids_properties.homogeneous_time = IDS_TIME_MODE_HOMOGENEOUS
    ids.time = np.array([10.0])

    ids.distribution.resize(1)
    dist = ids.distribution[0]
    dist.species.type = identifiers.species_reference_identifier.electron
    dist.species.type.index = 1

    num_markers = 100
    dist.markers.resize(len(ids.time))

    dist.markers[0].coordinate_identifier.resize(3)
    dist.markers[0].coordinate_identifier[0] = identifiers.coordinate_identifier.r
    dist.markers[0].coordinate_identifier[1] = identifiers.coordinate_identifier.phi
    dist.markers[0].coordinate_identifier[2] = identifiers.coordinate_identifier.z

    dist.markers[0].positions = np.column_stack(
        [
            np.linspace(1.0, 3.0, num_markers),  # r
            np.linspace(0.0, 2 * np.pi, num_markers),  # phi
            np.linspace(-1.0, 1.0, num_markers),  # z
        ]
    )
    dist.markers[0].weights = np.random.rand(num_markers)

    reader = DistributionsReader()
    reader._ids = ids
    reader.setup_ids()

    output = vtkMultiBlockDataSet()
    reader._selected = ["Electron"]
    reader._load_markers(output, time_idx=0)

    assert output.GetNumberOfBlocks() == 1

    block = output.GetBlock(0)
    assert block is not None
    assert block.GetNumberOfPoints() == num_markers

    point_data = block.GetPointData()
    assert point_data.HasArray("r")
    assert point_data.HasArray("phi")
    assert point_data.HasArray("z")
    assert point_data.HasArray("weights")
