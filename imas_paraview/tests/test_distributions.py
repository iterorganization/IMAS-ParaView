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
    d.species.type.name = "electron"
    d.species.type.index = 1
    assert reader._create_dist_name(d) == "electron"

    d = ids.distribution[1]
    d.species.type.name = "ion"
    d.species.type.index = 2
    d.species.ion.name = "D+"
    assert reader._create_dist_name(d) == "ion (D+)"

    d = ids.distribution[2]
    d.species.type.name = "ion"
    d.species.type.index = 3
    d.species.ion.name = "C+"
    d.species.ion.state.name = "C+2"
    assert reader._create_dist_name(d) == "ion (C+ (C+2))"

    d = ids.distribution[3]
    d.species.type.name = "neutral"
    d.species.type.index = 4
    d.species.neutral.name = "D"
    assert reader._create_dist_name(d) == "neutral (D)"

    d = ids.distribution[4]
    d.species.type.name = "neutral"
    d.species.type.index = 5
    d.species.neutral.name = "H"
    d.species.neutral.state.name = "excited"
    assert reader._create_dist_name(d) == "neutral (H (excited))"


def test_load_markers():
    ids = imas.IDSFactory(version="4.1.0").new("distributions")
    ids.ids_properties.homogeneous_time = IDS_TIME_MODE_HOMOGENEOUS
    ids.time = np.array([10.0])

    ids.distribution.resize(1)
    dist = ids.distribution[0]
    dist.species.type.name = "electron"
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
    reader._selected = ["electron"]
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

    assert point_data.HasArray("x")
    assert point_data.HasArray("y")
