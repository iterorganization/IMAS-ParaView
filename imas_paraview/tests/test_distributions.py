import imas
from imas.ids_defs import IDS_TIME_MODE_HOMOGENEOUS

from imas_paraview.plugins.distributions import DistributionsMarkersReader


def test_distribution_name():
    ids = imas.IDSFactory(version="4.1.0").new("distributions")
    reader = DistributionsMarkersReader()
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
