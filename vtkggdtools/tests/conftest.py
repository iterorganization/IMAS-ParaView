import imaspy
import pytest

from vtkggdtools.tests.fill_ggd import fill_ids


@pytest.fixture(
    params=[
        "distribution_sources",
        "distributions",
        "edge_profiles",
        "edge_sources",
        "edge_transport",
        "equilibrium",
        "mhd",
        "radiation",
        "tf",
        "transport_solver_numerics",
        "wall",
        "waves",
        # "runaway_electrons",
    ]
)
def ids_name(request):
    """Provides IDS names for testing."""

    return request.param


@pytest.fixture
def dummy_ids(ids_name):
    """Creates a dummy IDS object with a dummy grid and random GGD values for
    testing purposes."""

    ids = imaspy.IDSFactory().new(ids_name)
    fill_ids(ids)
    return ids
