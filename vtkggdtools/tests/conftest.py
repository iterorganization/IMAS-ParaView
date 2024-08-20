import imaspy
import pytest

from vtkggdtools.io.read_ps import SUPPORTED_IDS_NAMES
from vtkggdtools.tests.fill_ggd import fill_ids


@pytest.fixture(params=SUPPORTED_IDS_NAMES)
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
