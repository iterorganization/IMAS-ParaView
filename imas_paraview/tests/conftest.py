import logging
import os
import sys
from pathlib import Path

import imas
import pytest

from imas_paraview.plugins.vtkggdreader import SUPPORTED_IDS_NAMES
from imas_paraview.tests.fill_ggd import fill_ids


logger = logging.getLogger("imas_paraview")


def set_environment():
    """Set environment variables for PV_PLUGIN_PATH and PYTHONPATH."""
    current_path = Path(__file__).resolve().parent

    # Update PV_PLUGIN_PATH
    pv_plugin_path = os.environ.get("PV_PLUGIN_PATH", "")
    os.environ["PV_PLUGIN_PATH"] = f"{current_path.parent}/plugins:{pv_plugin_path}"

    # Update PYTHONPATH
    os.environ["PYTHONPATH"] = ":".join(sys.path)

    print("Setting the following environment variables:")
    print(f"PV_PLUGIN_PATH={os.environ['PV_PLUGIN_PATH']}")
    print(f"PYTHONPATH={os.environ['PYTHONPATH']}")


def pytest_sessionstart(session):
    """Set up the environment variables and generate test data to
    perform the integration tests on."""

    print("Setting up test environment...")
    set_environment()

    # Generate test data
    DD_VERSION = "3.42.0"
    ids = imas.IDSFactory(version=DD_VERSION).new("edge_profiles")
    fill_ids(ids, time_steps=10, grid_size=5)

    # Write test file as netCDF
    with imas.DBEntry("netcdf_testdb.nc", "w", dd_version=DD_VERSION) as dbentry:
        dbentry.put(ids)

    if not imas.backends.imas_core.imas_interface.has_imas:
        logger.warning(
            "IMAS-Core is not available, some integration tests are skipped."
        )
        return

    # Write test file as MDSPlus
    with imas.DBEntry(
        "imas:mdsplus?path=mdsplus_testdb", "w", dd_version=DD_VERSION
    ) as dbentry:
        dbentry.put(ids)

    # Write test file as HDF5
    with imas.DBEntry(
        "imas:hdf5?path=hdf5_testdb", "w", dd_version=DD_VERSION
    ) as dbentry:
        dbentry.put(ids)

    print("Test environment setup complete.")


@pytest.fixture(params=SUPPORTED_IDS_NAMES)
def ids_name(request):
    """Provides IDS names for testing."""

    return request.param


@pytest.fixture
def dummy_ids(ids_name):
    """Creates a dummy IDS object with a dummy grid and random GGD values for
    testing purposes."""

    ids = imas.IDSFactory().new(ids_name)
    fill_ids(ids)
    ids.validate()
    return ids


@pytest.fixture
def dummy_ids_five_steps(ids_name):
    """Creates a dummy IDS object containing five time steps with a dummy grid and
    random GGD values for testing purposes."""

    ids = imas.IDSFactory().new(ids_name)
    fill_ids(ids, time_steps=5)
    ids.validate()
    return ids
