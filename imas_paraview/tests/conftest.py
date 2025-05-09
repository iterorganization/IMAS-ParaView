import logging
import os
from pathlib import Path

import imas
import pytest

from imas_paraview.plugins.vtkggdreader import SUPPORTED_IDS_NAMES
from imas_paraview.tests.fill_ggd import fill_ids

logger = logging.getLogger("imas_paraview")


def set_environment():
    """Set environment variables for PV_PLUGIN_PATH and PYTHONPATH."""
    current_path = Path(__file__).resolve().parent

    # Set the LD_PRELOAD environment variable to handle Paraview bug
    hdf5_dir = os.getenv("HDF5_DIR")
    ld_preload_value = f"{hdf5_dir}/lib/libhdf5.so.310"
    os.environ["LD_PRELOAD"] = ld_preload_value

    # Update PV_PLUGIN_PATH
    pv_plugin_path = os.environ.get("PV_PLUGIN_PATH")
    if pv_plugin_path is None:
        os.environ["PV_PLUGIN_PATH"] = f"{current_path}/imas_paraview/plugins"
    else:
        os.environ["PV_PLUGIN_PATH"] = (
            f"{current_path}/imas_paraview/plugins:{pv_plugin_path}"
        )

    # Update PYTHONPATH
    python_path = os.environ.get("PYTHONPATH")
    if python_path is None:
        os.environ["PYTHONPATH"] = f"{current_path}"
    else:
        os.environ["PYTHONPATH"] = f"{current_path}:{python_path}"

    print("Setting the following environment variables:")
    print(f"PV_PLUGIN_PATH={os.environ['PV_PLUGIN_PATH']}")
    print(f"PYTHONPATH={os.environ['PYTHONPATH']}")
    print(f"LD_PRELOAD={os.environ['LD_PRELOAD']}")


def pytest_sessionstart(session):
    """Set up the environment variables and generate test data to
    perform the integration tests on."""

    if not imas.backends.imas_core.imas_interface.has_imas:
        logger.warning("IMAS-Core is not available, integration tests are skipped.")
        return
    print("Setting up test environment...")
    set_environment()

    # Generate test data
    ids = imas.IDSFactory(version="3.42.0").new("edge_profiles")
    fill_ids(ids, time_steps=10, grid_size=5)

    # Write test file as MDSPlus
    with imas.DBEntry("imas:mdsplus?path=mdsplus_testdb", "w") as dbentry:
        dbentry.put(ids)

    # Write test file as HDF5
    with imas.DBEntry("imas:hdf5?path=hdf5_testdb", "w") as dbentry:
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
