import pytest
from imas.ids_defs import ASCII_BACKEND, HDF5_BACKEND, MDSPLUS_BACKEND

from imas_paraview.imas_uri import uri_from_path, uri_from_pulse_run


def test_uri_from_path():
    # Test a couple of invalid paths
    assert uri_from_path(None) is None
    assert uri_from_path("") is None
    assert uri_from_path("xyz.test") is None
    assert uri_from_path("xyz") is None
    # Valid HDF5:
    assert uri_from_path("test/a.h5") == "imas:hdf5?path=test"
    assert uri_from_path("test/master.h5") == "imas:hdf5?path=test"
    # Valid MDSplus
    assert uri_from_path("test/ids_001.characteristics") == "imas:mdsplus?path=test"
    assert uri_from_path("test/ids_001.datafile") == "imas:mdsplus?path=test"
    assert uri_from_path("test/ids_001.tree") == "imas:mdsplus?path=test"


@pytest.mark.skip(reason="no IMAS-Core available")
def test_uri_from_pulse_run():
    # Missing values result in a None uri
    assert uri_from_pulse_run(-1, "a", 1, 2, "b", "c") is None
    assert uri_from_pulse_run(HDF5_BACKEND, "", 1, 2, "b", "c") is None
    assert uri_from_pulse_run(HDF5_BACKEND, "a", 1, 2, "", "c") is None
    assert uri_from_pulse_run(HDF5_BACKEND, "a", 1, 2, "b", "") is None

    # Valid parameters:
    assert (
        uri_from_pulse_run(HDF5_BACKEND, "db", 1, 2, "user", "3")
        == "imas:hdf5?user=user;pulse=1;run=2;database=db;version=3"
    )
    assert (
        uri_from_pulse_run(MDSPLUS_BACKEND, "db", 1, 2, "user", "4")
        == "imas:mdsplus?user=user;pulse=1;run=2;database=db;version=4"
    )
    assert (
        uri_from_pulse_run(ASCII_BACKEND, "db", 1, 2, "user", "3")
        == "imas:ascii?user=user;pulse=1;run=2;database=db;version=3"
    )
