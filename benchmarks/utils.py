import imaspy
import imaspy.exception

HDF5 = "HDF5"
MDSPLUS = "MDSplus"
MEMORY = "memory"
ASCII = "ASCII"

all_backends = [
    HDF5,
    MDSPLUS,
    MEMORY,
    ASCII,
]


def create_uri(backend, path):
    return f"imas:{backend}?path={path}"


def create_dbentry(backend, path):
    return imaspy.DBEntry(create_uri(backend, path), "w")
