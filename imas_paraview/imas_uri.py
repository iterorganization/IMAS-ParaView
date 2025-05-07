"""Helper functions for generating IMAS URIs."""

from pathlib import Path
from typing import Optional

import imas

BACKENDS = {
    imas.ids_defs.MDSPLUS_BACKEND: "mdsplus",
    imas.ids_defs.HDF5_BACKEND: "hdf5",
    imas.ids_defs.ASCII_BACKEND: "ascii",
}
"""Mapping of backend ID to the backend string used in the URI."""


def uri_from_path(path: Optional[str]) -> Optional[str]:
    """Create an IMAS URI from a local path.

    Args:
        path: Path to a local file. This should point to a valid file from the HDF5
            backend (with extension ``.h5``) or MDSplus backend (with extension
            ``.characteristics``, ``.datafile`` or ``.tree``).

    Returns:
        An IMAS URI pointing to the local file, or None if no URI could be constructed
        from the provided file.
    """
    uri = None
    if path:
        path = Path(path)
        if path.suffix == ".h5":
            uri = f"imas:hdf5?path={path.parent}"
        elif path.suffix in [".characteristics", ".datafile", ".tree"]:
            uri = f"imas:mdsplus?path={path.parent}"
        elif path.suffix == ".ids":
            uri = f"imas:ascii?path={path.parent}"
    return uri


def uri_from_pulse_run(
    backend: int, database: str, pulse: int, run: int, user: str, version: str
) -> Optional[str]:
    """Create an IMAS URI from 'legacy' parameters.

    Args:
        backend: Backend ID, e.g. HDF5_BACKEND
        database: Database name, e.g. "ITER"
        pulse: Pulse number of the database entry
        run: Run number of the database entry
        user: User name of the database
        version: Major version of the DD used by the the access layer

    Returns:
        An IMAS URI pointing to the data entry, or None if no valid URI could be
        constructed.
    """
    if backend not in BACKENDS:
        return None
    if not database or not user or not version:
        return None
    # TODO: use imas_core._al_lowlevel.al_build_uri_from_legacy_parameters instead?
    return (
        f"imas:{BACKENDS[backend]}?user={user};pulse={pulse};run={run};"
        f"database={database};version={version}"
    )
