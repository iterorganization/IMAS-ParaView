"""Launch ParaView with the IMAS-ParaView plugins and environment preconfigured.

This resolves the ParaView executable and the environment variables that ParaView
needs to discover the IMAS-ParaView plugins (``PV_PLUGIN_PATH``), to import the
installed Python packages (``PYTHONPATH``) and to load data through the imas_core
HDF5 backend (``LD_PRELOAD``), so users do not have to set them up by hand.
"""
from packaging.version import Version
from collections import namedtuple

import logging
import os
import shutil
import sys
import sysconfig
from pathlib import Path

import imas_paraview
import imas_core
import subprocess
import re

logger = logging.getLogger("imas_paraview_launcher")


def find_paraview_binary():
    """Locate the ParaView executable.

    Resolution order:

    1. The ``PARAVIEW_BINARY`` environment variable, if set.
    2. A ``paraview`` executable under the interpreter's base prefix. This is the
       case when the virtual environment was created from ParaView's bundled
       ``pvpython`` interpreter, so ``sys.base_prefix`` points at ParaView's root.
    3. A ``paraview`` executable found on ``PATH``.

    Returns the path to the executable as a :class:`~pathlib.Path`, or ``None`` if
    none could be found.
    """
    override = os.environ.get("PARAVIEW_BINARY")
    if override:
        return Path(override)

    candidate = Path(sys.base_prefix) / "bin" / "paraview"
    if candidate.is_file():
        return candidate

    on_path = shutil.which("paraview")
    if on_path:
        return Path(on_path)

    return None


def plugin_path():
    """Return the directory containing the IMAS-ParaView ParaView plugins."""
    return str(Path(imas_paraview.__path__[0]) / "plugins")


HDF5Version = namedtuple('HDF5Version', ("path", "version"))

def get_hdf5_versions_for_elf(elf_path):
    result = subprocess.run(['ldd', str(elf_path)], capture_output=True, text=True, check=True)

    match = re.search(r'=> (.*libhdf5.*\.so\.([\d\.]+)) ', result.stdout)
    
    hdf5_version = "0.0.0"
    hdf5_path = Path()
    if match:
        hdf5_path = Path(match.group(1))
        hdf5_version = match.group(2)
    
    

    return HDF5Version(hdf5_path, Version(hdf5_version))



def imas_hdf5_version():
    imas_module_dir = Path(imas_core.__file__).parent
    al_core_objects = list(imas_module_dir.glob("_al*.so"))
    assert len(al_core_objects) == 1
    al_core_so = al_core_objects[0]

    return get_hdf5_versions_for_elf(al_core_so)


def paraview_hdf5_version():
    path = find_paraview_binary()
    
    return get_hdf5_versions_for_elf(path.parent / "paraview-real")


def hdf5_preload(site_packages):
    """Return the imas_core bundled HDF5 library to ``LD_PRELOAD``, or ``None``.

    The bundled library's filename contains a build-specific hash, so it is
    resolved with a glob rather than hardcoded.
    """

    imas_version = imas_hdf5_version()
    paraview_version = paraview_hdf5_version()

    if imas_version.version != paraview_version.version:
    
        logger.warning(
                "IMAS HDF5 version (%s) does not match ParaView HDF5 version (%s).\n"
                "  IMAS HDF5: %s\n"
                "  ParaView HDF5: %s\n"
                "Using the IMAS version, which might break ParaView's native HDF5 "
                "handling.",
                imas_version.version,
                paraview_version.version,
                imas_version.path,
                paraview_version.path,
            )
        return str(imas_version.path)

    return None


def build_environment(env=None):
    """Return a copy of ``env`` with the variables ParaView needs prepended.

    Existing values are preserved by prepending the IMAS-ParaView entries.
    """
    env = dict(os.environ if env is None else env)
    site_packages = sysconfig.get_path("purelib")

    def prepend(name, value):
        if value:
            existing = env.get(name)
            env[name] = os.pathsep.join(p for p in [value, existing] if p)

    prepend("PYTHONPATH", site_packages)
    prepend("PV_PLUGIN_PATH", plugin_path())
    prepend("LD_PRELOAD", hdf5_preload(site_packages))
    return env


def launch(paraview_args=()):
    """Set up the environment and replace the current process with ParaView.

    ``paraview_args`` are passed through to the ParaView executable.
    """
    paraview = find_paraview_binary()
    if paraview is None:
        raise SystemExit(
            "Could not locate the 'paraview' executable. Set the PARAVIEW_BINARY "
            "environment variable to its full path, or make sure 'paraview' is on "
            "your PATH."
        )

    env = build_environment()
    argv = [str(paraview), *paraview_args]
    logger.debug("Launching ParaView: %s", argv)
    os.execve(str(paraview), argv, env)
