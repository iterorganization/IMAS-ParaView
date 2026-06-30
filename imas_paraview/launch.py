"""Launch ParaView with the IMAS-ParaView plugins and environment preconfigured.

This resolves the ParaView executable and the environment variables that ParaView
needs to discover the IMAS-ParaView plugins (``PV_PLUGIN_PATH``), to import the
installed Python packages (``PYTHONPATH``) and to load data through the imas_core
HDF5 backend (``LD_PRELOAD``), so users do not have to set them up by hand.
"""

import ctypes
import ctypes.util
import logging
import os
import shutil
import sys
import sysconfig
from collections import namedtuple
from pathlib import Path

import imas_core
from packaging.version import Version

import imas_paraview

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


HDF5Version = namedtuple("HDF5Version", ("path", "version"))


class Dl_info(ctypes.Structure):
    _fields_ = [
        ("dli_fname", ctypes.c_char_p),
        ("dli_fbase", ctypes.c_void_p),
        ("dli_sname", ctypes.c_char_p),
        ("dli_saddr", ctypes.c_void_p),
    ]


def get_function_origin(func):
    lib_path = ctypes.util.find_library("dl") or ctypes.util.find_library("c")
    if not lib_path:
        raise OSError("Could not find the system library containing dladdr.")

    sys_lib = ctypes.CDLL(lib_path)

    sys_lib.dladdr.argtypes = [ctypes.c_void_p, ctypes.POINTER(Dl_info)]
    sys_lib.dladdr.restype = ctypes.c_int

    func_ptr = ctypes.cast(func, ctypes.c_void_p)

    info = Dl_info()
    result = sys_lib.dladdr(func_ptr, ctypes.byref(info))

    if result != 0 and info.dli_fname:
        return info.dli_fname.decode("utf-8")
    else:
        return None


def get_hdf5_versions_for_elf(elf_path):
    hdf5_lib = ctypes.CDLL(elf_path)
    majnum = ctypes.c_uint()
    minnum = ctypes.c_uint()
    relnum = ctypes.c_uint()

    status = hdf5_lib.H5get_libversion(
        ctypes.byref(majnum), ctypes.byref(minnum), ctypes.byref(relnum)
    )
    if status < 0:
        raise RuntimeError("Failed to find HDF5 Version")

    return HDF5Version(
        Path(get_function_origin(hdf5_lib.H5get_libversion)),
        Version(f"{majnum.value}.{minnum.value}.{relnum.value}"),
    )


def imas_hdf5_version():
    al_core_so = imas_core._al_lowlevel.__file__

    return get_hdf5_versions_for_elf(al_core_so)


def hdf5_preload(site_packages):
    """Return the imas_core bundled HDF5 library to ``LD_PRELOAD``, or ``None``."""

    imas_version = imas_hdf5_version()

    logger.warning(
        "IMAS HDF5 version (%s) might not match ParaView HDF5 version.\n"
        "Using the IMAS version, which might break ParaView's native HDF5 "
        "handling.",
        imas_version.version,
    )
    return str(imas_version.path)


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
