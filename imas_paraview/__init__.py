import logging

import imas
from packaging import version

from imas_paraview.vtkhandler import VTKHandler

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "unknown"

# Setup logger
logger = logging.getLogger("imas_paraview")
handler = VTKHandler()
logger.addHandler(handler)
logger.setLevel(handler.get_level())

# Remove imas handler and add the VTK handler to make sure IMAS-Python logs are
# correctly displayed in Paraview's logs
imas_logger = logging.getLogger("imas")
if len(imas_logger.handlers) == 1:
    imas_logger.removeHandler(imas_logger.handlers[0])
imas_logger.addHandler(handler)
imas_logger.setLevel(handler.get_level())

# Check IMAS-Python version
if not hasattr(imas, "__version__"):
    raise RuntimeError(
        """
[ERROR] Detected an outdated version of the 'imas' module.

The installed 'imas' package appears to be an incompatible legacy
version of the high-level Python interface of the IMAS Access Layer.

To resolve this, remove / unload this version and re-install using:

    pip install imas-python

or load the appropriate environment module on your system, e.g.

    module load IMAS-Python

More info: https://pypi.org/project/imas-python/
"""
    )
imas_version = imas.__version__
required_imas_python_version = "2.0.0"
if version.parse(imas_version) < version.parse(required_imas_python_version):
    logger.warning(
        f"IMAS-Python version {imas_version} is lower than the recommended version "
        f"{required_imas_python_version}. Some features might not work as expected."
    )
