import logging

import imas
from packaging import version

from imas_paraview.vtkhandler import VTKHandler

from . import _version

__version__ = _version.get_versions()["version"]

# Setup logger for GGDVTK module
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
imas_version = imas.__version__
required_imas_python_version = "2.0.0"
if version.parse(imas_version) < version.parse(required_imas_python_version):
    logger.warning(
        f"IMAS-Python version {imas_version} is lower than the recommended version "
        f"{required_imas_python_version}. Some features might not work as expected."
    )
