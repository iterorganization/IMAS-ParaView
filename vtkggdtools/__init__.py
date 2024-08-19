import logging

from vtkggdtools.vtkhandler import VTKHandler

from . import _version

__version__ = _version.get_versions()["version"]

# Setup logger for GGDVTK module
logger = logging.getLogger("ggdvtk")
handler = VTKHandler()
logger.addHandler(handler)
logger.setLevel(handler.get_level())

# Remove imaspy handler and add the VTK handler to make sure IMASPy logs are correctly
# displayed in Paraview's logs
imaspy_logger = logging.getLogger("imaspy")
if len(imaspy_logger.handlers) == 1:
    imaspy_logger.removeHandler(imaspy_logger.handlers[0])
imaspy_logger.addHandler(handler)
imaspy_logger.setLevel(handler.get_level())
