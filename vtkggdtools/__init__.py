import logging

from vtkggdtools.vtkhandler import VTKHandler

from . import _version

__version__ = _version.get_versions()["version"]

logger = logging.getLogger(__name__)

# Add VTK handler
handler = VTKHandler()
logger.addHandler(handler)
logger.setLevel(handler.get_level())
