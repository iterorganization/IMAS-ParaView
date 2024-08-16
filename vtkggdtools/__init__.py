import logging

from vtkggdtools.vtkhandler import VTKHandler

from . import _version

__version__ = _version.get_versions()["version"]

# Create a logger object for your package/module
logger = logging.getLogger(__name__)

if not hasattr(logger, "_paraview_initialized"):
    logger._paraview_initialized = True

    handler = VTKHandler()
    handler.set_vtk_level(logging.DEBUG)
    logger.addHandler(handler)
