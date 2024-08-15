import logging

from vtkggdtools.vtkhandler import VTKHandler

from . import _version

__version__ = _version.get_versions()["version"]

# Create a logger object for your package/module
logger = logging.getLogger(__name__)

# Configure the logger to use VTKHandler and set the appropriate level
if not hasattr(logger, "_paraview_initialized"):
    logger._paraview_initialized = True

    # Add the VTKHandler from vtkggdtools
    handler = VTKHandler()
    logger.addHandler(handler)
