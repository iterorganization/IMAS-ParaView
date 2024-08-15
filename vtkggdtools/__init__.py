import logging

from . import _version

__version__ = _version.get_versions()["version"]
logging.basicConfig()

# Create a logger object for your package/module
logger = logging.getLogger(__name__)

# Set the default logging level (this can be overridden later)
logger.setLevel(logging.DEBUG)
