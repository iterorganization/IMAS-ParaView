import logging

from . import _version

__version__ = _version.get_versions()["version"]

logging.basicConfig()  # FIXME: setup proper logging
