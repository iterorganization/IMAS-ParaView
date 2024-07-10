import logging
import sys


def setup_logger():
    formatter = logging.Formatter(
        "%(filename)s:%(lineno)d   %(levelname)s| %(message)s"
    )
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(formatter)

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)

    return logger
