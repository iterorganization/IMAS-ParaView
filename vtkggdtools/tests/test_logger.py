import logging

import pytest

from vtkggdtools import VTKHandler  # Replace with the correct import path


@pytest.fixture
def logger():
    logger = logging.getLogger("test_logger")
    handler = VTKHandler()
    logger.addHandler(handler)
    return logger


@pytest.mark.parametrize(
    "level, expected_log_count",
    [
        (logging.DEBUG, 1),  # Expect 1 log record for DEBUG level
        (logging.INFO, 0),  # Expect 0 log records for INFO level
        (logging.WARNING, 0),  # Expect 0 log records for WARNING level
        (logging.ERROR, 0),  # Expect 0 log records for ERROR level
    ],
)
def test_debug_log(logger, caplog, level, expected_log_count):
    """Test logging behavior across different levels."""
    logger.setLevel(level)
    logger.debug("This is a debug message")
    assert len(caplog.records) == expected_log_count
    if expected_log_count > 0:
        assert caplog.records[0].levelname == "DEBUG"
        assert caplog.records[0].message == "This is a debug message"


@pytest.mark.parametrize(
    "level, expected_log_count",
    [
        (logging.DEBUG, 1),  # Expect 1 log record for DEBUG level
        (logging.INFO, 1),  # Expect 1 log records for INFO level
        (logging.WARNING, 0),  # Expect 0 log records for WARNING level
        (logging.ERROR, 0),  # Expect 0 log records for ERROR level
    ],
)
def test_info_log(logger, caplog, level, expected_log_count):
    """Test logging behavior across different levels."""
    logger.setLevel(level)
    logger.info("This is an info message")
    assert len(caplog.records) == expected_log_count
    if expected_log_count > 0:
        assert caplog.records[0].levelname == "INFO"
        assert caplog.records[0].message == "This is an info message"


@pytest.mark.parametrize(
    "level, expected_log_count",
    [
        (logging.DEBUG, 1),  # Expect 1 log record for DEBUG level
        (logging.INFO, 1),  # Expect 1 log records for INFO level
        (logging.WARNING, 1),  # Expect 1 log records for WARNING level
        (logging.ERROR, 0),  # Expect 0 log records for ERROR level
    ],
)
def test_warning_log(logger, caplog, level, expected_log_count):
    """Test logging behavior across different levels."""
    logger.setLevel(level)
    logger.warning("This is a warning message")
    assert len(caplog.records) == expected_log_count
    if expected_log_count > 0:
        assert caplog.records[0].levelname == "WARNING"
        assert caplog.records[0].message == "This is a warning message"


@pytest.mark.parametrize(
    "level, expected_log_count",
    [
        (logging.DEBUG, 1),  # Expect 1 log record for DEBUG level
        (logging.INFO, 1),  # Expect 1 log records for INFO level
        (logging.WARNING, 1),  # Expect 1 log records for WARNING level
        (logging.ERROR, 1),  # Expect 1 log records for ERROR level
    ],
)
def test_error_log(logger, caplog, level, expected_log_count):
    """Test logging behavior across different levels."""
    logger.setLevel(level)
    logger.error("This is an error message")
    assert len(caplog.records) == expected_log_count
    if expected_log_count > 0:
        assert caplog.records[0].levelname == "ERROR"
        assert caplog.records[0].message == "This is an error message"
