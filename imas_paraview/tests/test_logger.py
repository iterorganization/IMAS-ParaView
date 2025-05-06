import logging
from unittest.mock import Mock, patch

import pytest

from imas_paraview import VTKHandler


@pytest.fixture
def mock_output_window():
    with patch("imas_paraview.vtkhandler.win") as mock_vtk_win:
        mock_output_window = Mock()
        mock_vtk_win.GetInstance.return_value = mock_output_window
        yield mock_output_window


@pytest.fixture
def logger(request):
    # Ensure each test gets its own logger
    logger_name = f"test_logger_{request.node.name}"
    logger = logging.getLogger(logger_name)
    handler = VTKHandler()
    logger.addHandler(handler)
    return logger


def test_error_log(mock_output_window, logger):
    """Test logging behavior across different levels."""

    logger.setLevel(logging.ERROR)
    logger.error("This is an error message")

    # Ensure other text types were not called
    mock_output_window.DisplayWarningText.assert_not_called()
    mock_output_window.DisplayDebugText.assert_not_called()
    mock_output_window.DisplayText.assert_not_called()

    # Ensure the error text was called
    mock_output_window.DisplayErrorText.assert_called_once()


def test_warning_log(mock_output_window, logger):
    """Test logging behavior for WARNING level."""

    logger.setLevel(logging.WARNING)
    logger.warning("This is a warning message")

    # Ensure other text types were not called
    mock_output_window.DisplayErrorText.assert_not_called()
    mock_output_window.DisplayDebugText.assert_not_called()
    mock_output_window.DisplayText.assert_not_called()

    # Ensure the warning text was called
    mock_output_window.DisplayWarningText.assert_called_once()


def test_info_log(mock_output_window, logger):
    """Test logging behavior for INFO level."""

    logger.setLevel(logging.INFO)
    logger.info("This is an info message")

    # Ensure other text types were not called
    mock_output_window.DisplayErrorText.assert_not_called()
    mock_output_window.DisplayWarningText.assert_not_called()
    mock_output_window.DisplayDebugText.assert_not_called()

    # Ensure the general text was called
    mock_output_window.DisplayText.assert_called_once()


def test_debug_log(mock_output_window, logger):
    """Test logging behavior for DEBUG level."""

    logger.setLevel(logging.DEBUG)
    logger.debug("This is a debug message")

    # Ensure other text types were not called
    mock_output_window.DisplayErrorText.assert_not_called()
    mock_output_window.DisplayWarningText.assert_not_called()
    mock_output_window.DisplayText.assert_not_called()

    # Ensure the debug text was called
    mock_output_window.DisplayDebugText.assert_called_once()
