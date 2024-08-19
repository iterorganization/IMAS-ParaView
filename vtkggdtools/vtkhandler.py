# SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
# SPDX-License-Identifier: BSD-3-Clause

import logging

from vtkmodules.vtkCommonCore import vtkLogger
from vtkmodules.vtkCommonCore import vtkOutputWindow as win


class VTKHandler(logging.Handler):
    """
    A custom logging handler for VTK, intended to work independently from ParaView,
    adapted from paraview.detail.loghandler.py.

    This handler formats and logs messages using VTK's logging system and displays
    them in the VTK output window. It maps Python's standard logging levels to VTK's
    verbosity levels.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def emit(self, record):
        """Formats and logs a record using VTK's logging system and displays it in the
        VTK output window.

        Args:
            record: The log record to be processed.

        Raises:
            Exception: Any exceptions raised during the logging process will be
            handled by the `handleError` method.
        """
        try:
            msg = self.format(record)
            lvl = self.get_vtk_level(record.levelno)
            vtkLogger.Log(lvl, record.filename, record.lineno, msg)

            outputWindow = win.GetInstance()
            if outputWindow:
                # do not duplicate on standard output
                prevMode = outputWindow.GetDisplayMode()
                outputWindow.SetDisplayModeToNever()

                if lvl == vtkLogger.VERBOSITY_ERROR:
                    lvlText = "ERR: "
                    fullMsg = f"{record.filename}:{record.lineno} {lvlText}{msg}\n"
                    outputWindow.DisplayErrorText(fullMsg)
                elif lvl == vtkLogger.VERBOSITY_WARNING:
                    lvlText = "WARN: "
                    fullMsg = f"{record.filename}:{record.lineno} {lvlText}{msg}\n"
                    outputWindow.DisplayWarningText(fullMsg)
                elif lvl == vtkLogger.VERBOSITY_INFO:
                    lvlText = "INFO: "
                    fullMsg = f"{record.filename}:{record.lineno} {lvlText}{msg}\n"
                    outputWindow.DisplayText(fullMsg)
                elif lvl == vtkLogger.VERBOSITY_TRACE:
                    lvlText = "DEBUG: "
                    fullMsg = f"{record.filename}:{record.lineno} {lvlText}{msg}\n"
                    outputWindow.DisplayDebugText(fullMsg)
                else:
                    fullMsg = f"{record.filename}:{record.lineno} {msg}\n"
                    outputWindow.DisplayText(fullMsg)

                outputWindow.SetDisplayMode(prevMode)

        except Exception:
            self.handleError(record)

    def get_vtk_level(self, level):
        """Converts a Python logging level to a VTK verbosity level.

        Args:
            level: The Python logging level (e.g., logging.ERROR, logging.WARNING).

        Returns:
            The corresponding VTK verbosity level.
        """
        if level >= logging.ERROR:
            return vtkLogger.VERBOSITY_ERROR
        elif level >= logging.WARNING:
            return vtkLogger.VERBOSITY_WARNING
        elif level >= logging.INFO:
            return vtkLogger.VERBOSITY_INFO
        elif level >= logging.DEBUG:
            return vtkLogger.VERBOSITY_TRACE
        else:
            return vtkLogger.VERBOSITY_MAX

    def get_level(self, vtklevel=None):
        """Converts a VTK verbosity level to a Python logging level. If VTK verbosity
        level is provided, it defaults to using the current verbosity cutoff defined in
        the VTK logger.

        Args:
            vtklevel: The VTK verbosity level to convert.

        Returns:
            The corresponding Python logging level.
        """
        from vtkmodules.vtkCommonCore import vtkLogger

        if vtklevel is None:
            vtk_level = vtkLogger.GetCurrentVerbosityCutoff()

        if vtk_level == vtkLogger.VERBOSITY_ERROR:
            return logging.ERROR
        elif vtk_level == vtkLogger.VERBOSITY_WARNING:
            return logging.WARNING
        elif vtk_level == vtkLogger.VERBOSITY_INFO:
            return logging.INFO
        else:
            return logging.DEBUG
