import logging

from vtkmodules.vtkCommonCore import vtkLogger
from vtkmodules.vtkCommonCore import vtkOutputWindow as win


class VTKHandler(logging.Handler):
    def __init__(self, *args, **kwargs):
        super(VTKHandler, self).__init__(*args, **kwargs)

    def emit(self, record):
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
                else:
                    fullMsg = f"{record.filename}:{record.lineno} {msg}\n"
                    outputWindow.DisplayText(fullMsg)

                outputWindow.SetDisplayMode(prevMode)

        except Exception:
            self.handleError(record)

    def get_vtk_level(self, level):

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
