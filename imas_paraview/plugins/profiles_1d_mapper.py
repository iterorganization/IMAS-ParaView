import logging

import numpy as np
from paraview.util.vtkAlgorithm import smdomain, smhint, smproperty, smproxy
from vtk.util.numpy_support import vtk_to_numpy
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonDataModel import vtkPartitionedDataSetCollection, vtkTable

from imas_paraview.paraview_support.servermanager_tools import (
    arrayselectiondomain,
    arrayselectionstringvector,
)
from imas_paraview.progress import Progress

logger = logging.getLogger("imas_paraview")


@smproxy.filter(label="1D Profiles Mapper")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
@smproperty.input(name="1D Profile", port_index=1)
@smdomain.datatype(dataTypes=["vtkTable"], composite_data_supported=False)
@smproperty.input(name="Psi Grid", port_index=0)
@smdomain.datatype(
    dataTypes=["vtkPartitionedDataSetCollection"], composite_data_supported=True
)
class Profiles1DMapper(VTKPythonAlgorithmBase):
    """Filter that maps 1d profiles onto a 2D poloidal flux grid. This filter takes
    2 inputs:
    - A poloidal flux (Psi) from a GGD. i.e. from GGD Reader/JOREK Reader (input 0)
    - A 1D profile from the 1DProfilesReader. (input 1)
    The 1DProfilesReader must output at least the Psi grid coordinate, and can output
    any number of 1D profiles."""

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self,
            nInputPorts=2,
            nOutputPorts=1,
            outputType="vtkPartitionedDataSetCollection",
        )
        self._selected = []
        self._selectable = []

    def FillInputPortInformation(self, port, info):
        if port == 1:
            info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkTable")
        else:
            info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkPartitionedDataSetCollection")
        return 1

    def RequestData(self, request, inInfoVec, outInfo):
        input0 = dsa.WrapDataObject(
            vtkPartitionedDataSetCollection.GetData(inInfoVec[0], 0)
        )
        input1 = dsa.WrapDataObject(vtkTable.GetData(inInfoVec[1], 0))

        self._selectable = input1.RowData.keys()
        psi = input0.PointData["Psi [Wb]"]

        if isinstance(psi, dsa.VTKNoneArray):
            logger.warning(
                "The GGD Reader should output a poloidal flux GGD. Please select 'Psi' "
                "in the attribute selector window."
            )
            return 1
        psi_grid = psi.GetArrays()[0]
        psi_profiles = input1.RowData["Grid Psi"]

        if isinstance(psi_profiles, dsa.VTKNoneArray) and self._selected:
            logger.warning(
                "The 1DProfilesReader should output a poloidal flux grid. Please"
                " select 'Grid Psi' in the attribute selector window of the "
                "1DProfilesReader."
            )
            return 1

        output = dsa.WrapDataObject(vtkPartitionedDataSetCollection.GetData(outInfo))
        output.ShallowCopy(input0.VTKObject)

        progress = Progress(self.UpdateProgress)

        for profile_name in self._selected:
            profile = input1.RowData[profile_name]

            xp = vtk_to_numpy(psi_profiles)
            fp = vtk_to_numpy(profile)

            # numpy's interp requires xp to be strictly increasing, and the DD does not
            # enforce a specific ordering. So we invert the array if it is decreasing
            if xp[0] > xp[-1]:
                xp = xp[::-1]
                fp = fp[::-1]

            psi_grid_values = vtk_to_numpy(psi_grid)
            resample = np.interp(psi_grid_values, xp, fp, left=np.nan, right=np.nan)

            progress.increment(1 / len(self._selectable))
            output.PointData.append(resample, f"{profile_name} (resampled)")
        return 1

    @arrayselectiondomain(
        property_name="Profile1DArray",
        name="Profile1DSelector",
        label="Select 1D Profiles",
    )
    def SetProfile1D(self, array, status):
        """Select all or a subset of available 1d_profiles to load."""
        # Add an array to selected list
        if status == 1 and array not in self._selected:
            self._selected.append(array)
            self.Modified()

        # Remove an array from selected list
        if status == 0 and array in self._selected:
            self._selected.remove(array)
            self.Modified()

    @arrayselectionstringvector(
        property_name="Profile1DArray", attribute_name="Profile1D"
    )
    def _Profile1DArraySelector(self):
        pass

    def GetNumberOfProfile1DArrays(self):
        return len(self._selectable)

    def GetProfile1DArrayName(self, idx) -> str:
        return self._selectable[idx]

    def GetProfile1DArrayStatus(self, *args):
        return 1
