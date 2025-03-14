import numpy as np
import vtk
from paraview.util.vtkAlgorithm import smdomain, smproperty, smproxy
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonDataModel import vtkPartitionedDataSetCollection, vtkTable

from vtkggdtools.paraview_support.servermanager_tools import (
    arrayselectiondomain,
    arrayselectionstringvector,
)
from vtkggdtools.progress import Progress


@smproxy.filter(label="test_filter")
@smproperty.input(name="1D Profile", port_index=1)
@smdomain.datatype(dataTypes=["vtkTable"], composite_data_supported=False)
@smproperty.input(name="Psi Grid", port_index=0)
@smdomain.datatype(
    dataTypes=["vtkPartitionedDataSetCollection"], composite_data_supported=True
)
class ExampleTwoInputFilter(VTKPythonAlgorithmBase):
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

        psi_grid = input0.PointData["Psi [Wb]"].GetArrays()[0]
        psi_profiles = input1.RowData["Grid Psi"]
        self._selectable = input1.RowData.keys()

        output = dsa.WrapDataObject(vtkPartitionedDataSetCollection.GetData(outInfo))
        output.ShallowCopy(input0.VTKObject)

        # Create progress object to advance Paraview progress bar
        progress = Progress(self.UpdateProgress)
        for profile_name in self._selected:
            profile = input1.RowData[profile_name]

            piecewise_function = vtk.vtkPiecewiseFunction()
            # Interpolate to zero outside the core profiles range
            piecewise_function.ClampingOff()
            for i in range(psi_profiles.GetNumberOfTuples()):
                piecewise_function.AddPoint(
                    psi_profiles.GetTuple1(i), profile.GetTuple1(i)
                )

            # Perform resampling based on the psi values
            resample = np.zeros(psi_grid.GetSize())
            for i in range(psi_grid.GetSize()):
                psi_value = psi_grid.GetValue(i)
                resample[i] = piecewise_function.GetValue(psi_value)

            progress.increment(1 / len(self._selectable))
            output.PointData.append(resample, f"{profile_name} (resampled)")
        return 1

    @arrayselectiondomain(
        property_name="AttributeArray",
        name="AttributeArraySelector",
        label="Select attribute Arrays",
    )
    def P12_SetAttributeArray(self, array, status):
        """Select all or a subset of available GGD arrays to load."""
        # Add an array to selected list
        if status == 1 and array not in self._selected:
            self._selected.append(array)
            self.Modified()

        # Remove an array from selected list
        if status == 0 and array in self._selected:
            self._selected.remove(array)
            self.Modified()

    @arrayselectionstringvector(
        property_name="AttributeArray", attribute_name="Attribute"
    )
    def _AttributeArraySelector(self):
        pass

    def GetNumberOfAttributeArrays(self):
        return len(self._selectable)

    def GetAttributeArrayName(self, idx) -> str:
        return self._selectable[idx]

    def GetAttributeArrayStatus(self, *args):
        return 1
