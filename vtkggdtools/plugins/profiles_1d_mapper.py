import numpy as np
import vtk
from paraview.util.vtkAlgorithm import smdomain, smproperty, smproxy
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonDataModel import (
    vtkPartitionedDataSet,
    vtkPartitionedDataSetCollection,
    vtkPolyData,
    vtkTable,
)


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
        j_profiles = input1.RowData["J_total"]

        piecewise_function = vtk.vtkPiecewiseFunction()
        piecewise_function.ClampingOff()  # Interpolate to zero outside the core profiles range
        for i in range(psi_profiles.GetNumberOfTuples()):
            piecewise_function.AddPoint(
                psi_profiles.GetTuple1(i), j_profiles.GetTuple1(i)
            )

        # Perform resampling based on the psi values
        resample = np.zeros(psi_grid.GetSize())
        for i in range(psi_grid.GetSize()):
            psi_value = psi_grid.GetValue(i)
            resample[i] = piecewise_function.GetValue(psi_value)

        output = dsa.WrapDataObject(vtkPartitionedDataSetCollection.GetData(outInfo))
        output.ShallowCopy(input0.VTKObject)
        output.PointData.append(resample, "J_total (resampled)")
        return 1
