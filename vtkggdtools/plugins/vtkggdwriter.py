"""IMASPy version of the paraview plugin classes."""

import datetime
import getpass
import logging

import imaspy
import imaspy.ids_defs
import numpy as np
from paraview.util.vtkAlgorithm import smdomain, smhint, smproperty, smproxy
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonCore import vtkDataArray
from vtkmodules.vtkCommonDataModel import (
    vtkCellData,
    vtkDataObject,
    vtkPointSet,
    vtkUnstructuredGrid,
)

from vtkggdtools._version import get_versions
from vtkggdtools.io import write_geom, write_ps
from vtkggdtools.io.representables import GridSubsetRepresentable
from vtkggdtools.util import FauxIndexMap, create_first_grid

logger = logging.getLogger("vtkggdtools")


_ggd_types_xml = "".join(
    f'<Entry text="{member.name}" value="{member.index}" />'
    for member in imaspy.identifiers.ggd_identifier
)
_ggd_space_types_xml = "".join(
    f'<Entry text="{member.name}" value="{member.index}" />'
    for member in imaspy.identifiers.ggd_space_identifier
)


@smproxy.filter(label="GGD Writer")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkPointSet"], composite_data_supported=True)
class IMASPyGGDWriter(VTKPythonAlgorithmBase):
    """GGD Writer based on IMASPy"""

    def __init__(self):
        super().__init__(
            nInputPorts=1,
            nOutputPorts=1,
            inputType="vtkPointSet",
            outputType="vtkUnstructuredGrid",
        )

        # Properties
        self._uri = ""
        self._idsname = ""
        self._occurrence = 0

        self._grid_ggd_type = 0
        self._grid_ggd_space_type = 0

        # Cache
        self._dbentry = None

    @smproperty.stringvector(name="IMAS URI", default_values="")
    def SetURI(self, uri):
        if uri != self._uri:
            self._uri = uri
            if self._dbentry is not None:
                self._dbentry.close()
                self._dbentry = None
            self.Modified()

    @smproperty.stringvector(name="IDS", default_values="")
    def SetIDS(self, idsname):
        if idsname != self._idsname:
            self._idsname = idsname
            self.Modified()

    @smproperty.intvector(name="Occurrence", default_values=0)
    def SetOccurrence(self, occurrence):
        if occurrence != self._occurrence:
            self._occurrence = occurrence
            self.Modified()

    @smproperty.xml(
        f"""
        <IntVectorProperty command="SetGridGGDType"
                           default_values="4"
                           name="GridGGDType"
                           number_of_elements="1">
            <EnumerationDomain name="enum">
                {_ggd_types_xml}
            </EnumerationDomain>
            <Documentation>This property determines the grid ggd type.</Documentation>
        </IntVectorProperty>"""
    )
    def SetGridGGDType(self, val):
        if val != self._grid_ggd_type:
            self._grid_ggd_type = val
            self.Modified()

    @smproperty.xml(
        f"""
        <IntVectorProperty command="SetGridGGDSpaceType"
                           default_values="1"
                           name="GridGGDSpaceType"
                           number_of_elements="1">
            <EnumerationDomain name="enum">
                {_ggd_space_types_xml}
            </EnumerationDomain>
            <Documentation>This property determines the space type for the grid ggd
            geometry.</Documentation>
        </IntVectorProperty>"""
    )
    def SetGridGGDSpaceType(self, val):
        if val != self._grid_ggd_space_type:
            self._grid_ggd_space_type = val
            self.Modified()

    def FillInputPortInformation(self, port, info):
        info.Set(vtkDataObject.DATA_TYPE_NAME(), "vtkPointSet")
        return 1

    def FillOutputPortInformation(self, port, info):
        info.Set(vtkDataObject.DATA_TYPE_NAME(), "vtkUnstructuredGrid")
        return 1

    def RequestDataObject(self, request, inInfo, outInfo):
        output = vtkUnstructuredGrid()
        outInfo.GetInformationObject(0).Set(vtkDataObject.DATA_OBJECT(), output)
        outInfo.GetInformationObject(0).Set(
            vtkDataObject.DATA_EXTENT_TYPE(), output.GetExtentType()
        )
        return 1

    def RequestData(self, request, inInfo, outInfo):
        input0 = vtkPointSet.GetData(inInfo[0])
        output = vtkUnstructuredGrid.GetData(outInfo)
        # FIXME: the following line results in an AttributeError: 'NoneType' object has
        # no attribute 'GetCellData'. This was already the case before migrating to
        # IMASPy. To be investigated when working on the VTK -> GGD conversion logic.
        cell_data: vtkCellData = input0.GetCellData()

        if self._dbentry is None:
            self._dbentry = imaspy.DBEntry(self._uri, "a")

        # Fill IDS metadata
        ids_obj = self._dbentry.factory.new(self._idsname)
        ids_obj.ids_properties.homogeneous_time = 1
        ids_obj.ids_properties.provider = getpass.getuser()
        ids_obj.ids_properties.creation_date = str(datetime.datetime.today())

        ids_obj.code.name = "VTKGGDTools"
        ids_obj.code.version = get_versions()["version"]
        ids_obj.code.commit = get_versions()["full-revisionid"]
        ids_obj.code.repository = "https://git.iter.org/scm/imex/ggd-vtk.git"
        ids_obj.code.output_flag = np.array([0])

        # TODO: allow selecting other grids
        _aos_index_values = FauxIndexMap()
        grid_ggd = create_first_grid(ids_obj)

        # FIXME: time information?
        if not len(ids_obj.time):
            ids_obj.time.resize(1)
            ids_obj.time[0] = 0.0

        if logger is not ids_obj and logger is not None:
            logger.debug(f"grid_ggd: {grid_ggd}")
        else:
            logger.info("Failed to find a grid_ggd node.")
            return 0

        # We now have the grid_ggd
        # 1. Populate the space AoS
        # Implicit: nodes, edges, cells, volumes
        num_subsets = cell_data.GetNumberOfArrays() + 4
        grid_ggd.grid_subset.resize(num_subsets)
        grid_ggd.space.resize(1)
        space_idx = 0
        logger.info(f"Populating grid_ggd/space[{space_idx}]")

        space_name = imaspy.identifiers.ggd_space_identifier(
            self._grid_ggd_space_type
        ).name
        description = imaspy.identifiers.ggd_space_identifier[space_name].description

        grid_ggd.space[space_idx].identifier.name = space_name
        grid_ggd.space[space_idx].identifier.index = self._grid_ggd_space_type
        grid_ggd.space[space_idx].identifier.description = description

        # Fill up space with 0,1,2,3 d geometry from the first partition.
        # It is assumed all other partitions share same points.
        logger.info("Populating basic grid_subsets")
        grid_ggd_rep = write_geom.fill_grid_ggd_basic_geometry(
            input0, space_idx, grid_ggd
        )
        output.ShallowCopy(grid_ggd_rep.ugrid)
        subset_rep = GridSubsetRepresentable(
            ugrid=grid_ggd_rep.ugrid, num_subsets=num_subsets
        )

        # 2. Populate the grid_subset AoS
        # Every partitioned dataset goes into its own grid_subset.
        for subset_idx in range(4, num_subsets):
            data_array: vtkDataArray = cell_data.GetArray(subset_idx - 4)
            subset_name = data_array.GetName()
            write_geom.convert_vtk_dataset_to_grid_subset_geometry(
                grid_ggd_rep, subset_rep, space_idx, subset_idx, subset_name, grid_ggd
            )
        write_ps.write_plasma_state(
            self._idsname,
            ids_obj,
            _aos_index_values,
            space_idx,
            grid_ggd_rep,
            subset_rep,
            grid_ggd,
        )
        self._dbentry.put(ids_obj, self._occurrence)
        return 1

    def Write(self):
        self.Modified()
        self.Update()
