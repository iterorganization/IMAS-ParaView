"""IMASPy version of the paraview plugin classes.
"""

import datetime
import getpass

import imaspy
import numpy as np
from identifiers.ggd_identifier import ggd_identifier
from identifiers.ggd_space_identifier import ggd_space_identifier
from imaspy.exception import UnknownDDVersion
from paraview import logger
from paraview.util.vtkAlgorithm import smdomain, smhint, smproperty, smproxy
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonCore import vtkDataArray, vtkPoints
from vtkmodules.vtkCommonDataModel import (
    vtkCellData,
    vtkCompositeDataSet,
    vtkDataAssembly,
    vtkDataObject,
    vtkPartitionedDataSetCollection,
    vtkPointSet,
    vtkUnstructuredGrid,
)

from vtkggdtools._version import get_versions
from vtkggdtools.io import read_bezier, read_geom, read_ps, write_geom, write_ps
from vtkggdtools.io.representables import GridSubsetRepresentable
from vtkggdtools.util import create_first_grid, get_first_grid


@smproxy.source(label="IMASPy GGDReader")
@smhint.xml("""<ShowInMenu category="VTKGGDTools" />""")
class IMASPyGGDReader(VTKPythonAlgorithmBase):
    """GGD Reader based on IMASPy"""

    def __init__(self):
        super().__init__(
            nInputPorts=0,
            nOutputPorts=1,
            outputType="vtkPartitionedDataSetCollection",
        )
        # Properties
        self._uri = ""
        self._idsname = ""
        self._occurrence = 0
        self._n_plane = 0
        self._phi_start = 0
        self._phi_end = 0

        # Data caches
        self._dbentry = None
        self._ids = None

    @smproperty.stringvector(name="IMAS URI", default_values="")
    def SetURI(self, uri):
        if uri != self._uri:
            self._uri = uri
            if self._dbentry is not None:
                self._dbentry.close()
                self._ids = self._dbentry = None
            self.Modified()

    @smproperty.stringvector(name="IDS", default_values="")
    def SetIDS(self, idsname):
        if idsname != self._idsname:
            self._idsname = idsname
            self._ids = None
            self.Modified()

    @smproperty.intvector(name="Occurrence", default_values=0)
    def SetOccurrence(self, occurrence):
        if occurrence != self._occurrence:
            self._occurrence = occurrence
            self._ids = None
            self.Modified()

    @smproperty.intvector(name="N plane", default_values=0)
    def SetNPlane(self, val):
        if val != self._n_plane:
            self._n_plane = val
            self.Modified()

    @smproperty.doublevector(name="Phi range", default_values=[0, 0])
    @smdomain.doublerange(min=0, max=360.0)
    def SetPhiRange(self, val, val2):
        if val != self._phi_start or val2 != self._phi_end:
            self._phi_start = val
            self._phi_end = val2
            self.Modified()

    def FillOutputPortInformation(self, port, info):
        info.Set(vtkDataObject.DATA_TYPE_NAME(), "vtkPartitionedDataSetCollection")
        return 1

    def RequestDataObject(self, request, inInfo, outInfo):
        output = vtkPartitionedDataSetCollection()
        outInfo.GetInformationObject(0).Set(vtkDataObject.DATA_OBJECT(), output)
        outInfo.GetInformationObject(0).Set(
            vtkDataObject.DATA_EXTENT_TYPE(), output.GetExtentType()
        )
        return 1

    def _ensure_ids(self):
        if self._ids is None:
            if self._dbentry is None:
                self._dbentry = imaspy.DBEntry(self._uri, "r")

            # TODO: enable lazy loading
            logger.info("Loading IDS %s/%d ...", self._idsname, self._occurrence)
            lazy = False  # TODO: Test lazy loading before enabling
            try:
                ids = self._dbentry.get(
                    self._idsname, self._occurrence, autoconvert=False, lazy=lazy
                )
            except UnknownDDVersion:
                # Apparently IMASPy doesn't know the DD version that this IDS was
                # written with. Use the default DD version instead:
                ids = self._dbentry.get(self._idsname, self._occurrence, lazy=lazy)
            self._ids = ids
            logger.info("Done loading IDS.")

    def RequestData(self, request, inInfo, outInfo):
        self._ensure_ids()

        # TODO: allow selecting other grids
        _aos_index_values = FauxIndexMap()
        grid_ggd = get_first_grid(self._ids)

        # We now have the grid_ggd
        # Check if we have anything to read:
        if len(grid_ggd.grid_subset) < 1 and len(grid_ggd.space) < 1:
            logger.info("No points to read from grid_ggd")
            return 0

        output = vtkPartitionedDataSetCollection.GetData(outInfo)
        num_subsets = len(grid_ggd.grid_subset)
        points = vtkPoints()
        space_idx = 0

        read_geom.fill_vtk_points(grid_ggd, space_idx, points, self._idsname)
        assembly = vtkDataAssembly()
        output.SetDataAssembly(assembly)

        # Interpolate JOREK Fourier space
        # TODO: figure out if we can put this functionality in a post-processing step?
        if self._n_plane != 0:
            number_of_spaces = len(grid_ggd.space)
            if number_of_spaces > 1 and len(grid_ggd.space[0].coordinates_type) == 2:
                n_period = grid_ggd.space[1].geometry_type.index
                if n_period > 0:  # Fourier space with periodicity (JOREK)
                    logger.info(
                        f"Reading Bezier mesh with Fourier periodiciy {n_period}"
                    )
                    ugrid = read_bezier.convert_grid_subset_to_unstructured_grid(
                        self._idsname,
                        self._ids,
                        _aos_index_values,
                        self._n_plane,
                        self._phi_start,
                        self._phi_end,
                    )
                    output.SetPartition(0, 0, ugrid)
                    child = assembly.AddNode(self._idsname, 0)
                    assembly.AddDataSetIndex(child, 0)
                    output.GetMetaData(0).Set(vtkCompositeDataSet.NAME(), self._idsname)
                else:
                    logger.error(
                        f"Number of planes {self._n_plane} invalid for this "
                        f"{number_of_spaces} number of spaces"
                    )
            else:
                logger.error(
                    f"Number of planes {self._n_plane} invalid for this IDS type."
                    " Try using N = 0"
                )
            return 1

        def fill_grid_and_plasma_state(subset_idx, partition):
            subset = None if subset_idx < 0 else grid_ggd.grid_subset[subset_idx]
            ugrid = read_geom.convert_grid_subset_geometry_to_unstructured_grid(
                grid_ggd, subset_idx, points
            )
            read_ps.read_plasma_state(
                self._idsname, self._ids, _aos_index_values, subset_idx, ugrid
            )
            output.SetPartition(partition, 0, ugrid)
            label = str(subset.identifier.name) if subset else self._idsname
            child = assembly.AddNode(label.replace(" ", "_"), 0)
            assembly.AddDataSetIndex(child, partition)
            output.GetMetaData(partition).Set(vtkCompositeDataSet.NAME(), label)

        # Regular grid reading
        if num_subsets <= 1:
            logger.info("No subsets to read from grid_ggd")
            output.SetNumberOfPartitionedDataSets(1)
            fill_grid_and_plasma_state(-1, 0)
        elif self._idsname == "wall":
            # FIXME: what if num_subsets is 2 or 3?
            output.SetNumberOfPartitionedDataSets(num_subsets - 3)
            fill_grid_and_plasma_state(-1, 0)
            for subset_idx in range(4, num_subsets):
                fill_grid_and_plasma_state(subset_idx, subset_idx - 3)
                self.UpdateProgress(self.GetProgress() + 1 / num_subsets)
        else:
            output.SetNumberOfPartitionedDataSets(num_subsets)
            for subset_idx in range(num_subsets):
                fill_grid_and_plasma_state(subset_idx, subset_idx)
                self.UpdateProgress(self.GetProgress() + 1 / num_subsets)

        return 1


class FauxIndexMap:
    def __getitem__(self, item):
        return 0

    def get(self, name, default=None):
        return 0


_ggd_types_xml = "".join(
    f'<Entry text="{key}" value="{value.get("index")}" />'
    for key, value in ggd_identifier.items()
)
_ggd_space_types_xml = "".join(
    f'<Entry text="{key}" value="{value.get("index")}" />'
    for key, value in ggd_space_identifier.items()
)


@smproxy.filter(label="IMASPy GGDWriter")
@smhint.xml("""<ShowInMenu category="VTKGGDTools" />""")
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

        space_name = list(ggd_space_identifier.keys())[self._grid_ggd_space_type]
        grid_ggd.space[space_idx].identifier.name = space_name
        grid_ggd.space[space_idx].identifier.index = self._grid_ggd_space_type
        grid_ggd.space[space_idx].identifier.description = ggd_space_identifier.get(
            space_name
        ).get("description")

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
