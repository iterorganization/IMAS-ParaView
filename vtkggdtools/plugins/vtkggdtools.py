"""IMASPy version of the paraview plugin classes.
"""

import datetime
import getpass
import logging

import imaspy
import imaspy.ids_defs
import numpy as np
from identifiers.ggd_identifier import ggd_identifier
from identifiers.ggd_space_identifier import ggd_space_identifier
from imaspy.exception import UnknownDDVersion
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
from vtkggdtools.imas_uri import uri_from_path, uri_from_pulse_run
from vtkggdtools.io import read_bezier, read_geom, read_ps, write_geom, write_ps
from vtkggdtools.io.representables import GridSubsetRepresentable
from vtkggdtools.paraview_support.servermanager_tools import (
    doublevector,
    enumeration,
    genericdecorator,
    intvector,
    propertygroup,
    stringlistdomain,
    stringvector,
)
from vtkggdtools.util import FauxIndexMap, create_first_grid, get_grid_ggd

logger = logging.getLogger("vtkggdtools")

BACKENDS = {
    "MDSplus": imaspy.ids_defs.MDSPLUS_BACKEND,
    "HDF5": imaspy.ids_defs.HDF5_BACKEND,
    "ASCII": imaspy.ids_defs.ASCII_BACKEND,
}
"""Mapping of UI labels for each backend and their ID, used for the Backend dropdown."""
DEFAULT_BACKEND = imaspy.ids_defs.MDSPLUS_BACKEND
"""Default backend selected in the UI."""


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
        # URI properties
        self._uri_selection_mode = 1
        self._uri_input = ""
        self._uri_path = ""
        self._uri_backend = DEFAULT_BACKEND
        self._uri_pulse = 0
        self._uri_run = 0
        self._uri_database = "ITER"
        self._uri_user = "public"
        self._uri_version = "3"
        # IDS properties
        self._ids_and_occurrence = ""
        # Bezier interpolation properties
        self._n_plane = 0
        self._phi_start = 0
        self._phi_end = 0

        # URI is calculated from the possible inputs
        self._uri = ""
        self._uri_error = ""

        # Data caches
        self._dbentry = None
        self._ids_list = []
        self._ids = None

        # Load ggd_idx from paraview UI
        self._time_steps = []

    def _update_property(self, name, value, callback=None):
        """Convenience method to update a property when value changed."""
        if getattr(self, name) != value:
            setattr(self, name, value)
            if callback is not None:
                callback()
            self.Modified()

    def _calculate_uri(self) -> None:
        """Determine the IMAS URI from user input."""
        if self._uri_selection_mode == 1:
            # Use URI provided by user:
            uri = self._uri_input
        elif self._uri_selection_mode == 2:
            # Determine URI from selected file
            uri = uri_from_path(self._uri_path)
        elif self._uri_selection_mode == 3:
            uri = uri_from_pulse_run(
                self._uri_backend,
                self._uri_database,
                self._uri_pulse,
                self._uri_run,
                self._uri_user,
                self._uri_version,
            )
        else:
            raise ValueError(f"Invalid value for {self._uri_selection_mode=}")

        # Ensure uri is not None:
        uri = uri or ""

        if uri != self._uri:
            self._uri = uri
            if self._dbentry is not None:
                self._dbentry.close()
                self._ids = self._dbentry = None
            self._uri_error = ""
            if self._uri:
                # Try to open the DBEntry
                try:
                    self._dbentry = imaspy.DBEntry(self._uri, "r")
                except Exception as exc:
                    self._uri_error = str(exc)
            self._update_ids_list()
            self.Modified()

    def _update_ids_list(self) -> None:
        """Update the list of available IDSs in the selected Data Entry."""
        self._ids_list = []
        if self._dbentry is not None:
            for ids_name in read_ps.SUPPORTED_IDS_NAMES:
                for occurrence in self._dbentry.list_all_occurrences(ids_name):
                    val = ids_name if occurrence == 0 else f"{ids_name}/{occurrence}"
                    self._ids_list.append(val)

    # Note on property names:
    # - Properties are sorted in the GUI by name of the function in paraview < 5.13
    # - In paraview 5.13 and newer, this sorting is based on declaration order (see
    #   https://gitlab.kitware.com/paraview/paraview/-/commit/7ea6373ff67c59a7f43b63697aaf6f1ddee51cab)
    # - To properly sort the properties in paraview <= 5.12 we prefix each property with
    #   "P<groupID><propertyID>", e.g. P01_SetURI

    # Properties for setting the URI
    ####################################################################################

    @intvector(name="URISelection", label="", default_values="1")
    @enumeration("enum", {"Enter URI": 1, "Select file": 2, "Enter pulse, run, ..": 3})
    def P00_URISelectionMode(self, mode):
        """Select mode for locating the data entry: manually enter the URI, select a
        local HDF5/MDSplus file, or enter 'legacy' parameters."""
        self._update_property("_uri_selection_mode", mode, self._calculate_uri)

    @stringvector(name="Enter URI", default_values="")
    @genericdecorator(mode="visibility", property="URISelection", value="1")
    def P01_SetURI(self, uri):
        """The IMAS URI for the Data Entry to load."""
        self._update_property("_uri_input", uri, self._calculate_uri)

    @stringvector(name="URI from file", default_values="")
    @smdomain.filelist()
    @smhint.filechooser(
        extensions=["h5", "datafile"],
        file_description="IMAS HDF5/MDSplus backend files",
    )
    @genericdecorator(mode="visibility", property="URISelection", value="2")
    def P02_SetURIFromFileName(self, file):
        """Select a local file from the HDF5 or MDSplus backend."""
        self._update_property("_uri_path", file, self._calculate_uri)

    @intvector(name="Backend", default_values=DEFAULT_BACKEND)
    @enumeration("backend", BACKENDS)
    @genericdecorator(mode="visibility", property="URISelection", value="3")
    def P03_SetBackend(self, backend):
        """Select the IMAS backend that stores the data."""
        self._update_property("_uri_backend", backend, self._calculate_uri)

    @stringvector(name="Database", default_values="")
    @genericdecorator(mode="visibility", property="URISelection", value="3")
    def P04_SetDatabase(self, database):
        """Enter the Database name where the data is stored."""
        self._update_property("_uri_database", database, self._calculate_uri)

    @intvector(name="Pulse", default_values="0")
    @genericdecorator(mode="visibility", property="URISelection", value="3")
    def P05_SetPulse(self, pulse):
        """Pulse number."""
        self._update_property("_uri_pulse", pulse, self._calculate_uri)

    @intvector(name="Run", default_values="0")
    @genericdecorator(mode="visibility", property="URISelection", value="3")
    def P06_SetRun(self, run):
        """Run number."""
        self._update_property("_uri_run", run, self._calculate_uri)

    @stringvector(name="User", default_values=getpass.getuser())
    @genericdecorator(mode="visibility", property="URISelection", value="3")
    def P07_SetUser(self, user):
        """User name for the IMAS database."""
        self._update_property("_uri_user", user, self._calculate_uri)

    @stringvector(name="Version", default_values="3")
    @genericdecorator(mode="visibility", property="URISelection", value="3")
    def P08_SetVersion(self, version):
        """Major version of the DD used to store the data."""
        self._update_property("_uri_version", version, self._calculate_uri)

    @stringvector(name="Status", information_only=1, panel_visibility="default")
    @smhint.xml('<Widget type="one_liner_wrapped" />')
    def P09_GetURIStatus(self):
        """Status information related to the URI."""
        if not self._uri:
            return "No URI selected, press Apply to set URI."
        if self._uri_error:
            return f'Could not open URI "{self._uri}":\n{self._uri_error}'
        return (
            f'Successfully opened URI "{self._uri}".\n'
            f"This Data Entry contains {len(self._ids_list)} supported IDSs."
        )

    # Properties for setting the IDS name and occurrence
    ####################################################################################

    def _clear_ids(self):
        """Helper callback to clear any loaded IDS when idsname/occurrence change."""
        self._ids = None

    @stringvector(name="IDSAndOccurrence", label="IDS/Occurrence")
    @stringlistdomain("IDSList", name="ids_list", none_string="&lt;Select IDS&gt;")
    def P10_SetIDSAndOccurrence(self, value):
        """Select the IDS/occurrence to load.

        The dropdown is disabled when no Data Entry has been loaded (press Apply after
        selecting a URI to load the Data Entry), or when the loaded Data Entry contains
        no IDSs supported by this plugin.
        """
        if value not in self._ids_list:
            value = ""
        self._update_property("_ids_and_occurrence", value, self._clear_ids)

    @stringvector(name="IDSList", information_only=1)
    def P11_GetIDSList(self):
        """Return a list of IDSs with data inside the selected Data Entry."""
        return self._ids_list

    # Properties for handling time steps
    ####################################################################################

    @smproperty.doublevector(
        name="TimestepValues", information_only="1", si_class="vtkSITimeStepsProperty"
    )
    def GetTimestepValues(self):
        return self._time_steps

    # Properties for Bezier interpolation
    ####################################################################################

    @intvector(name="N plane", default_values=0)
    def P20_SetNPlane(self, val):
        self._update_property("_n_plane", val)

    @doublevector(name="Phi range", default_values=[0, 0])
    @smdomain.doublerange(min=0, max=360.0)
    def P21_SetPhiRange(self, val, val2):
        self._update_property("_phi_start", val)
        self._update_property("_phi_end", val2)

    # Property groups: sorted by name and must be alphabetically after the properties
    ####################################################################################

    @propertygroup(
        "Data entry URI",
        [
            "URISelection",
            "Enter URI",
            "URI from file",
            "Backend",
            "Database",
            "Pulse",
            "Run",
            "User",
            "Version",
            "Status",
        ],
    )
    def PG0_DataEntryGroup(self):
        """Dummy function to define a PropertyGroup."""

    @propertygroup("Select IDS", ["IDSAndOccurrence", "IDSList"])
    def PG1_IDSGroup(self):
        """Dummy function to define a PropertyGroup."""

    @propertygroup("Bezier interpolation settings", ["N plane", "Phi range"])
    def PG2_BezierGroup(self):
        """Dummy function to define a PropertyGroup."""

    # Implement VTK algorithm
    ####################################################################################

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
            idsname, _, occurrence = self._ids_and_occurrence.partition("/")
            occurrence = int(occurrence or 0)
            logger.info("Loading IDS %s/%d ...", idsname, occurrence)
            # TODO: Add option to turn off lazy loading
            lazy = True
            try:
                ids = self._dbentry.get(
                    idsname, occurrence, autoconvert=False, lazy=lazy
                )
            except UnknownDDVersion:
                # Apparently IMASPy doesn't know the DD version that this IDS was
                # written with. Use the default DD version instead:
                ids = self._dbentry.get(idsname, occurrence, lazy=lazy)
            self._ids = ids
            logger.info("Done loading IDS.")

    def RequestInformation(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence:
            return 1

        # Load IDS and available time steps
        self._ensure_ids()
        self._time_steps = self._ids.time

        # Pass time steps to Paraview
        executive = self.GetExecutive()
        outInfo = outInfo.GetInformationObject(0)
        outInfo.Remove(executive.TIME_STEPS())
        outInfo.Remove(executive.TIME_RANGE())
        for time_step in self._time_steps:
            outInfo.Append(executive.TIME_STEPS(), time_step)

        outInfo.Append(executive.TIME_RANGE(), self._time_steps[0])
        outInfo.Append(executive.TIME_RANGE(), self._time_steps[-1])
        return 1

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence:
            return 1

        # Retrieve time step from time selection in Paraview UI
        executive = self.GetExecutive()
        outInfo = outInfo.GetInformationObject(0)
        time_step = outInfo.Get(executive.UPDATE_TIME_STEP())
        time_step_idx = np.where(self._time_steps == time_step)[0]
        if len(time_step_idx) == 0:
            logger.warning("Selected invalid time step")
            return 1
        else:
            time_step_idx = time_step_idx[0]
            logger.debug(f"Selected time step: {self._time_steps[time_step_idx]}")

        # TODO: allow selecting other grids
        _aos_index_values = FauxIndexMap()
        grid_ggd = get_grid_ggd(self._ids, time_step_idx)

        # We now have the grid_ggd
        # Check if we have anything to read:
        if len(grid_ggd.grid_subset) < 1 or len(grid_ggd.space) < 1:
            logger.info("No points to read from grid_ggd")
            return 1

        output = vtkPartitionedDataSetCollection.GetData(outInfo)
        num_subsets = len(grid_ggd.grid_subset)
        points = vtkPoints()
        space_idx = 0
        idsname = self._ids.metadata.name

        read_geom.fill_vtk_points(grid_ggd, space_idx, points, idsname)
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
                        idsname,
                        self._ids,
                        _aos_index_values,
                        self._n_plane,
                        self._phi_start,
                        self._phi_end,
                    )
                    output.SetPartition(0, 0, ugrid)
                    child = assembly.AddNode(idsname, 0)
                    assembly.AddDataSetIndex(child, 0)
                    output.GetMetaData(0).Set(vtkCompositeDataSet.NAME(), idsname)
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

        def fill_grid_and_plasma_state(ps_reader, subset_idx, partition):
            subset = None if subset_idx < 0 else grid_ggd.grid_subset[subset_idx]
            ugrid = read_geom.convert_grid_subset_geometry_to_unstructured_grid(
                grid_ggd, subset_idx, points
            )
            ps_reader.read_plasma_state(subset_idx, ugrid)
            output.SetPartition(partition, 0, ugrid)
            label = str(subset.identifier.name) if subset else idsname
            child = assembly.AddNode(label.replace(" ", "_"), 0)
            assembly.AddDataSetIndex(child, partition)
            output.GetMetaData(partition).Set(vtkCompositeDataSet.NAME(), label)

        # Regular grid reading
        ps_reader = read_ps.PlasmaStateReader(self._ids, time_step_idx)
        if num_subsets <= 1:
            logger.info("No subsets to read from grid_ggd")
            output.SetNumberOfPartitionedDataSets(1)
            fill_grid_and_plasma_state(ps_reader, -1, 0)
        elif idsname == "wall":
            # FIXME: what if num_subsets is 2 or 3?
            output.SetNumberOfPartitionedDataSets(num_subsets - 3)
            fill_grid_and_plasma_state(ps_reader, -1, 0)
            for subset_idx in range(4, num_subsets):
                fill_grid_and_plasma_state(ps_reader, subset_idx, subset_idx - 3)
                self.UpdateProgress(self.GetProgress() + 1 / num_subsets)
        else:
            output.SetNumberOfPartitionedDataSets(num_subsets)
            for subset_idx in range(num_subsets):
                fill_grid_and_plasma_state(ps_reader, subset_idx, subset_idx)
                self.UpdateProgress(self.GetProgress() + 1 / num_subsets)

        return 1


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
