"""IMASPy version of the paraview plugin classes.
"""

import datetime
import getpass
import logging

import imaspy
import imaspy.ids_defs
import numpy as np
from paraview.util.vtkAlgorithm import smdomain, smhint, smproperty, smproxy
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonCore import vtkDataArray, vtkStringArray
from vtkmodules.vtkCommonDataModel import (
    vtkCellData,
    vtkDataObject,
    vtkPartitionedDataSetCollection,
    vtkPointSet,
    vtkUnstructuredGrid,
)

from vtkggdtools._version import get_versions
from vtkggdtools.convert import ggd_to_vtk
from vtkggdtools.imas_uri import uri_from_path, uri_from_pulse_run
from vtkggdtools.io import read_ps, write_geom, write_ps
from vtkggdtools.io.representables import GridSubsetRepresentable
from vtkggdtools.paraview_support.servermanager_tools import (
    arrayselectiondomain,
    arrayselectionstringvector,
    checkbox,
    doublevector,
    enumeration,
    genericdecorator,
    intvector,
    propertygroup,
    stringlistdomain,
    stringvector,
)
from vtkggdtools.progress import Progress
from vtkggdtools.util import FauxIndexMap, create_first_grid

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

        # GGD arrays to load
        self._selected_paths = []
        self._selectable_paths = []
        self._selectable_vector_paths = []
        self._selectable_scalar_paths = []

        self.grid_ggd = None
        self.ps_reader = None

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
            # Use URI provided by user, removing leading/trailing whitespaces
            uri = self._uri_input
            if uri is not None:
                uri = uri.strip()
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
                    self._selectable_paths = []
                    self._ids_list = []
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

    # Properties for setting the IDS name, occurrence and which GGD arrays to load
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

    @stringvector(name="IDSList", information_only=1, si_class="vtkSIDataArrayProperty")
    def P11_GetIDSList(self):
        """Return a list of IDSs with data inside the selected Data Entry."""
        arr = vtkStringArray()
        for val in self._ids_list:
            arr.InsertNextValue(val)
        return arr

    @arrayselectiondomain(
        property_name="GGDArray",
        name="GGDArraySelector",
        label="Select GGD Arrays",
    )
    def P12_SetGGDArray(self, array, status):
        """Select all or a subset of available GGD arrays to load."""
        # Add a GGD array to selected list
        if status == 1 and array not in self._selected_paths:
            self._selected_paths.append(array)
            self.Modified()

        # Remove a GGD array from selected list
        if status == 0 and array in self._selected_paths:
            self._selected_paths.remove(array)
            self.Modified()

    @arrayselectionstringvector(property_name="GGDArray", attribute_name="GGD")
    def _GGDArraySelector(self):
        pass

    def GetNumberOfGGDArrays(self):
        return len(self._selectable_paths)

    def GetGGDArrayName(self, idx):
        return self._name_from_idspath(self._selectable_paths[idx])

    def GetGGDArrayStatus(self, *args):
        return 1

    @checkbox(
        name="LazyLoading",
        label="Preload Data",
        default_values="0",
    )
    def P13_SetLazyLoading(self, val):
        """Turn on to preload the entire IDS beforehand, if this is left off, the data
        is loaded on-demand through lazy loading. It is recommended to leave this off
        if you are only loading a small subset of the data. If you want to load most of
        the data, it is recommended to turn this on."""
        if val == 0:
            status = "enabled"
            self.lazy = True
        else:
            status = "disabled"
            self.lazy = False
        logger.info(f"Lazy Loading is {status}.")
        self.Modified()

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

    @propertygroup(
        "Select IDS", ["IDSAndOccurrence", "IDSList", "GGDArraySelector", "LazyLoading"]
    )
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

    def RequestInformation(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence:
            return 1

        # Load IDS and available time steps
        # TODO: Add support for IDSs with heterogeneous time mode
        self._ensure_ids()
        if (
            self._ids.ids_properties.homogeneous_time
            != imaspy.ids_defs.IDS_TIME_MODE_HOMOGENEOUS
        ):
            logger.error(
                "Only IDSs with homogeneous time mode are currently supported."
            )
            return 1
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
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        # Retrieve the selected time step and GGD arrays
        selected_scalar_paths, selected_vector_paths = self._get_selected_ggd_paths()
        time_idx = self._get_selected_time_step(outInfo)
        if time_idx is None:
            logger.warning("Selected invalid time step")
            return 1

        # Create progress object to advance Paraview progress bar
        progress = Progress(self.UpdateProgress, self.GetProgress)

        # Convert GGD of IDS to VTK format
        output = ggd_to_vtk(
            self._ids,
            time_idx=time_idx,
            scalar_paths=selected_scalar_paths,
            vector_paths=selected_vector_paths,
            n_plane=self._n_plane,
            phi_start=self._phi_start,
            phi_end=self._phi_end,
            outInfo=outInfo,
            progress=progress,
        )
        if output is None:
            logger.warning("Could not convert GGD to VTK.")
        return 1

    def _ensure_ids(self):
        """
        Loads the IDS if not already loaded. Once loaded, initializes plasma state
        reader and populates scalar and vector paths for selection.
        """
        if self._ids is None:
            idsname, _, occurrence = self._ids_and_occurrence.partition("/")
            occurrence = int(occurrence or 0)
            logger.info("Loading IDS %s/%d ...", idsname, occurrence)

            self._ids = self._dbentry.get(
                idsname,
                occurrence,
                autoconvert=False,
                lazy=self.lazy,
                ignore_unknown_dd_version=True,
            )

            # Load paths from IDS
            self.ps_reader = read_ps.PlasmaStateReader(self._ids)
            (
                self._selectable_scalar_paths,
                self._selectable_vector_paths,
                self._filled_scalar_paths,
                self._filled_vector_paths,
            ) = self.ps_reader.load_paths_from_ids()
            self._selectable_paths = (
                self._selectable_vector_paths + self._selectable_scalar_paths
            )

    def _name_from_idspath(self, path):
        """Converts an IDSPath to a string by removing 'ggd' and capitalizing each part
        of the path, with the parts separated by spaces.

        Example:
            If path is IDSPath('ggd/electrons/pressure'), the function returns
            "Electrons Pressure"

        Args:
            path: The IDSPath object to convert into a formatted string

        Returns:
            A formatted string of the IDSPath
        """
        path_list = list(path.parts)
        if "ggd" in path_list:
            path_list.remove("ggd")
        for i in range(len(path_list)):
            path_list[i] = path_list[i].capitalize()

        name = " ".join(path_list)

        # If GGDs are not filled in the first time step, add (?) to their name and add
        # zero-width space to move them to the bottom of the list.
        if (
            path not in self._filled_scalar_paths
            and path not in self._filled_vector_paths
        ):
            name = f"\u200B{name} (?)"
        return name

    def _get_selected_ggd_paths(self):
        """Retrieve the IDSPaths of the selected scalar and vector GGD arrays.

        Returns:
            selected_scalar_paths: List of IDSPaths of selected scalar GGD arrays
            selected_vector_paths: List of IDSPaths of selected vector GGD arrays
        """

        # Determine if selected GGD arrays are scalar or vector arrays
        selected_scalar_paths = [
            obj
            for obj in self._selectable_scalar_paths
            if self._name_from_idspath(obj) in self._selected_paths
        ]
        selected_vector_paths = [
            obj
            for obj in self._selectable_vector_paths
            if self._name_from_idspath(obj) in self._selected_paths
        ]
        return selected_scalar_paths, selected_vector_paths

    def _get_selected_time_step(self, outInfo):
        """Retrieves the selected time step index based on the time selection widget in
        the ParaView UI.

        Args:
            outInfo: An information object that holds details about the update request

        Returns:
            time_step_idx: The index of the selected time step. If the selected time
            step cannot be found in the IDS, returns None.
        """
        # Retrieve time step from time selection widget in Paraview UI
        executive = self.GetExecutive()
        outInfo = outInfo.GetInformationObject(0)
        time_step = outInfo.Get(executive.UPDATE_TIME_STEP())
        time_step_idx = np.where(self._time_steps == time_step)[0]

        # Paraview (v5.12.1) contains a bug where, if you load a dataset with only a
        # single timestep, it automatically loads 10 synthetic timesteps ranging from
        # timestep upto time step + 1. The issue can be found here:
        # https://gitlab.kitware.com/paraview/paraview/-/issues/22360
        # A workaround is as follows:
        # 1. load the single timestep data
        # 2. open "View > Time Manager"
        # 3. uncheck and check again the checkbox "Time Sources"
        # For now, return None if such a timestep is chosen
        if len(time_step_idx) == 0:
            return None
        else:
            time_step_idx = time_step_idx[0]
            logger.debug(f"Selected time step: {self._time_steps[time_step_idx]}")

        return time_step_idx


_ggd_types_xml = "".join(
    f'<Entry text="{member.name}" value="{member.index}" />'
    for member in imaspy.identifiers.ggd_identifier
)
_ggd_space_types_xml = "".join(
    f'<Entry text="{member.name}" value="{member.index}" />'
    for member in imaspy.identifiers.ggd_space_identifier
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
