"""IMASPy version of the paraview plugin classes.
"""

import getpass
import logging
from abc import ABC, abstractmethod

import imaspy
import imaspy.ids_defs
from paraview.util.vtkAlgorithm import smdomain, smhint, smproperty
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonCore import vtkStringArray

from vtkggdtools.imas_uri import uri_from_path, uri_from_pulse_run
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

logger = logging.getLogger("vtkggdtools")

BACKENDS = {
    "MDSplus": imaspy.ids_defs.MDSPLUS_BACKEND,
    "HDF5": imaspy.ids_defs.HDF5_BACKEND,
    "ASCII": imaspy.ids_defs.ASCII_BACKEND,
}
"""Mapping of UI labels for each backend and their ID, used for the Backend dropdown."""
DEFAULT_BACKEND = imaspy.ids_defs.MDSPLUS_BACKEND
"""Default backend selected in the UI."""


class GGDVTKPluginBase(VTKPythonAlgorithmBase, ABC):
    """GGD Reader based on IMASPy"""

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        # Copy all function definitions to the subclass
        # This allows paraview to see that they are properties that should be exposed
        for name, value in vars(GGDVTKPluginBase).items():
            if (
                callable(value)
                and not name.startswith("__")
                and not name == "GetAttributeArrayName"
                and not name == "_ensure_ids"
            ):
                setattr(cls, name, value)

    def __init__(self, output_type, supported_ids):
        super().__init__(
            nInputPorts=0,
            nOutputPorts=1,
            outputType=output_type,
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
        self._supported_ids = supported_ids

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

        # Values to fill the array selector with
        self._selectable = []
        self._selected = []

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
                    self._selectable = []
                    self._ids_list = []
            self._update_ids_list()
            self.Modified()

    def _update_ids_list(self) -> None:
        """Update the list of available IDSs in the selected Data Entry."""
        self._ids_list = []
        if self._dbentry is not None:
            for ids_name in self._supported_ids:
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

    @abstractmethod
    def GetAttributeArrayName(self, idx) -> str:
        pass

    def GetAttributeArrayStatus(self, *args):
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
        "Select IDS",
        ["IDSAndOccurrence", "IDSList", "AttributeArraySelector", "LazyLoading"],
    )
    def PG1_IDSGroup(self):
        """Dummy function to define a PropertyGroup."""

    @propertygroup("Bezier interpolation settings", ["N plane", "Phi range"])
    def PG2_BezierGroup(self):
        """Dummy function to define a PropertyGroup."""

    # Implement VTK algorithm
    ####################################################################################

    def RequestInformation(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence:
            return 1

        # Load IDS and available time steps
        idsname, _, _ = self._ids_and_occurrence.partition("/")
        if idsname not in self._ids_list:
            logger.warning("Could not find the selected IDS.")
            self._selectable = []
            return 1
        self._ensure_ids()

        # TODO: Add support for IDSs with heterogeneous time mode
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

    @abstractmethod
    def _ensure_ids(self):
        pass

    def _get_selected_time_step(self, outInfo):
        """Retrieves the selected time step based on the time selection widget in
        the ParaView UI.

        Args:
            outInfo: An information object that holds details about the update request

        Returns:
            time_step: The selected time step. If the selected time step cannot be found
            in the IDS, returns None.
        """
        # Retrieve time step from time selection widget in Paraview UI
        executive = self.GetExecutive()
        outInfo = outInfo.GetInformationObject(0)
        time_step = outInfo.Get(executive.UPDATE_TIME_STEP())

        # Paraview (v5.12.1) contains a bug where, if you load a dataset with only a
        # single timestep, it automatically loads 10 synthetic timesteps ranging from
        # timestep upto time step + 1. The issue can be found here:
        # https://gitlab.kitware.com/paraview/paraview/-/issues/22360
        # A workaround is as follows:
        # 1. load the single timestep data
        # 2. open "View > Time Manager"
        # 3. uncheck and check again the checkbox "Time Sources"
        # For now, return None if such a timestep is chosen

        # Check if it exists in list of time steps
        if time_step in self._time_steps:
            logger.debug(f"Selected time step in Paraview: {time_step}")
            return time_step
        else:
            return None
