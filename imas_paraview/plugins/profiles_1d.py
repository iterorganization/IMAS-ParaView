"""Plugin to view profiles_1D nodes"""

import logging
from dataclasses import dataclass

import numpy as np
from imas.ids_struct_array import IDSStructArray
from imas.ids_structure import IDSStructure
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.vtkCommonDataModel import vtkTable

from imas_paraview.ids_util import create_name_recursive, get_object_by_name
from imas_paraview.plugins.base_class import GGDVTKPluginBase
from imas_paraview.util import find_closest_indices

logger = logging.getLogger("imas_paraview")

PROFILES_1D_IDS_NAMES = ["core_profiles", "core_sources"]


@dataclass
class Profile_1d:
    """Data class that stores 1d profiles, along with its name and coordinate array."""

    name: str
    coordinates: np.ndarray
    profile: np.ndarray


@smproxy.source(label="1D Profiles Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class Profiles1DReader(GGDVTKPluginBase, is_time_dependent=True):
    """profiles_1d reader based on IMAS-Python"""

    def __init__(self):
        super().__init__("vtkTable", PROFILES_1D_IDS_NAMES)

    def GetAttributeArrayName(self, idx) -> str:
        return self._selectable[idx].name

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        # Retrieve the selected time step and profiles
        time = self._get_selected_time_step(outInfo)
        if time is None:
            logger.warning("Selected invalid time step")
            return 1

        index_list = find_closest_indices([time], self._ids.time)
        time_idx = index_list[0]
        self.update_available_profiles(time_idx)

        if len(self._selected) > 0:
            output = vtkTable.GetData(outInfo)
            self._load_profiles(output)
        return 1

    def _load_profiles(self, output):
        """Creates vtkDoubleArrays for the selected profiles and stores them into
        a vtkTable. If multiple profiles are selected, it is checked if their
        coordinates match, if they do they can be plotted in the same 1d plot.

        Args:
            output: vtkTable output of the plugin containing the profiles and their
                coordinates as columns.
        """
        prev_x = None
        for profile_name in self._selected:
            profile = get_object_by_name(self._selectable, profile_name)
            if profile is None:
                logger.warning(
                    f"Could not find a matching profile with name {profile_name}"
                )
                continue

            if len(profile.profile) == 0:
                logger.warning(f"The selected profile {profile.name} is empty.")
                continue
            if len(profile.coordinates) != len(profile.profile):
                logger.warning(
                    "The length of the linked coordinate array does not match."
                )
                continue

            logger.info(f"Selected {profile_name}.")
            y_values = self._create_vtk_double_array(profile.profile, profile.name)
            output.AddColumn(y_values)

            # If multiple profiles are selected, check if their coordinates match
            if prev_x is None:
                x_values = self._create_vtk_double_array(
                    profile.coordinates, profile.coordinates.metadata.name
                )
                output.AddColumn(x_values)
            else:
                if profile.coordinates is not prev_x:
                    raise RuntimeError(
                        "The X values for the selected profiles do not match. "
                        "Select the profiles one by one instead."
                    )
            prev_x = profile.coordinates

    def _create_vtk_double_array(self, values, name):
        """Creates a vtkDoubleArray with the given name and values.

        Args:
            values: values to store in the array
            name: name to give to the array

        Returns:
            vtkDoubleArray with given name and values
        """
        vtk_array = numpy_to_vtk(values, deep=1)
        vtk_array.SetName(name)
        return vtk_array

    def update_available_profiles(self, time_idx):
        """
        Searches through the profiles_1d node at the current time step for
        available profiles to select. Which profiles show in the array domain selector
        is based on whether the "Show All" checkbox is enabled.
        """

        if self._ids is not None:
            self._filled_profiles = []
            self._all_profiles = []
            if self._ids.metadata.name == "core_profiles":
                for profile_node in self._ids.profiles_1d[time_idx]:
                    self._recursive_find_profiles(profile_node)
            elif self._ids.metadata.name == "core_sources":
                for source in self._ids.source:
                    for profile_node in source.profiles_1d[time_idx]:
                        self._recursive_find_profiles(profile_node)
            else:
                raise NotImplementedError(
                    "Currently only the 1D profiles of the 'core_profiles' and "
                    "'core_sources' are supported"
                )
        if self.show_all:
            self._selectable = self._get_profiles(self._all_profiles)
        else:
            self._selectable = self._get_profiles(self._filled_profiles)

    def setup_ids(self):
        """
        Select which profiles to show in the array domain selector, based
        on whether the "Show All" checkbox is enabled.
        """
        # WARN: The selected time cannot be fetched during the RequestInformation
        # step, so we take the first time index here to update the domain selection
        # array. This causes some issues if there are profiles in the IDS which are
        # not filled in later time steps. For example, if the ion for
        # profiles_1d[0]/ion[0].density has name "A", but for
        # profiles_1d[1]/ion[0].density the ion has the name "B", this profile is
        # lost. To mitigate this, ensure that the ions which will be used are
        # all defined at the first time step. Alternatively, a current workaround
        # for this is to press the "Show All" checkbox twice (i.e. enable and then
        # disable), when a later time step is selected. This will cause the UI
        # to update and show the available profiles at the selected time step.
        time_idx = 0
        self.update_available_profiles(time_idx)

    def _get_profiles(self, profiles):
        """Filters and processes a list of profiles to extract those containing
        coordinates and creates a list of Profile_1d objects with generated names.

        Args:
            profiles: A list of profile objects to be processed.

        Returns:
            A list of `Profile_1d` objects created from the provided profiles that
                  contain valid coordinates.
        """
        profile_list = []
        for profile in profiles:
            # Only store the profile if it contains coordinates
            if profile.metadata.coordinate1.references:
                path = profile.metadata.coordinate1.references[0]
                coordinates = path.goto(profile)
                name = create_name_recursive(profile)

                profile = Profile_1d(name, coordinates, profile)
                profile_list.append(profile)
        return profile_list

    def _recursive_find_profiles(self, node):
        """Recursively traverses through the IDS node searching for filled 1d profiles.
        If a filled profile is found, it is appended to self._filled_profiles.

        Args:
            node: the node to search through.
        """
        if isinstance(node, IDSStructure) or isinstance(node, IDSStructArray):
            for subnode in node:
                self._recursive_find_profiles(subnode)
        else:
            if node.metadata.ndim == 1:
                if node.has_value:
                    self._filled_profiles.append(node)
                self._all_profiles.append(node)
