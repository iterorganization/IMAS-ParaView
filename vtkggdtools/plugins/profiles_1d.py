"""IMASPy plugin to view profiles_1D nodes"""

import logging
from dataclasses import dataclass

import numpy as np
from imaspy.ids_struct_array import IDSStructArray
from imaspy.ids_structure import IDSStructure
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.vtkCommonDataModel import vtkTable

from vtkggdtools.ids_util import create_name_recursive
from vtkggdtools.plugins.base_class import GGDVTKPluginBase
from vtkggdtools.util import find_closest_indices

logger = logging.getLogger("vtkggdtools")

PROFILES_1D_IDS_NAMES = ["core_profiles", "core_sources"]


@dataclass
class Profile_1d:
    """Data class that stores 1d profiles, along with its name and coordinate array."""

    name: str
    coordinates: np.ndarray
    profile: np.ndarray


@smproxy.source(label="profiles_1d Reader")
@smhint.xml("""<ShowInMenu category="VTKGGDTools" />""")
class IMASPyProfiles1DReader(GGDVTKPluginBase):
    """profiles_1d reader based on IMASPy"""

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
            profile = self.find_profile(profile_name, self._selectable)
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

    def find_profile(self, name, profiles_list):
        """Return the profile in the list that has the provided name.
        If no match is found, None is returned instead.

        Args:
            name: Name of the 1d profile to search for
            profiles_list: list of Profile_1d dataclasses to search through

        Returns:
            Profile_1d dataclass containing the given name
        """
        for profile in profiles_list:
            if profile.name == name:
                return profile
        return None

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

    def request_information(self):
        """
        Select which profiles to show in the array domain selector, based
        on whether the "Show All" checkbox is enabled.
        """
        if self._ids is not None:
            self._filled_profiles = []
            self._all_profiles = []
            time_idx = 0
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
        Called after an IDS is loaded for the first time. Intentionally left empty.
        """
        pass

    def _get_profiles(self, profiles):
        """Filters and processes a list of profiles to extract those containing coordinates and
        creates a list of Profile_1d objects with generated names.

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
