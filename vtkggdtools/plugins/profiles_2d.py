"""IMASPy plugin to view profiles_2d nodes"""

import logging
from dataclasses import dataclass

import numpy as np
import vtk
from imaspy.ids_struct_array import IDSStructArray
from imaspy.ids_structure import IDSStructure
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from vtkggdtools.ids_util import create_name_recursive, get_object_by_name
from vtkggdtools.plugins.base_class import GGDVTKPluginBase
from vtkggdtools.util import find_closest_indices

logger = logging.getLogger("vtkggdtools")

# TODO: waves and distributions are currently not supported as their profiles_2d are
# stored differently
PROFILES_2D_IDS_NAMES = [
    "core_profiles",
    "equilibrium",
    "plasma_initiation",
    "plasma_profiles",
]


@dataclass
class Profile_2d:
    """Data class that stores a profiles_2d node, along with its name."""

    name: str
    profile: np.ndarray


@smproxy.source(label="2D Profiles Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class IMASPyProfiles2DReader(GGDVTKPluginBase, is_time_dependent=True):
    """profiles_2d reader based on IMASPy"""

    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", PROFILES_2D_IDS_NAMES)
        self.r = np.array([])
        self.z = np.array([])
        self._filled_profiles = []
        self._all_profiles = []

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
            output = vtkMultiBlockDataSet.GetData(outInfo)
            self._load_profiles(output)
        return 1

    def update_available_profiles(self, time_idx):
        """Searches through the profiles_2d node at the current time step for
        available profiles to select. Which profiles are shown in the array domain
        selector is based on whether the "Show All" checkbox is enabled.

        Args:
            time_idx: The time index of the 2D profile
        """
        if self._ids is None:
            logger.error("Could not find the IDS.")
            return

        if self._ids.metadata.name == "equilibrium":
            if not time_idx < len(self._ids.time_slice):
                logger.error(f"There is no profiles_2d at time step {time_idx}")
                return
            profiles_2d = self._ids.time_slice[time_idx].profiles_2d
        else:
            if not time_idx < len(self._ids.profiles_2d):
                logger.error(f"There is no profiles_2d at time step {time_idx}")
                return
            profiles_2d = self._ids.profiles_2d[time_idx]

        self._search_valid_profile(profiles_2d)

        if not self._filled_profiles:
            logger.error("Could not find a valid profiles_2d node.")
            return

        if self.show_all:
            self._selectable = self._get_selectable_profiles(self._all_profiles)
        else:
            self._selectable = self._get_selectable_profiles(self._filled_profiles)

    def _search_valid_profile(self, profiles_2d):
        """Looks for valid profiles within a profiles_2d node and stores them into a
        list.

        Args:
            profiles_2d: The profiles_2d node at a specific time step.
        """
        for profile in profiles_2d:
            self.r = np.array([])
            self.z = np.array([])
            self._filled_profiles = []
            self._all_profiles = []

            if profile.r.has_value and profile.r.shape == profile.z.shape:
                # Use the 2D R and Z arrays to construct the grid
                self.r = profile.r
                self.z = profile.z
            else:
                # Check if this is a rectangular grid, and use dim1/dim2
                grid_type = profile.grid_type.index
                if grid_type != 1:
                    logger.debug(
                        f"Found a grid type with identifier index of {grid_type}. "
                        "Only rectangular profiles (index = 1) are supported."
                    )
                    continue
                self.r = np.tile(profile.grid.dim1, (len(profile.grid.dim2), 1))
                self.z = np.tile(profile.grid.dim2, (len(profile.grid.dim1), 1)).T

            self._recursively_find_profiles(profile)

            # Only load the first encountered valid profile
            if self._filled_profiles and len(self.r) > 0 and len(self.z) > 0:
                break

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

    def _get_selectable_profiles(self, profiles):
        """Creates a list of Profile_2d objects with generated names.

        Args:
            profiles: A list of profile objects to be processed.

        Returns:
            A list of `Profile_2d` objects created from the provided profiles
        """
        profile_list = []
        for profile in profiles:
            name = create_name_recursive(profile)
            profile = Profile_2d(name, profile)
            profile_list.append(profile)
        return profile_list

    def _recursively_find_profiles(self, node):
        """Recursively traverses through the IDS node searching for filled 2d profiles.
        If a filled profile is found, it is appended to self._filled_profiles, all
        profiles are appended to self._all_profiles, filled or not.

        Args:
            node: the node to search through.
        """

        if node.metadata.name in ("grid", "r", "z"):
            return
        elif isinstance(node, IDSStructure) or isinstance(node, IDSStructArray):
            for subnode in node:
                self._recursively_find_profiles(subnode)
        else:
            if node.metadata.ndim == 2:
                if node.has_value:
                    self._filled_profiles.append(node)
                self._all_profiles.append(node)

    def _load_profiles(self, output):
        """Go through the list of selected profiles, and load each of them in a
        separate vtkUnstructuredGrid object, which are all combined into a
        vtkMultiBlockDataSet.

        Args:
            output: The vtkMultiBlockDataSet containing the limiter contours.
        """
        for i, profile_name in enumerate(self._selected):
            profile = get_object_by_name(self._selectable, profile_name)
            if profile is None:
                raise ValueError(f"Could not find {profile_name}")

            logger.info(f"Selected {profile.name}")
            vtk_ugrid = self._create_ugrid(profile)
            output.SetBlock(i, vtk_ugrid)

    def _create_ugrid(self, profile):
        """Create a vtkUnstructuredGrid of the given profile.

        Args:
            profile: The profile to create a ugrid for.

        Returns:
            The created unstructured grid.
        """

        r_flat = self.r.ravel()
        z_flat = self.z.ravel()
        values_flat = profile.profile.ravel()

        points = np.column_stack((r_flat, np.zeros_like(r_flat), z_flat))
        vtk_points = vtk.vtkPoints()
        vtk_points.SetData(numpy_to_vtk(points.astype(np.float64), deep=True))
        vtk_scalars = numpy_to_vtk(values_flat.astype(np.float64), deep=True)
        vtk_scalars.SetName(profile.name)

        ugrid = vtk.vtkUnstructuredGrid()
        ugrid.SetPoints(vtk_points)
        ugrid.GetPointData().SetScalars(vtk_scalars)

        num_points = points.shape[0]
        cells = vtk.vtkCellArray()

        for i in range(num_points):
            vertex = vtk.vtkVertex()
            vertex.GetPointIds().SetId(0, i)
            cells.InsertNextCell(vertex)

        ugrid = vtk.vtkUnstructuredGrid()
        ugrid.SetPoints(vtk_points)
        ugrid.GetPointData().SetScalars(vtk_scalars)
        ugrid.SetCells(vtk.VTK_VERTEX, cells)
        return ugrid
