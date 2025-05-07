import logging
from dataclasses import dataclass

import numpy as np
import vtk
from imas.ids_struct_array import IDSStructArray
from imas.ids_structure import IDSStructure
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.vtkCommonDataModel import (
    vtkPartitionedDataSet,
    vtkPartitionedDataSetCollection,
)

from imas_paraview.ids_util import create_name_recursive, get_object_by_name
from imas_paraview.plugins.base_class import GGDVTKPluginBase
from imas_paraview.util import find_closest_indices

logger = logging.getLogger("imas_paraview")

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
class Profiles2DReader(GGDVTKPluginBase, is_time_dependent=True):
    """profiles_2d reader based on IMAS-Python"""

    def __init__(self):
        super().__init__("vtkPartitionedDataSetCollection", PROFILES_2D_IDS_NAMES)
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
            output = vtkPartitionedDataSetCollection.GetData(outInfo)
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
            for profile in self._ids.time_slice[time_idx].profiles_2d:
                self._search_valid_profile(profile)
                # Exit once a valid profiles_2d node has been found
                if self._filled_profiles:
                    break
        else:
            if not time_idx < len(self._ids.profiles_2d):
                logger.error(f"There is no profiles_2d at time step {time_idx}")
                return
            profile = self._ids.profiles_2d[time_idx]
            self._search_valid_profile(profile)

        if not self._filled_profiles:
            logger.error("Could not find a valid profiles_2d node.")
            return

        if self.show_all:
            self._selectable = self._get_selectable_profiles(self._all_profiles)
        else:
            self._selectable = self._get_selectable_profiles(self._filled_profiles)

    def _search_valid_profile(self, profile):
        """Looks for valid profiles within a profiles_2d node and stores them into a
        list.

        Args:
            profile: The profile node at a specific time step.
        """
        self.r = np.array([])
        self.z = np.array([])
        self._filled_profiles = []
        self._all_profiles = []

        if (
            hasattr(profile, "r")
            and hasattr(profile, "z")
            and profile.r.has_value
            and profile.r.shape == profile.z.shape
        ):
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
                return
            self.r = np.tile(profile.grid.dim1, (len(profile.grid.dim2), 1))
            self.z = np.tile(profile.grid.dim2, (len(profile.grid.dim1), 1)).T

        # If the grid is not filled, continue to the next profiles_2d
        if len(self.r) == 0 or len(self.z) == 0:
            return

        self._recursively_find_profiles(profile)

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
        vtkPartitionedDataSetCollection.

        Args:
            output: The vtkPartitionedDataSetCollection containing the profiles.
        """

        # Since profiles share the same grid, it can be converted once
        vtk_points = self._create_vtkpoints()

        for i, profile_name in enumerate(self._selected):
            profile = get_object_by_name(self._selectable, profile_name)
            if profile is None:
                raise ValueError(f"Could not find {profile_name}")

            logger.info(f"Selected {profile.name}")

            vtk_scalars = self._create_vtkscalars(profile)
            vtk_ugrid = self._create_ugrid(vtk_points, vtk_scalars)

            partitioned_dataset = vtkPartitionedDataSet()
            partitioned_dataset.SetPartition(0, vtk_ugrid)

            output.SetPartitionedDataSet(i, partitioned_dataset)

    def _create_vtkpoints(self):
        """Create vtkPoints containing the radial and height coordinates of the profile.

        Returns:
            The VTK points.
        """
        r_flat = self.r.ravel()
        z_flat = self.z.ravel()
        points = np.column_stack((r_flat, np.zeros_like(r_flat), z_flat))

        vtk_points = vtk.vtkPoints()
        vtk_points.SetData(numpy_to_vtk(points.astype(np.float64), deep=True))
        return vtk_points

    def _create_vtkscalars(self, profile):
        """Create VTK array containing the values of the profile at the grid.

        Args:
            profile: The profile to load the values from

        Returns:
            The converted VTK array.
        """
        values_flat = profile.profile.ravel()
        vtk_scalars = numpy_to_vtk(values_flat.astype(np.float64), deep=True)
        units = profile.profile.metadata.units
        vtk_scalars.SetName(f"{profile.name} [{units}]")
        return vtk_scalars

    def _create_ugrid(self, vtk_points, vtk_scalars):
        """Create a vtkUnstructuredGrid of the given profile.

        Args:
            profile: The profile to create a ugrid for.

        Returns:
            The created unstructured grid.
        """

        ugrid = vtk.vtkUnstructuredGrid()
        ugrid.GetPointData().SetScalars(vtk_scalars)

        num_points = vtk_points.GetNumberOfPoints()
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
