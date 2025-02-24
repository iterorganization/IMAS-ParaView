"""IMASPy plugin to view profiles_2D nodes"""

import logging
from dataclasses import dataclass

import numpy as np
import vtk
from imaspy.ids_structure import IDSStructure
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from vtkggdtools.ids_util import create_name_recursive, get_object_by_name
from vtkggdtools.plugins.base_class import GGDVTKPluginBase
from vtkggdtools.util import find_closest_indices

logger = logging.getLogger("vtkggdtools")

PROFILES_2D_IDS_NAMES = [
    "core_profiles",
    "equilibrium",
    "plasma_initiation",
    "plasm_profiles",
]


@dataclass
class Profile_2d:
    """Data class that stores 2d profiles, along with its name."""

    name: str
    profile: np.ndarray


@smproxy.source(label="2D Profiles Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class IMASPyProfiles2DReader(GGDVTKPluginBase, is_time_dependent=True):
    """profiles_2d reader based on IMASPy"""

    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", PROFILES_2D_IDS_NAMES)

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

    def _load_profiles(self, output):
        """Go through the list of selected limiters, and load each of them in a
        separate vtkPolyData object, which are all combined into a vtkMultiBlockDataSet.

        Args:
            output: The vtkMultiBlockDataSet containing the limiter contours.
        """
        for i, profile_name in enumerate(self._selected):
            profile = get_object_by_name(self._selectable, profile_name)
            if profile is None:
                raise ValueError(f"Could not find {profile_name}")

            logger.info(f"Selected {profile.name}")
            vtk_poly = self._create_ugrid(profile)
            output.SetBlock(i, vtk_poly)

    def _create_ugrid(self, profile):
        """Create a contour based on the r,z coordinates in the limiter.
        The r,z-coordinates are stored as vtkPoints, and connected using vtkLines, which
        are both stored in a vtkPolyData object. If the contour is closed, the start
        and end points are connected.

        Args:
            limiter: limiter object containing contour

        Returns:
            vtkPolyData containing contour data
        """

        r_flat = self.r.flatten()
        z_flat = self.z.flatten()
        values_flat = profile.profile.flatten()

        points = vtk.vtkPoints()
        ugrid = vtk.vtkUnstructuredGrid()
        scalars = vtk.vtkDoubleArray()
        scalars.SetName(profile.name)

        point_ids = []
        for i in range(len(r_flat)):
            pid = points.InsertNextPoint(r_flat[i], 0, z_flat[i])
            scalars.InsertNextValue(values_flat[i])
            point_ids.append(pid)

        ugrid.SetPoints(points)
        ugrid.GetPointData().SetScalars(scalars)
        cells = vtk.vtkCellArray()
        for pid in point_ids:
            vertex = vtk.vtkVertex()
            vertex.GetPointIds().SetId(0, pid)
            cells.InsertNextCell(vertex)

        ugrid.SetCells(vtk.VTK_VERTEX, cells)

        return ugrid

    def update_available_profiles(self, time_idx):
        """
        Searches through the profiles_1d node at the current time step for
        available profiles to select. Which profiles show in the array domain selector
        is based on whether the "Show All" checkbox is enabled.
        """
        self.r = None
        self.z = None
        self._filled_profiles = []
        self._all_profiles = []
        if self._ids is None:
            logger.error("Could not find the IDS.")
            return

        if self._ids.metadata.name == "equilibrium":
            profiles_2d = self._ids.time_slice[time_idx].profiles_2d
        else:
            profiles_2d = self._ids.profiles_2d[time_idx]

        for profile in profiles_2d:
            grid_type = profile.grid_type.index
            if grid_type != 1:
                logger.warning(
                    f"Found a grid type with identifier index of {grid_type}. "
                    "Only rectangular profiles (index = 1) are supported."
                )
                continue

            self._get_filled_quantities(profile)

            # Only load the first encountered valid profile
            if self._filled_profiles:
                break

        if not self._filled_profiles:
            logger.error("Could not find a valid profiles_2d node.")
            return

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
            name = create_name_recursive(profile)
            profile = Profile_2d(name, profile)
            profile_list.append(profile)
        return profile_list

    def _get_filled_quantities(self, profile):
        for plasma_quantity in profile:
            if isinstance(plasma_quantity, IDSStructure):
                continue

            if plasma_quantity.metadata.name == "r":
                self.r = plasma_quantity
                continue
            elif plasma_quantity.metadata.name == "z":
                self.z = plasma_quantity
                continue
            if plasma_quantity.has_value:
                self._filled_profiles.append(plasma_quantity)
            self._all_profiles.append(plasma_quantity)

            # Check if r- and z- coordinates are filled, otherwise take dim1 and dim2
            # as r and z, respectively.
            if self._filled_profiles:
                if self.r is None:
                    if profile.grid.dim1:
                        logger.info(
                            "Could not find 'r' node in the profiles_2d. "
                            "Setting dim1 as the radial coordinate."
                        )
                        self.r = np.tile(
                            profile.grid.dim1[:, np.newaxis],
                            (1, len(profile.grid.dim2)),
                        )
                    else:
                        logger.error(
                            "Could not find 'r' or 'grid.dim1' node in the profiles_2d. "
                        )
                        self._filled_profiles = self._all_profiles = []
                if self.z is None:
                    if profile.grid.dim2:
                        logger.info(
                            "Could not find 'z' node in the profiles_2d. "
                            "Setting dim2 as the radial coordinate."
                        )
                        self.z = np.tile(
                            profile.grid.dim2[np.newaxis, :],
                            (len(profile.grid.dim1), 1),
                        )

                    else:
                        logger.error(
                            "Could not find 'z' or 'grid.dim2' node in the profiles_2d. "
                        )
                        self._filled_profiles = self._all_profiles = []

                return
