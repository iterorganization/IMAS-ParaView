"""Plugin to view arbitrary 0D time-dependent data from any IDS."""

import logging

import imas
import numpy as np
from imas.ids_data_type import IDSDataType
from imas.ids_defs import (
    IDS_TIME_MODE_HETEROGENEOUS,
    IDS_TIME_MODE_INDEPENDENT,
)
from imas.ids_metadata import IDSType
from imas.ids_struct_array import IDSStructArray
from imas.ids_structure import IDSStructure
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.vtkCommonDataModel import vtkTable

from imas_paraview.ids_util import create_name_recursive
from imas_paraview.plugins.base_class import GGDVTKPluginBase

logger = logging.getLogger("imas_paraview")

# This plugin is generic and should work with any IDS
SUPPORTED_IDS_NAMES = imas.IDSFactory().ids_names()


@smproxy.source(label="0D Time-Dependent Data Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class TimeDependent0DReader(GGDVTKPluginBase, is_time_dependent=True):
    """Reader for arbitrary 0D time-dependent data from any IDS."""

    def __init__(self):
        super().__init__("vtkTable", SUPPORTED_IDS_NAMES)
        self.selectable_map = {}
        self._filled_quantities = []

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        # Retrieve the selected time step
        time = self._get_selected_time_step(outInfo)
        if time is None:
            logger.warning("Selected invalid time step")
            return 1

        if len(self._selected) > 0:
            output = vtkTable.GetData(outInfo)
            self._load_time_dependent_data(output, time)
        return 1

    def setup_ids(self):
        """Scan the IDS for all 0D time-dependent quantities and populate
        the selection list.
        """
        if self._ids is None:
            return

        logger.info(
            "Scanning IDS '%s' for time-dependent 0D data...", self._ids.metadata.name
        )

        self._filled_quantities = []
        self._recursively_find_time_dependent_quantities(self._ids)

        logger.info(
            f"Found {len(self._filled_quantities)} filled time-dependent quantities"
        )

        self.selectable_map = self._get_quantity_names()
        self._selectable = list(self.selectable_map)

    def _get_quantity_names(self):
        selectable_map = {}

        for node in self._filled_quantities:
            name = f"{create_name_recursive(node)} [{node.metadata.units}]"
            selectable_map[name] = node

        return selectable_map

    def _recursively_find_time_dependent_quantities(self, node):
        metadata = node.metadata
        # Time and GGD quantities
        if metadata.name in ("time", "grid_ggd", "grids_ggd", "ggd", "description_ggd"):
            return

        if isinstance(node, IDSStructure) or isinstance(node, IDSStructArray):
            for subnode in node:
                self._recursively_find_time_dependent_quantities(subnode)
                # Only scan the first time slice
                if subnode.metadata.name == "time_slice":
                    break
        elif (
            metadata.data_type in (IDSDataType.FLT, IDSDataType.INT)
            and metadata.type == IDSType.DYNAMIC
            and node.has_value
            and metadata.ndim in [0, 1]
        ):
            self._filled_quantities.append(node)

    def _load_time_dependent_data(self, output, selected_time):
        """Load the selected time-dependent quantities up to the selected time.

        Args:
            output: vtkTable to populate with data
            selected_time: The time step selected in ParaView
        """
        # Get the time array for this IDS
        time_mode = self._ids.ids_properties.homogeneous_time

        if time_mode == IDS_TIME_MODE_INDEPENDENT:
            logger.warning("This IDS has no time-dependent data")
            return
        elif time_mode == IDS_TIME_MODE_HETEROGENEOUS:
            logger.warning("Heterogeneous IDSs are not supported")
            return
        else:
            time_array = self._ids.time

        if len(time_array) == 0:
            logger.warning("Time array is empty")
            return

        time_indices = np.where(time_array <= selected_time)[0]
        if len(time_indices) == 0:
            logger.warning(f"No data available up to time {selected_time}")
            return

        output_times = time_array[time_indices]

        time_vtk = numpy_to_vtk(output_times, deep=1)
        time_vtk.SetName("Time [s]")
        output.AddColumn(time_vtk)

        for quantity_name in self._selected:
            node = self.selectable_map[quantity_name]
            quantity_values = self._get_quantity_values(node, output_times)

            # Create VTK array
            data_vtk = numpy_to_vtk(quantity_values, deep=1)
            data_vtk.SetName(quantity_name)
            output.AddColumn(data_vtk)

            logger.info(f"Loaded {len(quantity_values)} points for '{quantity_name}'")

    def _get_quantity_values(self, node, output_times):
        if node.metadata.ndim == 1:
            quantity_values = node[: min(len(output_times), len(node))]
        else:
            full_path = imas.util.get_full_path(node)
            parts = full_path.split("time_slice", 1)
            time_slice_path = parts[0] + "time_slice"
            path_in_slice = parts[1].lstrip("[0]/")

            quantity_values = []
            for i in range(len(output_times)):
                node_in_slice = self._ids[time_slice_path][i][path_in_slice]
                quantity_values.append(node_in_slice)
        return quantity_values
