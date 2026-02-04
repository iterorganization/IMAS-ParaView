"""Plugin to visualize 0D time-dependent data from any IDS."""

import logging
from dataclasses import dataclass
from typing import Optional

import imas
import numpy as np
from imas.ids_base import IDSBase
from imas.ids_data_type import IDSDataType
from imas.ids_defs import IDS_TIME_MODE_HOMOGENEOUS
from imas.ids_metadata import IDSType
from imas.ids_struct_array import IDSStructArray
from imas.ids_structure import IDSStructure
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.vtkCommonDataModel import vtkTable

from imas_paraview.ids_util import create_name_recursive, is_time_dependent_aos
from imas_paraview.plugins.base_class import GGDVTKPluginBase

logger = logging.getLogger("imas_paraview")

# This plugin is generic and should work with any IDS
SUPPORTED_IDS_NAMES = imas.IDSFactory().ids_names()


@dataclass
class FilledQuantity:
    """Stores information about a filled time-dependent quantity inside a time-dependent
    AoS.

    For example, the node `equilibrium['time_slice[0]/global_quantities/ip']`
    will contain time_slice = `equilibrium['time_slice']`,
    and `remaining_path = "global_quantities/ip"`
    """

    node: IDSBase  # The filled time-dependent IDS node
    time_slice: Optional[IDSBase] = None  # Parent time-dependent AoS (e.g., time_slice)
    remaining_path: Optional[str] = None  # Path from time_slice to node


# TODO: rename reader
@smproxy.source(label="Scalar Time Trace Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class ScalarTimeTraceReader(GGDVTKPluginBase, is_time_dependent=True):
    """Reader for visualizing time-dependent scalar data from any IDS."""

    def __init__(self):
        super().__init__("vtkTable", SUPPORTED_IDS_NAMES)
        self._filled_quantities_map = {}

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

        self._filled_quantities_map = {}
        self._recursively_find_time_dependent_quantities(self._ids)

        logger.info(
            "Found %d filled time-dependent quantities",
            len(self._filled_quantities_map),
        )
        self._selectable = list(self._filled_quantities_map)

    def _recursively_find_time_dependent_quantities(
        self, node, time_slice=None, path_from_time_slice=""
    ):
        metadata = node.metadata
        # Skip time and GGD quantities
        if metadata.name in ("time", "grid_ggd", "grids_ggd", "ggd", "description_ggd"):
            return

        parent = imas.util.get_parent(node)
        if parent is not None and is_time_dependent_aos(parent):
            # This becomes our new time slice reference point
            time_slice = parent
            path_from_time_slice = ""

        if isinstance(node, (IDSStructure, IDSStructArray)):
            for i, subnode in enumerate(node):
                if isinstance(node, IDSStructArray):
                    new_path = f"{path_from_time_slice}[{i}]"
                else:
                    new_path = (
                        f"{path_from_time_slice}/{subnode.metadata.name}"
                        if path_from_time_slice
                        else subnode.metadata.name
                    )
                self._recursively_find_time_dependent_quantities(
                    subnode, time_slice=time_slice, path_from_time_slice=new_path
                )
                if is_time_dependent_aos(node):  # Only scan the first time slice
                    break
        elif (
            metadata.data_type in (IDSDataType.FLT, IDSDataType.INT)
            and metadata.type == IDSType.DYNAMIC
            and node.has_value
            and (metadata.ndim == 0 or (metadata.ndim == 1 and time_slice is None))
        ):
            name = f"{create_name_recursive(node)} [{node.metadata.units}]"
            self._filled_quantities_map[name] = FilledQuantity(
                node, time_slice, path_from_time_slice
            )

    def _load_time_dependent_data(self, output, selected_time):
        """Load the selected time-dependent quantities up to the selected time.

        Args:
            output: vtkTable to populate with data
            selected_time: The time step selected in ParaView
        """
        if self._ids.ids_properties.homogeneous_time != IDS_TIME_MODE_HOMOGENEOUS:
            logger.warning("Only IDSs homogeneous time-mode are supported.")
            return
        time_array = self._ids.time

        if len(time_array) == 0:
            logger.warning("Time array is empty")
            return

        time_indices = np.where(time_array <= selected_time)[0]
        if len(time_indices) == 0:
            logger.warning("No data available up to time %f", selected_time)
            return

        output_times = time_array[time_indices]

        time_vtk = numpy_to_vtk(output_times, deep=1)
        time_vtk.SetName("Time [s]")
        output.AddColumn(time_vtk)

        for quantity_name in self._selected:
            node = self._filled_quantities_map[quantity_name]
            quantity_values = self._get_quantity_values(node, len(output_times))

            # Create VTK array
            data_vtk = numpy_to_vtk(quantity_values, deep=1)
            data_vtk.SetName(quantity_name)
            output.AddColumn(data_vtk)

            logger.info(
                "Loaded %d points for '%s'", len(quantity_values), quantity_name
            )

    def _get_quantity_values(self, quantity, n_steps):
        if quantity.node.metadata.ndim == 0:
            time_slice = quantity.time_slice
            remaining_path = quantity.remaining_path
            quantity_values = [time_slice[i][remaining_path] for i in range(n_steps)]
        else:
            quantity_values = quantity.node[:n_steps]
        return quantity_values
