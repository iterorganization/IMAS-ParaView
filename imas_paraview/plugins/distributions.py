import logging

import numpy as np
import vtk
from imas import identifiers
from imas.ids_defs import IDS_TIME_MODE_HOMOGENEOUS
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.vtkCommonCore import vtkPoints, vtkStringArray
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkMultiBlockDataSet,
    vtkPolyData,
)

from imas_paraview.paraview_support.servermanager_tools import (
    doublevector,
    propertygroup,
    stringlistdomain,
    stringvector,
)
from imas_paraview.plugins.base_class import GGDVTKPluginBase
from imas_paraview.util import find_closest_indices, pol_to_cart

logger = logging.getLogger("imas_paraview")

SUPPORTED_IDS_NAMES = ["distributions"]


@smproxy.source(label="Distributions Markers Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class DistributionsMarkersReader(GGDVTKPluginBase, is_time_dependent=True):
    _NONE_LABEL = "(none)"

    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", SUPPORTED_IDS_NAMES)
        self.selectable_map = {}
        self._available_coords = set()
        self._x_axis = "x"
        self._y_axis = "y"
        self._z_axis = "z"
        self._x_scale = 1.0
        self._y_scale = 1.0
        self._z_scale = 1.0

    @stringvector(
        name="AxisCoordList", information_only=1, si_class="vtkSIDataArrayProperty"
    )
    def P15_GetAxisCoordList(self):
        """Return the list of available coordinates for axis dropdowns."""
        arr = vtkStringArray()
        for name in self._available_coords:
            arr.InsertNextValue(name)
        return arr

    @stringvector(name="XAxis", label="X-Axis", default_values="x")
    @stringlistdomain("AxisCoordList", name="x_axis_list")
    def P16_SetXAxis(self, coord_name):
        """Select which coordinate to map to the X axis."""
        if self._x_axis != coord_name:
            self._x_axis = coord_name
            self.Modified()

    @stringvector(name="YAxis", label="Y-Axis", default_values="y")
    @stringlistdomain("AxisCoordList", name="y_axis_list")
    def P17_SetYAxis(self, coord_name):
        """Select which coordinate to map to the Y axis."""
        if self._y_axis != coord_name:
            self._y_axis = coord_name
            self.Modified()

    @stringvector(name="ZAxis", label="Z-Axis", default_values="z")
    @stringlistdomain("AxisCoordList", name="z_axis_list")
    def P18_SetZAxis(self, coord_name):
        """Select which coordinate to map to the Z axis."""
        if self._z_axis != coord_name:
            self._z_axis = coord_name
            self.Modified()

    @doublevector(name="XScale", label="Scale X-axis", default_values="1.0")
    def P19_SetXScale(self, value):
        if getattr(self, "_x_scale", 1.0) != value:
            self._x_scale = value
            self.Modified()

    @doublevector(name="YScale", label="Scale Y-axis", default_values="1.0")
    def P20_SetYScale(self, value):
        if getattr(self, "_y_scale", 1.0) != value:
            self._y_scale = value
            self.Modified()

    @doublevector(name="ZScale", label="Scale Z-axis", default_values="1.0")
    def P21_SetZScale(self, value):
        if getattr(self, "_z_scale", 1.0) != value:
            self._z_scale = value
            self.Modified()

    @propertygroup(
        "Axis Coordinate Mapping",
        [
            "AxisCoordList",
            "XAxis",
            "YAxis",
            "ZAxis",
            "XScale",
            "YScale",
            "ZScale",
        ],
    )
    def PG3_AxisGroup(self):
        """Dummy function to define a PropertyGroup."""

    def setup_ids(self):
        """
        Scan the distribution array to find distributions containing marker data,
        populate the selection UI, and collect available coordinate names for the
        axis dropdowns.
        """
        assert self._ids is not None, "IDS cannot be empty during setup."
        self.selectable_map = {}
        coord_names_seen = {self._NONE_LABEL}

        for i, dist in enumerate(self._ids.distribution):
            dist_name = self._create_dist_name(dist)
            if not dist_name:
                dist_name = f"distribution {i}"
            if len(dist.markers) == 0:
                logger.warning("'%s' does not contain any markers, skipping", dist_name)
                continue

            self.selectable_map[dist_name] = dist

            # Collect coordinate names from the first time slice, this assumes the
            # marker coordinates stay the same over time
            markers_slice = dist.markers[0]
            for coord in markers_slice.coordinate_identifier:
                name = str(coord.name)
                if name not in coord_names_seen:
                    coord_names_seen.add(name)

            # If r and phi exist, make sure x and y are available coordinates, they
            # will be calculated upon loading
            if "r" in coord_names_seen and "phi" in coord_names_seen:
                coord_names_seen.add("x")
                coord_names_seen.add("y")

        self._selectable = list(self.selectable_map.keys())
        self._available_coords = coord_names_seen

        self._x_axis = "x" if "x" in coord_names_seen else self._NONE_LABEL
        self._y_axis = "y" if "y" in coord_names_seen else self._NONE_LABEL
        self._z_axis = "z" if "z" in coord_names_seen else self._NONE_LABEL

    def _create_dist_name(self, dist):
        """Generate a name based on the species of the distribution.

        Args:
            dist: distribution IDSStructure
        """
        species = dist.species
        type_index = species.type.index
        ref_id = identifiers.species_reference_identifier

        dist_name = ""
        if dist.species.type.name:
            dist_name = dist.species.type.name

        if type_index in [ref_id.ion.index, ref_id.ion_state.index]:
            ion_name = species.ion.name
            if type_index == ref_id.ion_state.index:
                ion_name = f"{ion_name} ({species.ion.state.name})"
            dist_name = f"{dist_name} ({ion_name})"
        elif type_index in [ref_id.neutral.index, ref_id.neutral_state.index]:
            neutral_name = species.neutral.name
            if type_index == ref_id.neutral_state.index:
                neutral_name = f"{neutral_name} ({species.neutral.state.name})"
            dist_name = f"{dist_name} ({neutral_name})"

        return str(dist_name)

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1
        if self._ids.ids_properties.homogeneous_time != IDS_TIME_MODE_HOMOGENEOUS:
            logger.error("Only IDSs with homogeneous time-mode are supported.")
            return 1
        if len(self._ids.time) == 0:
            logger.error("The IDS does not have a filled time array.")
            return 1

        time = self._get_selected_time_step(outInfo)
        if time is None:
            return 1

        index_list = find_closest_indices([time], self._ids.time)
        time_idx = index_list[0]

        if len(self._selected) > 0:
            output = vtkMultiBlockDataSet.GetData(outInfo)
            self._load_markers(output, time_idx)

        return 1

    def _load_markers(self, output, time_idx):
        """Load selected distributions into the vtkMultiBlockDataSet.

        Args:
            output: The vtkMultiBlockDataSet to populate.
            time_idx: The current time index.
        """
        for block_id, dist_name in enumerate(self._selected):
            dist = self.selectable_map[dist_name]
            markers = dist.markers[time_idx]
            logger.info("Loading %s", dist_name)

            vtk_poly = self._create_vtk_markers(markers)
            output.SetBlock(block_id, vtk_poly)
            output.GetMetaData(block_id).Set(vtkMultiBlockDataSet.NAME(), dist_name)

    def _create_vtk_markers(self, markers):
        """Convert marker positions to vtkPolyData.

        Args:
            markers: The markers structure to convert to vtk.

        Returns:
            vtkPolyData containing marker data
        """
        column_map = {
            str(coord.name): i for i, coord in enumerate(markers.coordinate_identifier)
        }

        num_pts = markers.positions.shape[0]

        xyz = np.column_stack(
            [
                self._resolve_axis(self._x_axis, markers, column_map) * self._x_scale,
                self._resolve_axis(self._y_axis, markers, column_map) * self._y_scale,
                self._resolve_axis(self._z_axis, markers, column_map) * self._z_scale,
            ]
        )

        vtk_pts = vtkPoints()
        vtk_pts.SetData(numpy_to_vtk(xyz, deep=True))

        cells = np.column_stack([np.ones(num_pts), np.arange(num_pts)]).ravel()
        verts = vtkCellArray()
        verts.SetCells(
            num_pts, numpy_to_vtk(cells, deep=True, array_type=vtk.VTK_ID_TYPE)
        )

        poly = vtkPolyData()
        poly.SetPoints(vtk_pts)
        poly.SetVerts(verts)

        point_data = poly.GetPointData()

        weights_array = numpy_to_vtk(markers.weights, deep=True)
        weights_array.SetName("weights")
        point_data.AddArray(weights_array)
        point_data.SetActiveScalars("weights")

        for name, idx in column_map.items():
            arr = numpy_to_vtk(markers.positions[:, idx], deep=True)
            arr.SetName(name)
            point_data.AddArray(arr)

        return poly

    def _resolve_axis(self, axis_name, markers, column_map):
        """Resolve a coordinate axis for marker positions.

        Args:
            axis_name: Name of the axis to resolve.
            markers: The markers structure.
            column_map: Dictionary mapping the coordinate name to column index in
                markers.positions structure.

        Returns:
            nd-array containing the marker.positions coordinate.
        """
        num_pts = markers.positions.shape[0]
        if axis_name == self._NONE_LABEL:
            return np.zeros(num_pts)

        column_index = column_map.get(axis_name)
        if column_index is not None:
            return markers.positions[:, column_index]

        # Calculate x and y coordinates from r and phi, if they are available
        if axis_name in ["x", "y"] and "r" in column_map and "phi" in column_map:
            r_vals = markers.positions[:, column_map["r"]]
            phi_vals = markers.positions[:, column_map["phi"]]
            x_calc, y_calc = pol_to_cart(r_vals, phi_vals)
            if axis_name == "x":
                return x_calc
            elif axis_name == "y":
                return y_calc
