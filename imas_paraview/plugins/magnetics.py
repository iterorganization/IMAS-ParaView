"""Plugin to visualize geometries from the magnetics IDS.

Supports:
- flux_loop
- b_field_pol_probe
- b_field_phi_probe
- rogowski_coil
"""

import logging
from dataclasses import dataclass
from typing import Literal

from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from imas_paraview.ids_util import cyl_vector_has_value
from imas_paraview.paraview_support.servermanager_tools import (
    doublevector,
    propertygroup,
)
from imas_paraview.plugins.base_class import GGDVTKPluginBase
from imas_paraview.util import (
    angles_to_vectors,
    create_vtk_arrows,
    ensure_unique_name,
    points_to_vtkpoly,
    pol_to_cart,
)

logger = logging.getLogger("imas_paraview")

SUPPORTED_IDS_NAMES = ["magnetics"]


@dataclass
class FluxLoop:
    positions: list


@dataclass
class RogowskiCoil:
    positions: list


@dataclass
class BFieldProbe:
    r: float
    phi: float
    z: float
    poloidal_angle: float
    toroidal_angle: float
    kind: Literal["pol", "phi"]


@smproxy.source(label="Magnetics Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class MagneticsReader(GGDVTKPluginBase):
    def __init__(self):
        super().__init__("vtkMultiBlockDataSet", SUPPORTED_IDS_NAMES)
        self.selectable_map = {}
        self.b_probe_length = 1.0
        """Length of the arrow drawn for each B-field probe."""

    @doublevector(
        label="Probe Arrow Length (m)",
        name="b_probe_length",
        default_values=1.0,
    )
    def P99_SetSensorLength(self, val):
        """Total length of the arrow used to visualize each B-field probe.
        The arrow tip is placed b_probe_length away from the coil centre along
        the sensor normal axis *n*."""
        self._update_property("b_probe_length", val)

    @propertygroup("Magnetics Reader Settings", ["b_probe_length"])
    def PG3_MagneticsReaderGroup(self):
        """Dummy function to define a PropertyGroup."""

    def setup_ids(self):
        """Scan the magnetics IDS and populate the selectable map with all
        flux loops, B-field probes and Rogowski coils that carry geometry."""
        assert self._ids is not None, "IDS cannot be empty during setup."
        self.selectable_map = {}

        self._populate_flux_loops()
        self._populate_b_field_probes()
        self._populate_rogowski_coils()

        self._selectable = list(self.selectable_map.keys())

    def _populate_flux_loops(self):
        """Populate selectable_map with flux loops."""
        for i, loop in enumerate(self._ids.flux_loop):
            name = self._create_name("Flux loop", loop, i)
            if len(loop.position) == 0:
                logger.warning("'%s' has no position points, skipping.", name)
                continue
            self.selectable_map[name] = FluxLoop(loop.position)

    def _populate_rogowski_coils(self):
        """Populate selectable_map with Rogowski coils."""
        for i, coil in enumerate(self._ids.rogowski_coil):
            name = self._create_name("Rogowski coil", coil, i)
            if len(coil.position) < 2:
                logger.warning("'%s' has fewer than 2 points, skipping.", name)
                continue
            self.selectable_map[name] = RogowskiCoil(coil.position)

    def _populate_b_field_probes(self):
        """Populate selectable_map with poloidal and toroidal B-field probes."""
        for i, probe in enumerate(self._ids.b_field_pol_probe):
            name = self._create_name("Poloidal field probe", probe, i)
            self._add_b_field_probe(name, probe, "pol")

        # Toroidal field probes were renamed in DD 3.42.0
        b_field_phi_probes = (
            getattr(self._ids, "b_field_phi_probe", None) or self._ids.b_field_tor_probe
        )

        for i, probe in enumerate(b_field_phi_probes):
            name = self._create_name("Toroidal field probe", probe, i)
            self._add_b_field_probe(name, probe, "phi")

    def _add_b_field_probe(self, name, probe, kind):
        """Add a B-field probe to the selectable map if it has valid position data.

        Args:
            name: Name of the probe.
            probe: B-field probe IDS node.
            kind: Either "pol" or "phi" for poloidal or toroidal probe.
        """
        if not cyl_vector_has_value(probe.position):
            logger.warning(
                "'%s' does not have position coordinates filled, skipping.", name
            )
            return

        self.selectable_map[name] = BFieldProbe(
            r=probe.position.r,
            phi=probe.position.phi,
            z=probe.position.z,
            poloidal_angle=(
                probe.poloidal_angle if probe.poloidal_angle.has_value else 0.0
            ),
            toroidal_angle=(
                probe.toroidal_angle if probe.toroidal_angle.has_value else 0.0
            ),
            kind=kind,
        )

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        if self._selected:
            output = vtkMultiBlockDataSet.GetData(outInfo)
            self._convert_to_vtk(output)
        return 1

    def _convert_to_vtk(self, output: vtkMultiBlockDataSet):
        """Convert each selected device to VTK geometry.

        Args:
            output: vtkMultiBlockDataSet containing separate block for each device.
        """
        for block_id, name in enumerate(self._selected):
            selected = self.selectable_map[name]

            if isinstance(selected, (FluxLoop, RogowskiCoil)):
                vtk_geom = self._create_loop(selected)
            else:  # Bfield probe
                vtk_geom = self._create_b_field_probe_arrow(selected)

            output.SetBlock(block_id, vtk_geom)

    def _create_b_field_probe_arrow(self, data):
        """Create a VTK arrow representation of a B-field probe.

        Args:
            data: BFieldProbe dataclass containing position and angle information.

        Returns:
            VTK polydata representing the B-field probe as an arrow.
        """
        position, direction = angles_to_vectors(
            data.r, data.phi, data.z, data.poloidal_angle, data.toroidal_angle
        )
        return create_vtk_arrows(
            position, direction, scaling_factor=self.b_probe_length
        )

    def _create_loop(self, data):
        """Create a closed VTK polyline from loop positions.

        Args:
            data: FluxLoop or RogowskiCoil dataclass containing positions.

        Returns:
            VTK polydata representing the closed loop.
        """
        positions = data.positions
        points = []
        for pos in positions:
            x, y = pol_to_cart(pos.r, pos.phi)
            points.append((x, y, pos.z))

        return points_to_vtkpoly(points, is_closed=True, is_filled=False)

    def _create_name(self, device_type, device, index):
        """Create a unique name for a diagnostic device.

        Args:
            device_type: Type of device.
            device: Device IDS node.
            index: Index of the device.

        Returns:
            A unique name string for the device.
        """
        name = (
            ensure_unique_name(
                f"{device_type} ({device.name})", list(self.selectable_map.keys())
            )
            if device.name
            else f"{device_type} #{index}"
        )
        return name
