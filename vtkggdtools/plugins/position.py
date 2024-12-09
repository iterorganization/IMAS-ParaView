"""IMASPy plugin to view beam structures of ec_launchers IDS in Paraview."""

import logging
from dataclasses import dataclass

import numpy as np
from imaspy.ids_structure import IDSStructure
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.vtkCommonCore import vtkPoints

from vtkggdtools.ids_util import get_object_by_name
from vtkggdtools.paraview_support.servermanager_tools import doublevector, propertygroup
from vtkggdtools.plugins.base_class import GGDVTKPluginBase
from vtkggdtools.util import find_closest_indices, points_to_vtkpoly, pol_to_cart

logger = logging.getLogger("vtkggdtools")

SUPPORTED_IDS_NAMES = ["barometry", "langmuir_probes"]


@dataclass
class PositionStructure:
    """Data class that stores ``beam`` IDS structures, along with its name."""

    name: str
    position_structure: IDSStructure


@smproxy.source(label="Beam Reader")
@smhint.xml("""<ShowInMenu category="VTKGGDTools" />""")
class IMASPyBeamReader(GGDVTKPluginBase):
    """Beam reader based on IMASPy"""

    def __init__(self):
        super().__init__("vtkPoints", ["ec_launchers"])
        self.distance = 10

    @doublevector(label="Beam Distance", name="beam_distance", default_values=10)
    def P99_SetDistance(self, val):
        """Sets the beam's distance from the launching position."""
        self._update_property("distance", val)

    @propertygroup("Beam settings", ["beam_distance"])
    def PG3_BeamGroup(self):
        """Dummy function to define a PropertyGroup."""

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

        if len(self._selected) > 0:
            output = vtkPoints.GetData(outInfo)
            self._load_position(output)
        return 1

    def request_information(self):
        """
        Abstract method which is called during the RequestInformation pipeline step.
        Intentionally left empty.
        """

    def setup_ids(self):
        """
        Select which beams to show in the array domain selector. The IDS is is searched
        for beam structures and their names are added to the array domain selector.
        """
        assert self._ids is not None, "IDS cannot be empty during setup."

        self._selectable = []

        if self._ids.metadata.name == "barometry":
            aos = self._ids.embedded
        elif self._ids.metadata.name == "langmuir_probes":
            aos = self._ids.gauge
        else:
            raise NotImplementedError(
                f"Unable to find load position for {self._ids.metadata.name}."
            )

        for i, structure in enumerate(aos):
            name = structure.name
            if name == "":
                name = f" {i}"
                logger.warning(
                    "Found a channel without a name, " f"it will be loaded as {name}."
                )

            selectable = PositionStructure(str(name), structure)
            self._selectable.append(selectable)

    def _load_position(self, output):
        """Go through the list of selected beams, and load each of them in a
        separate vtkPolyData object, which are all combined into a vtkMultiBlockDataSet.

        Args:
            output: The vtkMultiBlockDataSet containing the beams.
        """
        for name in self._selected:
            pos_struct = get_object_by_name(self._selectable, name)
            if pos_struct is None:
                continue

            pos = pos_struct.position
            if pos_struct is None:
                raise ValueError(f"Could not find {name}")

            logger.info(f"Selected {pos_struct.name}")
            pos_cart = (*pol_to_cart(pos.r, pos.phi), pos.z)
            output.InsertNextPoint(*pos_cart)
