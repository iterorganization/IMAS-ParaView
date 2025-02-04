"""IMASPy plugin to view positional data of barometry gauges and embedded langmuir
probes."""

import logging
from dataclasses import dataclass

from imaspy.ids_structure import IDSStructArray, IDSStructure
from paraview.util.vtkAlgorithm import smhint, smproxy
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkPolyData

from vtkggdtools.ids_util import get_object_by_name
from vtkggdtools.plugins.base_class import GGDVTKPluginBase
from vtkggdtools.util import pol_to_cart

logger = logging.getLogger("vtkggdtools")

SUPPORTED_IDS_NAMES = ["barometry", "langmuir_probes", "magnetics"]


@dataclass
class PositionStructure:
    """Data class that stores positional IDS structures, along with its name. For
    barometry the structures are gauges, for langmuire probes, the structures are
    embedded probes"""

    name: str
    position_structure: IDSStructure


@smproxy.source(label="Position Reader")
@smhint.xml("""<ShowInMenu category="IMAS Tools" />""")
class IMASPyPositionReader(GGDVTKPluginBase):
    """Positions reader based on IMASPy"""

    def __init__(self):
        super().__init__("vtkPolyData", SUPPORTED_IDS_NAMES)

    def GetAttributeArrayName(self, idx) -> str:
        return self._selectable[idx].name

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        if len(self._selected) > 0:
            output = vtkPolyData.GetData(outInfo)
            self._load_position(output)
        return 1

    def setup_ids(self):
        """
        Select which probes or gauges to show in the array domain selector. The names
        from the gauges or probes are to the array domain selector.
        """
        assert self._ids is not None, "IDS cannot be empty during setup."

        self._selectable = []

        if self._ids.metadata.name == "barometry":
            aos_list = [self._ids.gauge]
        elif self._ids.metadata.name == "langmuir_probes":
            aos_list = [self._ids.embedded]
        elif self._ids.metadata.name == "magnetics":
            aos_list = [
                self._ids.flux_loop,
                self._ids.bpol_probe,
                self._ids.b_field_pol_probe,
                self._ids.b_field_tor_probe,
                self._ids.rogowski_coil,
                self._ids.b_field_phi_probe,
            ]

        else:
            raise NotImplementedError(
                f"Unable to find load position for {self._ids.metadata.name}."
            )
        for aos in aos_list:
            for i, structure in enumerate(aos):
                name = structure.name
                if name == "":
                    name = f"device {i}"
                    logger.warning(
                        f"Found a device without a name, it will be loaded as {name}."
                    )

                if hasattr(structure, "identifier"):
                    identifier = structure.identifier
                    if not identifier == "":
                        name = f"{name} / {identifier}"

                selectable = PositionStructure(name, structure)
                self._selectable.append(selectable)

    def _load_position(self, output):
        """Go through the list of selected position structures, and load each of them
        in a vtkPoints object, which are all saved into a vtkPolyData object.

        Args:
            output: The vtkPolyData containing the positions.
        """
        points = vtkPoints()
        for name in self._selected:
            pos_struct = get_object_by_name(self._selectable, name)

            if pos_struct is None:
                raise ValueError(f"Could not find {name}")

            pos = pos_struct.position_structure.position
            logger.info(f"Selected {pos_struct.name}")

            if isinstance(pos, IDSStructArray):
                for pos_struct in pos:
                    pos_cart = (
                        *pol_to_cart(pos_struct.r, pos_struct.phi),
                        pos_struct.z,
                    )
                    points.InsertNextPoint(*pos_cart)
            else:
                pos_cart = (*pol_to_cart(pos.r, pos.phi), pos.z)
                points.InsertNextPoint(*pos_cart)

        output.SetPoints(points)
