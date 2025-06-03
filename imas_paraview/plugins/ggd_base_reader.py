import logging

from vtkmodules.vtkCommonDataModel import vtkDataObject, vtkPartitionedDataSetCollection

from imas_paraview.convert import InterpSettings
from imas_paraview.io import read_ps
from imas_paraview.plugins.base_class import GGDVTKPluginBase
from imas_paraview.progress import Progress

logger = logging.getLogger("imas_paraview")


class GGDBaseReader(GGDVTKPluginBase):
    """Base class for GGDReader and the JorekGGDReader."""

    def __init__(self, supported_ids):
        super().__init__("vtkPartitionedDataSetCollection", supported_ids)
        # GGD arrays to load
        self._all_vector_paths = []
        self._all_scalar_paths = []

    def GetAttributeArrayName(self, idx) -> str:
        return self._name_from_idspath(self._selectable[idx])

    def RequestDataObject(self, request, inInfo, outInfo):
        output = vtkPartitionedDataSetCollection()
        outInfo.GetInformationObject(0).Set(vtkDataObject.DATA_OBJECT(), output)
        outInfo.GetInformationObject(0).Set(
            vtkDataObject.DATA_EXTENT_TYPE(), output.GetExtentType()
        )
        return 1

    def RequestData(self, request, inInfo, outInfo):
        if self._dbentry is None or not self._ids_and_occurrence or self._ids is None:
            return 1

        # Retrieve the selected time step and GGD arrays
        selected_scalar_paths, selected_vector_paths = self._get_selected_ggd_paths()
        time = self._get_selected_time_step(outInfo)
        if time is None:
            logger.warning("Selected invalid time step")
            return 1

        # Create progress object to advance Paraview progress bar
        progress = Progress(self.UpdateProgress)

        if self._n_plane < self._MIN_N_PLANE:
            raise ValueError(
                f"The number of Bezier planes cannot be less than {self._MIN_N_PLANE}."
            )

        plane_config = InterpSettings(
            n_plane=self._n_plane, phi_start=self._phi_start, phi_end=self._phi_end
        )
        output = self.converter.ggd_to_vtk(
            time=time,
            scalar_paths=selected_scalar_paths,
            vector_paths=selected_vector_paths,
            plane_config=plane_config,
            outInfo=outInfo,
            progress=progress,
        )

        if output is None:
            logger.warning("Could not convert GGD to VTK.")
        return 1

    def request_information(self):
        """
        Select which GGD arrays to show in the array domain selector, based
        on whether the "Show All" checkbox is enabled.
        """
        if self._ids is not None:
            if self.show_all:
                self._selectable = self._all_vector_paths + self._all_scalar_paths
            else:
                self._selectable = self._filled_vector_paths + self._filled_scalar_paths
        pass

    def setup_ids(self):
        """
        Initializes plasma state reader and populates scalar and vector paths for
        the array selection domaindomain selection.
        """
        if self._ids is not None:
            # Load paths from IDS
            ps_reader = read_ps.PlasmaStateReader(self._ids)
            (
                self._all_scalar_paths,
                self._all_vector_paths,
                self._filled_scalar_paths,
                self._filled_vector_paths,
            ) = ps_reader.load_paths_from_ids()

    def _name_from_idspath(self, path):
        """Converts an IDSPath to a string by removing 'ggd' and capitalizing each part
        of the path, with the parts separated by spaces.

        Example:
            If path is IDSPath('ggd/electrons/pressure'), the function returns
            "Electrons Pressure"

        Args:
            path: The IDSPath object to convert into a formatted string

        Returns:
            A formatted string of the IDSPath
        """
        path_list = list(path.parts)
        if "ggd" in path_list:
            path_list.remove("ggd")
        if "description_ggd" in path_list:
            path_list.remove("description_ggd")
        for i in range(len(path_list)):
            path_list[i] = path_list[i].capitalize()

        name = " ".join(path_list)

        # If GGDs are not filled in the first time step, add (?) to their name. This
        # notifies that the GGD array is not filled and will most likely not contain
        # any data. We do not remove these GGD arrays from the selector window entirely,
        # because later time steps could still contain data.
        if (
            path not in self._filled_scalar_paths
            and path not in self._filled_vector_paths
        ):
            name = f"{name} (?)"
        return name

    def _get_selected_ggd_paths(self):
        """Retrieve the IDSPaths of the selected scalar and vector GGD arrays.

        Returns:
            selected_scalar_paths: List of IDSPaths of selected scalar GGD arrays
            selected_vector_paths: List of IDSPaths of selected vector GGD arrays
        """

        # Determine if selected GGD arrays are scalar or vector arrays
        selected_scalar_paths = [
            obj
            for obj in self._all_scalar_paths
            if self._name_from_idspath(obj) in self._selected
        ]
        selected_vector_paths = [
            obj
            for obj in self._all_vector_paths
            if self._name_from_idspath(obj) in self._selected
        ]
        return selected_scalar_paths, selected_vector_paths
