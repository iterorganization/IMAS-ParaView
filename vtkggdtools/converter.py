import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from vtk import vtkXMLPartitionedDataSetCollectionWriter
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import (
    vtkCompositeDataSet,
    vtkDataAssembly,
    vtkPartitionedDataSetCollection,
)

from vtkggdtools.io import read_bezier, read_geom, read_ps
from vtkggdtools.util import FauxIndexMap, get_grid_ggd

logger = logging.getLogger("vtkggdtools")


@dataclass
class InterpSettings:
    """Helper data class containing Bezier interpolation settings."""

    n_plane: int = 0
    phi_start: float = 0.0
    phi_end: float = 0.0


class Converter:
    def __init__(self, ids):
        self.ids = ids
        self.time_idx = None
        self.grid_ggd = None
        self.output = None
        self.ps_reader = None
        self.ugrids = []

    def write_to_xml(self, output_path: Path, index_list=[0]):
        """Convert an IDS to VTK format and write it to disk using the XML output
        writer.

        Args:
        ids: The IDS to be converted.
        output: The name of the output directory.
        index_list: A list of time indices to convert. By default only the first time
            step is converted.
        """
        logger.info(f"Creating a output directory at {output_path}")
        output_path.mkdir(parents=True, exist_ok=True)
        output_path = output_path / self.ids.metadata.name
        if any(index >= len(self.ids.time) for index in index_list):
            raise RuntimeError("A provided index is out of bounds.")

        for index in index_list:
            logger.info(f"Converting time step {self.ids.time[index]}...")
            vtk_object = self.ggd_to_vtk(time_idx=index)
            self._write_vtk_to_xml(vtk_object, Path(f"{output_path}_{index}"))

    def ggd_to_vtk(
        self,
        time=None,
        time_idx=None,
        scalar_paths=None,
        vector_paths=None,
        plane_config: InterpSettings = InterpSettings(),
        outInfo=None,
        progress=None,
        ugrids=None,
    ):
        """Converts the GGD of an IDS to VTK format.

        Args:

            time: Time step to convert. Defaults to converting the first time step.
            time_idx: Time index to convert. Defaults to converting the first time step.
            scalar_paths: A list of IDSPaths of GGD scalar arrays to convert. Defaults
                to None, in which case all scalar arrays are converted.
            vector_paths: A list of IDSPaths of GGD vector arrays to convert. Defaults
                to None, in which case all vectors arrays are converted.
            plane_config: Data class containing the interpolation settings.
            outInfo: Source outInfo information object. Defaults to None.
            progress: Progress indicator for Paraview. Defaults to None.

        Returns:
            vtkPartitionedDataSetCollection containing the converted GGD data.
        """
        self.points = vtkPoints()
        self.assembly = vtkDataAssembly()
        self.time_idx = self._resolve_time_idx(time_idx, time)
        self.input_ugrids = ugrids
        if self.time_idx is None:
            return None, None

        self.grid_ggd = get_grid_ggd(self.ids, self.time_idx)

        if not self._is_grid_valid():
            return None, None

        self._setup_vtk_context(outInfo)

        if plane_config.n_plane != 0:
            self._interpolate_jorek(plane_config)
            return self.output, None

        self.ps_reader = read_ps.PlasmaStateReader(self.ids)
        self.ps_reader.load_arrays_from_path(self.time_idx, scalar_paths, vector_paths)
        self._fill_grid_and_plasma_state(progress)

        return self.output, self.ugrids

    def _resolve_time_idx(self, time_idx, time):
        if time is not None and time_idx is not None:
            logger.error("The time and time index cannot be provided at the same time.")
            return None
        elif time_idx is not None:
            if time_idx >= len(self.ids.time):
                logger.error("The requested index cannot be found in the IDS.")
                return None
            return time_idx
        elif time is not None:
            return self._get_nearest_time_idx(time)
        else:
            time_idx = len(self.ids.time) // 2
            logger.info(
                "No time or time index provided, so converting the middle time "
                f"step: t = {self.ids.time[time_idx]} at index {time_idx}."
            )
            return time_idx

    def _get_nearest_time_idx(self, time):
        """Finds the index of the nearest time step in the IDS time array that is less
        than or equal to the provided time value.

        Args:
            time: Timestep to retrieve. If it is None, the first time step is retrieved.

        Returns:
            Index of the nearest time step
        """
        candidates = self.ids.time[self.ids.time <= time]
        if candidates.size == 0:
            logger.warning(
                "No time steps found that are less than or equal to the provided time."
                " Converting the first time step instead."
            )
            return 0

        nearest_time = candidates.max()
        time_idx = np.where(self.ids.time == nearest_time)[0][0]
        logger.info(
            f"Converting timestep: t = {self.ids.time[time_idx]} at index = {time_idx}"
        )
        return time_idx

    def _setup_vtk_context(self, outInfo):
        """Setup the VTK context with output and assembly."""
        if outInfo is None:
            self.output = vtkPartitionedDataSetCollection()
        else:
            self.output = vtkPartitionedDataSetCollection.GetData(outInfo)

        read_geom.fill_vtk_points(self.grid_ggd, 0, self.points, self.ids.metadata.name)
        self.output.SetDataAssembly(self.assembly)

    def _is_grid_valid(self):
        """Validates if the grid is properly loaded."""
        if self.grid_ggd is None:
            logger.warning("Could not load a valid GGD grid.")
            return False
        if not hasattr(self.grid_ggd, "space") or len(self.grid_ggd.space) < 1:
            logger.warning("The grid_ggd does not contain a space.")
            return False
        return True

    def _interpolate_jorek(self, plane_config: InterpSettings):
        """Interpolate JOREK Fourier space.

        Args:
            plane_config: Data class containing the interpolation settings.
        """
        aos_index_values = FauxIndexMap()
        n_period = self.grid_ggd.space[1].geometry_type.index
        if n_period > 0:
            ugrid = read_bezier.convert_grid_subset_to_unstructured_grid(
                self.ids.metadata.name,
                self.ids,
                aos_index_values,
                plane_config.n_plane,
                plane_config.phi_start,
                plane_config.phi_end,
            )
            self.output.SetPartition(0, 0, ugrid)
            child = self.assembly.AddNode(self.ids.metadata.name, 0)
            self.assembly.AddDataSetIndex(child, 0)
            self.output.GetMetaData(0).Set(
                vtkCompositeDataSet.NAME(), self.ids.metadata.name
            )
        else:
            logger.error("Invalid plane configuration for the given IDS type.")

    def _fill_grid_and_plasma_state(self, progress):
        """Fill grid and plasma state data."""
        num_subsets = len(self.grid_ggd.grid_subset)

        if num_subsets <= 1:
            logger.info("No subsets to read from grid_ggd")
            self.output.SetNumberOfPartitionedDataSets(1)
            ugrid = self._get_ugrid(-1, 0)
            self.ps_reader.read_plasma_state(-1, ugrid)
        elif self.ids.metadata.name == "wall":
            # FIXME: what if num_subsets is 2 or 3?
            self.output.SetNumberOfPartitionedDataSets(num_subsets - 3)
            ugrid = self._get_ugrid(-1, 0)
            self.ps_reader.read_plasma_state(-1, ugrid)

            for subset_idx in range(4, num_subsets):
                ugrid = self._get_ugrid(subset_idx, subset_idx - 3)
                self.ps_reader.read_plasma_state(subset_idx, ugrid)

                if progress:
                    progress.increment(1.0 / num_subsets)
        else:
            self.output.SetNumberOfPartitionedDataSets(num_subsets)
            for subset_idx in range(num_subsets):
                ugrid = self._get_ugrid(subset_idx, subset_idx)
                self.ps_reader.read_plasma_state(subset_idx, ugrid)
                if progress:
                    progress.increment(1.0 / num_subsets)

    def _get_ugrid(self, subset_idx, partition):
        if self.input_ugrids is None:
            ugrid = self._fill_grid(subset_idx, partition)
        else:
            ugrid = self._fill_grid(
                subset_idx, partition, ugrid=self.input_ugrids[subset_idx]
            )
        self.ugrids.append(ugrid)
        return ugrid

    def _fill_grid(self, subset_idx, partition, ugrid=None):
        """Read GGD data from the IDS and convert it to VTK data."""
        subset = None if subset_idx < 0 else self.grid_ggd.grid_subset[subset_idx]
        if ugrid is None:
            ugrid = read_geom.convert_grid_subset_geometry_to_unstructured_grid(
                self.grid_ggd, subset_idx, self.points
            )
        self.output.SetPartition(partition, 0, ugrid)
        label = str(subset.identifier.name) if subset else self.ids.metadata.name
        child = self.assembly.AddNode(label.replace(" ", "_"), 0)
        self.assembly.AddDataSetIndex(child, partition)
        self.output.GetMetaData(partition).Set(vtkCompositeDataSet.NAME(), label)
        return ugrid

    def _write_vtk_to_xml(self, vtk_object, output_file):
        """Writes the VTK object to disk."""
        if vtk_object is None:
            logger.error("Could not convert GGD to VTK file.")
            raise RuntimeError

        logger.info(f"Writing VTK file to {output_file}...")
        writer = vtkXMLPartitionedDataSetCollectionWriter()
        writer.SetInputData(vtk_object)
        writer.SetFileName(output_file.with_suffix(".vtpc"))
        writer.Write()
        logger.info(f"Successfully wrote VTK object to {output_file}")
