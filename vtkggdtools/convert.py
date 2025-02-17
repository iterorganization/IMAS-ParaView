import logging
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from vtk import vtkXMLPartitionedDataSetCollectionWriter
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import (
    vtkCompositeDataSet,
    vtkDataAssembly,
    vtkPartitionedDataSetCollection,
)

from vtkggdtools.io import read_geom, read_jorek, read_ps
from vtkggdtools.util import find_closest_indices, get_grid_ggd

logger = logging.getLogger("vtkggdtools")


@dataclass
class InterpSettings:
    """Data class containing Fourier interpolation settings."""

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
        self.get_grids_at_time = lru_cache(maxsize=32)(self.get_grids_at_time)

    def write_to_xml(self, output_path: Path, index_list=[0]):
        """Convert an IDS to VTK format and write it to disk using the XML output
        writer.

        Args:
            output_path: Path for the output directory.
            index_list: A list of time indices to convert. By default only the first
                time step is converted.
        """
        logger.info(f"Creating a output directory at {output_path}")
        output_path.mkdir(parents=True, exist_ok=True)
        output_path = output_path / self.ids.metadata.name
        if any(index >= len(self.ids.time) for index in index_list):
            raise RuntimeError("A provided index is out of bounds.")

        for index in index_list:
            logger.info(f"Converting time step {self.ids.time[index]}...")
            vtk_object = self.ggd_to_vtk(time_idx=index)
            if vtk_object is None:
                logger.warning(f"Could not convert GGD at time index {index} to VTK.")
                continue
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
            outInfo: Paraview's Source outInfo information object.
            progress: Progress indicator for Paraview.

        Returns:
            vtkPartitionedDataSetCollection containing the converted GGD data.
        """
        self.points = vtkPoints()
        self.assembly = vtkDataAssembly()
        self.time_idx = self._resolve_time_idx(time_idx, time)

        self.progress = progress

        if self.time_idx is None:
            return None

        self.grid_ggd = get_grid_ggd(self.ids, self.time_idx)

        if not self._is_grid_valid():
            return None

        if outInfo is None:
            self.output = vtkPartitionedDataSetCollection()
        else:
            self.output = vtkPartitionedDataSetCollection.GetData(outInfo)
        self.output.SetDataAssembly(self.assembly)

        self.ps_reader = read_ps.PlasmaStateReader(self.ids)
        self.ps_reader.load_arrays_from_path(self.time_idx, scalar_paths, vector_paths)
        if plane_config.n_plane != 0:
            self._interpolate_jorek(plane_config)
        else:
            self._fill_grid_and_plasma_state()

        return self.output

    def _resolve_time_idx(self, time_idx, time):
        """Resolves the appropriate time index based on the given time index or time
        value.

        Args:
            time_idx: The index of the time step to convert.
            time: The time value to convert.

        Returns:
            The resolved time index, or None if an error occurs.
        """
        if time is not None and time_idx is not None:
            logger.error("The time and time index cannot be provided at the same time.")
            return None
        elif time_idx is not None:
            if time_idx >= len(self.ids.time):
                logger.error("The requested index cannot be found in the IDS.")
                return None
            return time_idx
        elif time is not None:
            indices = find_closest_indices([time], self.ids.time)
            if len(indices) == 0:
                logger.warning(
                    "No time steps found that are less than or equal to the provided "
                    "time. Converting the first time step instead."
                )
                time_idx = 0
            else:
                time_idx = indices[0]
            logger.info(
                f"Converting timestep: t = {self.ids.time[time_idx]} at index = "
                f"{time_idx}"
            )
            return time_idx
        else:
            time_idx = len(self.ids.time) // 2
            logger.info(
                "No time or time index provided, so converting the middle time "
                f"step: t = {self.ids.time[time_idx]} at index {time_idx}."
            )
            return time_idx

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
        n_period = self.grid_ggd.space[1].geometry_type.index
        if n_period > 0:
            ugrid = read_jorek.convert_grid_subset_to_unstructured_grid(
                self.ids, self.time_idx, plane_config, self.ps_reader
            )
            self.output.SetPartition(0, 0, ugrid)
            child = self.assembly.AddNode(self.ids.metadata.name, 0)
            self.assembly.AddDataSetIndex(child, 0)
            self.output.GetMetaData(0).Set(
                vtkCompositeDataSet.NAME(), self.ids.metadata.name
            )
        else:
            logger.error("Invalid plane configuration for the given IDS type.")

    def _fill_grid_and_plasma_state(self):
        """Fills the VTK output object with the GGD grid and GGD array values."""
        num_subsets = len(self.grid_ggd.grid_subset)

        ugrids = self.get_grids_at_time(self.time_idx)
        if num_subsets <= 1:
            logger.info("No subsets to read from grid_ggd")
            self.output.SetNumberOfPartitionedDataSets(1)
            self._set_partition(0, ugrids[-1], -1)
            self.ps_reader.read_plasma_state(-1, ugrids[-1])
        elif self.ids.metadata.name == "wall":
            # FIXME: what if num_subsets is 2 or 3?
            self.output.SetNumberOfPartitionedDataSets(num_subsets - 3)
            self.partition = 0
            self._set_partition(0, ugrids[-1], -1)
            self.ps_reader.read_plasma_state(-1, ugrids[-1])

            for subset_idx in range(4, num_subsets):
                self._set_partition(subset_idx - 3, ugrids[subset_idx], subset_idx)
                self.ps_reader.read_plasma_state(subset_idx, ugrids[subset_idx])
                if self.progress:
                    self.progress.increment(0.5 / num_subsets)
        else:
            self.output.SetNumberOfPartitionedDataSets(num_subsets)
            for subset_idx in range(num_subsets):
                self._set_partition(subset_idx, ugrids[subset_idx], subset_idx)
                self.ps_reader.read_plasma_state(subset_idx, ugrids[subset_idx])
                if self.progress:
                    self.progress.increment(0.5 / num_subsets)

    def get_grids_at_time(self, time_idx):
        """Fetches the unstructured grids at a certain time index. Note that this
        function is cached using lru_cache based on `time_idx`.

        Args:
            time_idx: Time index of the grid to get, is used for creation the hash of
                the cache entry.
        """
        logger.info("No cache found, loading the grid from the IDS.")
        num_subsets = len(self.grid_ggd.grid_subset)
        read_geom.fill_vtk_points(
            self.grid_ggd, 0, self.points, self.ids.metadata.name, self.progress
        )
        ugrids = {}
        for subset_idx in range(-1, num_subsets):
            ugrids[subset_idx] = (
                read_geom.convert_grid_subset_geometry_to_unstructured_grid(
                    self.grid_ggd, subset_idx, self.points, self.progress
                )
            )
        return ugrids

    def _set_partition(self, partition, ugrid, subset_idx):
        """Sets a partition in the output dataset and updates the assembly structure.

        Args:
            partition: Partition index in the output dataset.
            ugrid: Unstructured grid data.
            subset_idx: Index of the subset in `self.grid_ggd.grid_subset`, if negative
                the metadata name is used.
        """
        subset = None if subset_idx < 0 else self.grid_ggd.grid_subset[subset_idx]
        label = str(subset.identifier.name) if subset else self.ids.metadata.name

        self.output.SetPartition(partition, 0, ugrid)
        child = self.assembly.AddNode(label.replace(" ", "_"), 0)
        self.assembly.AddDataSetIndex(child, partition)

        self.output.GetMetaData(partition).Set(vtkCompositeDataSet.NAME(), label)

    def _write_vtk_to_xml(self, vtk_object, output_path):
        """Writes a VTK object to XML file, and a directory containing files for each
        grid subset. The directory structure looks as follows:
        .
        ├── test
        |   ├── test_0_0.vtu
        |   ├── test_1_0.vtu
        |   ├── test_2_0.vtu
        |   ├── ...
        └── test.vtpc

        Args:
            vtk_object: The VTK object to be written to disk.
            output_path: The output file path where the VTK objects will be saved.
        """
        if vtk_object is None:
            logger.error("Cannot write None object to XML.")
            return

        logger.info(f"Writing VTK file to {output_path}...")
        writer = vtkXMLPartitionedDataSetCollectionWriter()
        writer.SetInputData(vtk_object)
        writer.SetFileName(output_path.with_suffix(".vtpc"))
        writer.Write()
        logger.info(f"Successfully wrote VTK object to {output_path}")
