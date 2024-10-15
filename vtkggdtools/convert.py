import logging
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


def convert_to_xml(
    ids,
    output,
    index=None,
    index_range=None,
    time=None,
    time_range=None,
    all_times=False,
):
    """_summary_

    Args:
        ids: _description_
        output: _description_
        time_options: _description_
    """
    # Check if at most one of the time options are provided
    if (
        sum(param is not None for param in [index, index_range, time, time_range])
        + all_times
        > 1
    ):
        raise RuntimeError(
            "You can only provide one of the following parameters: index, index_range, "
            "time, time_range or all_times."
        )
    vtk_object = None

    # Create parent directory and point output path there
    logger.info(f"Creating a output directory at {output}")
    output.mkdir(parents=True, exist_ok=True)
    output = output / output.stem

    # Convert a single time step
    if time is not None:
        logger.info(f"Converting time step closest to {time}.")
        vtk_object = ggd_to_vtk(ids, time=time)
        _write_vtk(vtk_object, output)
    elif index is not None:
        logger.info(f"Converting time step at index {index}.")
        vtk_object = ggd_to_vtk(ids, time_idx=index)
        _write_vtk(vtk_object, output)

    # Convert a slice of time steps
    elif time_range:
        min_idx = np.argmin(np.abs(ids.time - time_range[0]))
        max_idx = np.argmin(np.abs(ids.time - time_range[1]))
        sliced_indices = range(min_idx, max_idx)
        for idx in sliced_indices:
            vtk_object = ggd_to_vtk(ids, time_idx=idx)
            _write_vtk(vtk_object, Path(f"{output}_{ids.time[idx]}"))
    elif index_range:
        if index_range[1] >= len(ids.time):
            raise RuntimeError(
                f"Index {index_range[1]} is out of bounds of the IDS time array "
                f"of size {len(ids.time)}"
            )
        sliced_indices = range(index_range[0], index_range[1] + 1)
        for idx in sliced_indices:
            vtk_object = ggd_to_vtk(ids, time_idx=idx)
            _write_vtk(vtk_object, Path(f"{output}_{idx}"))

    # Convert all time steps
    elif all_times:
        for idx in range(len(ids.time)):
            vtk_object = ggd_to_vtk(ids, time_idx=idx)
            _write_vtk(vtk_object, Path(f"{output}_{idx}"))
    else:
        raise RuntimeError("No time was selected to convert. Aborting.")


def ggd_to_vtk(
    ids,
    *,
    time=None,
    time_idx=None,
    scalar_paths=None,
    vector_paths=None,
    n_plane=0,
    phi_start=0,
    phi_end=0,
    outInfo=None,
    progress=None,
):
    """Converts the GGD of an IDS to VTK format.

    Args:
        ids: The IDS to convert to VTK.
        time: Time step to convert. Defaults to converting the first time step.
        scalar_paths: A list of IDSPaths of GGD scalar arrays to convert. Defaults
            to None, in which case all scalar arrays are converted.
        vector_paths: A list of IDSPaths of GGD vector arrays to convert. Defaults
            to None, in which case all vectors arrays are converted.
        n_plane: Number of toroidal planes to be generated if 3D axysimetric. Defaults
            to 0.
        phi_start: Start phi plane in degrees. Defaults to 0.
        phi_end: End plane at phi in degrees. Defaults to 0.
        outInfo: Source outInfo information object. Defaults to None.
        progress: Progress indicator for Paraview. Defaults to None.

    Returns:
        vtkPartitionedDataSetCollection containing the converted GGD data.
    """
    if time and time_idx:
        logger.error("The time and time index can not be provided at the same time.")
        return None
    if time_idx and time_idx >= len(ids.time):
        logger.error("The requested index can not be found in the IDS.")
        return None
    # Retrieve GGD grid from IDS
    if time:
        time_idx = _get_nearest_time_idx(ids, time)
    elif not time and not time_idx:
        logger.info(f"Converting first timestep: t = {ids.time[time_idx]}")
        time_idx = 0

    grid_ggd = get_grid_ggd(ids, time_idx)

    # Check if grid is valid
    if grid_ggd is None:
        logger.warning("Could not load a valid GGD grid.")
        return None
    if not hasattr(grid_ggd, "space") or len(grid_ggd.space) < 1:
        logger.warning("The grid_ggd does not contain a space.")
        return None
    num_subsets = len(grid_ggd.grid_subset)

    # Create output VTK object
    if outInfo is None:
        output = vtkPartitionedDataSetCollection()
    else:
        output = vtkPartitionedDataSetCollection.GetData(outInfo)

    points = vtkPoints()
    space_idx = 0
    ids_name = ids.metadata.name

    read_geom.fill_vtk_points(grid_ggd, space_idx, points, ids_name)
    assembly = vtkDataAssembly()
    output.SetDataAssembly(assembly)

    if n_plane != 0:
        _interpolate_jorek(ids, grid_ggd, n_plane, phi_start, phi_end, output, assembly)
        return output

    # Load the GGD arrays from the selected GGD paths
    ps_reader = read_ps.PlasmaStateReader(ids)
    ps_reader.load_arrays_from_path(time_idx, scalar_paths, vector_paths)

    if num_subsets <= 1:
        logger.info("No subsets to read from grid_ggd")
        output.SetNumberOfPartitionedDataSets(1)
        _fill_grid_and_plasma_state(
            ids_name, grid_ggd, ps_reader, -1, 0, points, output, assembly
        )
    elif ids_name == "wall":
        # FIXME: what if num_subsets is 2 or 3?
        output.SetNumberOfPartitionedDataSets(num_subsets - 3)
        _fill_grid_and_plasma_state(
            ids_name, grid_ggd, ps_reader, -1, 0, points, output, assembly
        )
        for subset_idx in range(4, num_subsets):
            _fill_grid_and_plasma_state(
                ids_name,
                grid_ggd,
                ps_reader,
                subset_idx,
                subset_idx - 3,
                points,
                output,
                assembly,
            )
            if progress:
                progress.increment(1.0 / num_subsets)
    else:
        output.SetNumberOfPartitionedDataSets(num_subsets)
        for subset_idx in range(num_subsets):
            _fill_grid_and_plasma_state(
                ids_name,
                grid_ggd,
                ps_reader,
                subset_idx,
                subset_idx,
                points,
                output,
                assembly,
            )
            if progress:
                progress.increment(1.0 / num_subsets)

    logger.info("Finished loading IDS.")
    return output


def _write_vtk(vtk_object, output):
    """_summary_

    Args:
        vtk_object: _description_
        output: _description_

    Raises:
        RuntimeError: _description_
    """
    if vtk_object is None:
        logger.info("Could not convert GGD to VTK file.")
        raise RuntimeError
    logger.info("Writing VTK file to disk...")
    writer = vtkXMLPartitionedDataSetCollectionWriter()
    writer.SetInputData(vtk_object)
    output_file = output.with_suffix(".vtpc")
    writer.SetFileName(output_file)
    writer.Write()
    logger.info(f"Successfully wrote VTK object to {output_file}.")


def _get_nearest_time_idx(ids, time):
    """Finds the index of the nearest time step in the IDS time array to the provided
    time value.

    Args:
        ds: The IDS object.
        time: Timestep to retrieve. If it is None, the first time step is retrieved.

    Returns:
        Index of the nearest time step
    """
    time_idx = np.argmin(np.abs(ids.time - time))
    logger.info(f"Converting timestep: t = {ids.time[time_idx]} at index = {time_idx}")
    return time_idx


def _interpolate_jorek(ids, grid_ggd, n_plane, phi_start, phi_end, output, assembly):
    """Interpolate JOREK Fourier space.

    Args:
        ids: The IDS object.
        grid_ggd: The grid_ggd IDS node.
        n_plane: Number of toroidal planes to be generated if 3D axysimetric.
        phi_start: Start phi plane.
        phi_end: End plane at phi in degrees.
        output: vtkPartitionedDataSetCollection containing the converted GGD data.
        assembly: vtkDataAssembly containing the hierarchical hierarchical organization
            of items in the vtkPartitionedDataSetCollection.
    """
    # TODO: allow selecting other grids for Bezier
    aos_index_values = FauxIndexMap()
    # Interpolate JOREK Fourier space
    # TODO: figure out if we can put this functionality in a post-processing step?
    ids_name = ids.metadata.name
    number_of_spaces = len(grid_ggd.space)
    if number_of_spaces > 1 and len(grid_ggd.space[0].coordinates_type) == 2:
        n_period = grid_ggd.space[1].geometry_type.index
        if n_period > 0:  # Fourier space with periodicity (JOREK)
            logger.info(f"Reading Bezier mesh with Fourier periodiciy {n_period}")
            ugrid = read_bezier.convert_grid_subset_to_unstructured_grid(
                ids_name,
                ids,
                aos_index_values,
                n_plane,
                phi_start,
                phi_end,
            )
            output.SetPartition(0, 0, ugrid)
            child = assembly.AddNode(ids_name, 0)
            assembly.AddDataSetIndex(child, 0)
            output.GetMetaData(0).Set(vtkCompositeDataSet.NAME(), ids_name)
        else:
            logger.error(
                f"Number of planes {n_plane} invalid for this "
                f"{number_of_spaces} number of spaces"
            )
            output = None
    else:
        logger.error(
            f"Number of planes {n_plane} invalid for this IDS type." " Try using N = 0"
        )
        output = None


def _fill_grid_and_plasma_state(
    ids_name,
    grid_ggd,
    ps_reader,
    subset_idx,
    partition,
    vtk_grid_points,
    output,
    assembly,
):
    """Read GGD data from the IDS and convert it to VTK data.

    Args:
        ids_name: Name of the IDS object.
        grid_ggd: The grid_ggd IDS node.
        ps_reader: Plasmastate reader class responsible for reading the IDS data.
        subset_idx: Index of the grid subset.
        partition: Index of the VTK partition.
        vtk_grid_points: The point coordinates corresponding to 1d objects in
            the subset elements.
        output: vtkPartitionedDataSetCollection containing the converted GGD data.
        assembly: vtkDataAssembly containing the hierarchical hierarchical organization
            of items in the vtkPartitionedDataSetCollection.
    """
    subset = None if subset_idx < 0 else grid_ggd.grid_subset[subset_idx]
    ugrid = read_geom.convert_grid_subset_geometry_to_unstructured_grid(
        grid_ggd, subset_idx, vtk_grid_points
    )
    ps_reader.read_plasma_state(subset_idx, ugrid)
    output.SetPartition(partition, 0, ugrid)
    label = str(subset.identifier.name) if subset else ids_name
    child = assembly.AddNode(label.replace(" ", "_"), 0)
    assembly.AddDataSetIndex(child, partition)
    output.GetMetaData(partition).Set(vtkCompositeDataSet.NAME(), label)
