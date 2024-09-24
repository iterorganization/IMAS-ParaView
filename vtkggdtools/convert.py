# ggd_to_vtk.py

import logging

from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import (
    vtkCompositeDataSet,
    vtkDataAssembly,
    vtkPartitionedDataSetCollection,
)

from vtkggdtools.io import read_geom, read_ps
from vtkggdtools.util import get_grid_ggd

logger = logging.getLogger("vtkggdtools")


def ggd_to_vtk(
    ids,
    time_step_idx,
    scalar_paths_to_load=None,
    vector_paths_to_load=None,
    outInfo=None,
):
    ps_reader = read_ps.PlasmaStateReader(ids)

    # TODO: allow selecting other grids for Bezier

    # TODO: make time actual time not index
    grid_ggd = get_grid_ggd(ids, time_step_idx)

    # Check if grid is valid
    if grid_ggd is None:
        logger.warning("Could not load a valid GGD grid.")
        return None
    if not hasattr(grid_ggd, "space") or len(grid_ggd.space) < 1:
        logger.warning("The grid_ggd does not contain a space.")
        return None

    if outInfo is None:
        output = vtkPartitionedDataSetCollection()
    else:
        output = vtkPartitionedDataSetCollection.GetData(outInfo)
    num_subsets = len(grid_ggd.grid_subset)
    points = vtkPoints()
    space_idx = 0
    idsname = ids.metadata.name

    read_geom.fill_vtk_points(grid_ggd, space_idx, points, idsname)
    assembly = vtkDataAssembly()
    output.SetDataAssembly(assembly)

    # TODO add jorek functionality
    # TODO fix updateprogress

    # _aos_index_values = FauxIndexMap()
    # # Interpolate JOREK Fourier space
    # # TODO: figure out if we can put this functionality in a post-processing step?
    # if self._n_plane != 0:
    #     number_of_spaces = len(grid_ggd.space)
    #     if number_of_spaces > 1 and len(grid_ggd.space[0].coordinates_type) == 2:
    #         n_period = grid_ggd.space[1].geometry_type.index
    #         if n_period > 0:  # Fourier space with periodicity (JOREK)
    #             logger.info(f"Reading Bezier mesh with Fourier periodiciy {n_period}")
    #             ugrid = read_bezier.convert_grid_subset_to_unstructured_grid(
    #                 idsname,
    #                 ids,
    #                 _aos_index_values,
    #                 self._n_plane,
    #                 self._phi_start,
    #                 self._phi_end,
    #             )
    #             output.SetPartition(0, 0, ugrid)
    #             child = assembly.AddNode(idsname, 0)
    #             assembly.AddDataSetIndex(child, 0)
    #             output.GetMetaData(0).Set(vtkCompositeDataSet.NAME(), idsname)
    #         else:
    #             logger.error(
    #                 f"Number of planes {self._n_plane} invalid for this "
    #                 f"{number_of_spaces} number of spaces"
    #             )
    #     else:
    #         logger.error(
    #             f"Number of planes {self._n_plane} invalid for this IDS type."
    #             " Try using N = 0"
    #         )
    #     return 1

    # Load the GGD arrays from the selected GGD paths
    ps_reader.load_arrays_from_path(
        time_step_idx, scalar_paths_to_load, vector_paths_to_load
    )

    if num_subsets <= 1:
        logger.info("No subsets to read from grid_ggd")
        output.SetNumberOfPartitionedDataSets(1)
        _fill_grid_and_plasma_state(
            -1, 0, points, output, assembly, grid_ggd, ids, ps_reader
        )
    elif idsname == "wall":
        # FIXME: what if num_subsets is 2 or 3?
        output.SetNumberOfPartitionedDataSets(num_subsets - 3)
        _fill_grid_and_plasma_state(
            -1, 0, points, output, assembly, grid_ggd, ids, ps_reader
        )
        for subset_idx in range(4, num_subsets):
            _fill_grid_and_plasma_state(
                subset_idx,
                subset_idx - 3,
                points,
                output,
                assembly,
                grid_ggd,
                ids,
                ps_reader,
            )
            # self.UpdateProgress(self.GetProgress() + 1 / num_subsets)
    else:
        output.SetNumberOfPartitionedDataSets(num_subsets)
        for subset_idx in range(num_subsets):
            _fill_grid_and_plasma_state(
                subset_idx,
                subset_idx,
                points,
                output,
                assembly,
                grid_ggd,
                ids,
                ps_reader,
            )
            # self.UpdateProgress(self.GetProgress() + 1 / num_subsets)

    logger.info("Finished loading IDS.")
    return output


def _fill_grid_and_plasma_state(
    subset_idx, partition, points, output, assembly, grid_ggd, ids, ps_reader
):
    subset = None if subset_idx < 0 else grid_ggd.grid_subset[subset_idx]
    ugrid = read_geom.convert_grid_subset_geometry_to_unstructured_grid(
        grid_ggd, subset_idx, points
    )
    ps_reader.read_plasma_state(subset_idx, ugrid)
    output.SetPartition(partition, 0, ugrid)
    label = str(subset.identifier.name) if subset else ids.metadata.name
    child = assembly.AddNode(label.replace(" ", "_"), 0)
    assembly.AddDataSetIndex(child, partition)
    output.GetMetaData(partition).Set(vtkCompositeDataSet.NAME(), label)
