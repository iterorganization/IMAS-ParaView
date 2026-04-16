"""
A collection of methods to read a grid_ggd IDS node into a VTK dataset.
These methods copy contents from the grid_ggd/space and grid_ggd/grid_subset
children into distinct vtkUnstructuredGrid objects.
"""

import logging

import numpy as np
from imas import identifiers
from imas.ids_structure import IDSStructure
from vtkmodules.vtkCommonCore import vtkIdList, vtkPoints
from vtkmodules.vtkCommonDataModel import (
    VTK_EMPTY_CELL,
    VTK_HEXAHEDRON,
    VTK_LINE,
    VTK_POLY_LINE,
    VTK_POLYGON,
    VTK_POLYHEDRON,
    VTK_PYRAMID,
    VTK_QUAD,
    VTK_TETRA,
    VTK_TRIANGLE,
    VTK_VERTEX,
    VTK_WEDGE,
    vtkUnstructuredGrid,
)

logger = logging.getLogger("imas_paraview")


def convert_grid_subset_geometry_to_unstructured_grid(
    grid_ggd, subset_idx: int, vtk_grid_points, progress=None
) -> vtkUnstructuredGrid:
    """Copy the elements found in given grid_ggd/grid_subset IDS node into a
    vtkUnstructuredGrid instance. This method uses the supplied point coordinates in
    the form of a vtkPoints instance.
    Args:
        grid_ggd: a grid_ggd ids node
        subset_idx: an index into grid_ggd/grid_subset
        vtk_grid_points: the point coordinates corresponding to 1d objects in
        the subset elements.
        progress: Progress indicator for Paraview.


    Returns:
        The vtkUnstructuredGrid containing the given points.
    """
    grid = vtkUnstructuredGrid()
    grid.SetPoints(vtk_grid_points)
    if subset_idx >= 0:
        _fill_vtk_cell_array_from_gs(grid_ggd, subset_idx, grid, progress)
    else:
        _fill_vtk_cell_array_from_gs2(grid_ggd, grid, progress)
    return grid


def fill_vtk_points(
    grid_ggd, space_idx: int, points: vtkPoints, ids_name: str, progress=None
) -> None:
    """Populate the vtkPoints data structure with coordinates from the grid_ggd/space
    IDS node.

    Args:
        grid_ggd: a grid_ggd ids node.
        space_idx: an index into the grid_ggd/space AoS.
        points: the vtk points instance.
        ids_name: The name of the IDS.
        progress: Progress indicator for Paraview.
    """
    num_objects0d = len(grid_ggd.space[space_idx].objects_per_dimension[0].object)
    logger.info(
        "Reading %d points from grid_ggd/space[%d]/objects_per_dimension[0]",
        num_objects0d,
        space_idx,
    )

    if len(grid_ggd.space[space_idx].objects_per_dimension[0].object[0].geometry) == 0:
        raise RuntimeError("Geometry of object is empty.")

    coordinates_type = grid_ggd.space[space_idx].coordinates_type
    # coordinates_type changed from INT_1D to an AoS of identifiers in DD4.0.0
    if isinstance(coordinates_type[0], IDSStructure):
        coord_indices = [int(ct.index) for ct in coordinates_type]
    else:
        coord_indices = [int(ct) for ct in coordinates_type]

    coord_id = identifiers.coordinate_identifier
    X, Y, Z = coord_id.x.value, coord_id.y.value, coord_id.z.value
    R, PHI = coord_id.r.value, coord_id.phi.value

    # Old version of GGD Fortran library (<=1.12.0) did not set coordinate identifiers
    # correctly, setting the points in r,phi-coordinates instead of r,z. For this case
    # we overwrite the coordinate identifiers.
    # This issue has been fixed in the following commit:
    # https://github.com/iterorganization/GGD/commit/23af2f113e550fa6e8d05c982ddae53bf29c1cf1 # noqa: E501
    if grid_ggd.space[space_idx].geometry_type.name == "Poloidal" and coord_indices == [
        R,
        PHI,
    ]:
        logger.warning(
            "The geometry type was set to 'Poloidal' but the coordinate identifiers "
            "were set to (r, phi). They have been interpreted as (r, z) instead."
        )
        coord_indices = [R, Z]

    supported = {X, Y, Z, R, PHI}
    unsupported = set(coord_indices) - supported
    if unsupported:
        logger.error(
            "Unsupported coordinate types in '%s', space[%d]: %s. They will be ignored",
            ids_name,
            space_idx,
            unsupported,
        )

    # Map coordinate identifier to position in geometry array
    coord_pos = {coord_index: idx for idx, coord_index in enumerate(coord_indices)}

    objects = grid_ggd.space[space_idx].objects_per_dimension[0].object
    points.Allocate(num_objects0d, 0)
    for obj in objects:
        geom = obj.geometry
        geom_len = len(geom)

        x = 0.0
        y = 0.0
        z = 0.0

        if X in coord_pos and coord_pos[X] < geom_len:
            x = geom[coord_pos[X]]
        if Y in coord_pos and coord_pos[Y] < geom_len:
            y = geom[coord_pos[Y]]
        if Z in coord_pos and coord_pos[Z] < geom_len:
            z = geom[coord_pos[Z]]

        # Handle cylindrical coordinates
        if R in coord_pos and coord_pos[R] < geom_len:
            r = geom[coord_pos[R]]
            phi = (
                geom[coord_pos[PHI]]
                if PHI in coord_pos and coord_pos[PHI] < geom_len
                else 0.0
            )
            x = r * np.cos(phi)
            y = r * np.sin(phi)

        points.InsertNextPoint(x, y, z)

        if progress:
            progress.increment(0.5 / num_objects0d)


def _fill_vtk_cell_array_from_gs2(
    grid_ggd, ugrid: vtkUnstructuredGrid, progress=None
) -> None:
    """_fill_vtk_cell_array_from_gs() for wall IDS.

    Args:
        grid_ggd: a grid_ggd ids node.
        ugrid: the vtk unstructured grid instance.
        progress: Progress indicator for Paraview.
    """
    grid = grid_ggd.space[0].objects_per_dimension
    target_indices = [2, 3]

    total_cells = 0
    for i in target_indices:
        if i < len(grid):
            total_cells += len(grid[i].object)
    if total_cells == 0:
        logger.info("No 2D or 3D objects found in grid_ggd/space[0].")
        return

    ugrid.AllocateEstimate(total_cells, 10)

    for obj_dimension in target_indices:
        if obj_dimension >= len(grid):
            continue

        objects = grid[obj_dimension].object

        if len(objects) == 0:
            continue

        logger.info(
            "Reading %d elements from space[0]/objects_per_dimension[%d]",
            len(objects),
            obj_dimension,
        )

        for obj in objects:
            if progress:
                progress.increment(0.5 / total_cells)

            obj_nodes = obj.nodes

            pt_ids = [val - 1 for val in obj_nodes]
            npts = len(pt_ids)
            cell_type = _get_vtk_cell_type(obj_dimension, npts)
            ugrid.InsertNextCell(cell_type, npts, pt_ids)


def _fill_vtk_cell_array_from_gs(
    grid_ggd, subset_idx: int, ugrid: vtkUnstructuredGrid, progress=None
) -> None:
    """Populate the cells in the vtk unstructured grid instance with elements from the
    grid_ggd/grid_subset IDS node.

    Args:
        grid_ggd: a grid_ggd ids node.
        subset_idx: an index into the grid_ggd/grid_subset AoS.
        ugrid: the vtk unstructured grid instance.
        progress: Progress indicator for Paraview.
    """
    grid_subset = grid_ggd.grid_subset[subset_idx]
    num_gs_el = len(grid_subset.element)

    if hasattr(grid_subset, "identifier"):
        logger.info(
            "Reading %d elements from %s", num_gs_el, grid_subset.identifier.name
        )
    else:
        logger.info(
            "Reading %d elements from grid_ggd/grid_subset[%d]", num_gs_el, subset_idx
        )

    ugrid.AllocateEstimate(num_gs_el, 10)
    object_3d_pt_ids = vtkIdList()

    for element in grid_subset.element:
        if progress:
            progress.increment(0.5 / num_gs_el)
        for object_ in element.object:
            obj_space = object_.space - 1
            obj_index = object_.index - 1
            obj_dimension = object_.dimension - 1
            obj_nodes = (
                grid_ggd.space[obj_space]
                .objects_per_dimension[obj_dimension]
                .object[obj_index]
                .nodes
            )
            obj_boundary = (
                grid_ggd.space[obj_space]
                .objects_per_dimension[obj_dimension]
                .object[obj_index]
                .boundary
            )

            # offset by -1 as fortran indexing used in IMAS( 1,...n)
            pt_ids = [val - 1 for val in obj_nodes]
            npts = len(pt_ids)
            cell_type = _get_vtk_cell_type(obj_dimension, npts)

            if cell_type != VTK_POLYHEDRON:
                ugrid.InsertNextCell(cell_type, npts, pt_ids)

            elif cell_type == VTK_POLYHEDRON:
                num_faces = len(obj_boundary)
                object_3d_pt_ids.Reset()
                object_3d_pt_ids.InsertNextId(num_faces)

                for f in range(num_faces):
                    object_2d_idx = obj_boundary[f].index - 1
                    object_2d = (
                        grid_ggd.space[obj_space]
                        .objects_per_dimension[2]
                        .object[object_2d_idx]
                    )
                    # offset by -1
                    object_2d_pt_ids = [val - 1 for val in object_2d.nodes]
                    num_face_pts = len(object_2d_pt_ids)

                    # the format for 3d cell point ids is
                    #   [numFace0Pts, p0, p1, .., numFace1Pts, p0, p1, ...]
                    object_3d_pt_ids.InsertNextId(num_face_pts)

                    for pt in object_2d_pt_ids:
                        object_3d_pt_ids.InsertNextId(pt)

                # the other overload (int, int, (int,...)) raises TypeError,
                # so use the one with vtkIdList
                ugrid.InsertNextCell(cell_type, object_3d_pt_ids)


def _get_vtk_cell_type(dimension: int, npts: int) -> int:
    """Determines a suitable VTK cell type from given cell dimensionality and number of
    points for that cell.

    Args:
        dimension: the number of dimensions for a cell.
        npts: the number of points for a cell.

    Returns:
        VTK Cell Type
    """
    if dimension == 0:
        return VTK_VERTEX

    elif dimension == 1:
        if npts == 2:
            return VTK_LINE
        else:
            return VTK_POLY_LINE

    elif dimension == 2:
        if npts == 3:
            return VTK_TRIANGLE
        elif npts == 4:
            return VTK_QUAD
        else:
            return VTK_POLYGON

    elif dimension == 3:
        vtk_map = {
            4: VTK_TETRA,
            5: VTK_PYRAMID,
            6: VTK_WEDGE,
            8: VTK_HEXAHEDRON,
        }
        return vtk_map.get(npts, VTK_POLYHEDRON)

    else:
        return VTK_EMPTY_CELL
