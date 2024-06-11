"""
A collection of methods to read a grid_ggd IDS node into a VTK dataset.
These methods copy contents from the grid_ggd/space and grid_ggd/grid_subset
children into distinct vtkUnstructuredGrid objects.
"""

from typing import Any, Callable

from paraview import logger as pvlog
from vtkmodules.vtkCommonCore import vtkIdList, vtkPoints
from vtkmodules.vtkCommonDataModel import (
    VTK_EMPTY_CELL,
    VTK_LINE,
    VTK_POLY_LINE,
    VTK_POLYGON,
    VTK_POLYHEDRON,
    VTK_QUAD,
    VTK_TRIANGLE,
    VTK_VERTEX,
    vtkUnstructuredGrid,
)


def convert_grid_subset_geometry_to_unstructured_grid(
    grid_ggd, subset_idx: int, vtk_grid_points
) -> vtkUnstructuredGrid:
    """
    Copy the elements found in given grid_ggd/grid_subset IDS node into a
    vtkUnstructuredGrid instance.
    This method uses the supplied point coordinates in the form of a vtkPoints instance.
    :param grid_ggd: a grid_ggd ids node
    :param subset_idx: an index into grid_ggd/grid_subset
    :param vtk_grid_points: the point coordinates corresponding to 1d objects in
        the subset elements.
    :return:
    """
    grid = vtkUnstructuredGrid()
    grid.SetPoints(vtk_grid_points)
    if subset_idx >= 0:
        _fill_vtk_cell_array_from_gs(grid_ggd, subset_idx, grid)
    else:
        _fill_vtk_cell_array_from_gs2(grid_ggd, subset_idx, grid)
    return grid


def fill_vtk_points(grid_ggd, space_idx: int, points: vtkPoints, ids_name: str) -> None:
    """
    Populate the vtkPoints data structure with coordinates from the grid_ggd/space IDS
    node.
    :param grid_ggd: a grid_ggd ids node.
    :param space_idx: an index into the grid_ggd/space AoS.
    :param points: the vtk points instance.
    :return: None
    """
    num_objects0d = len(grid_ggd.space[space_idx].objects_per_dimension[0].object)
    pvlog.info(f"Reading {num_objects0d} points from grid_ggd/space[{space_idx}]")

    s = 1  # scale objects from mm to m
    if grid_ggd.space[space_idx].objects_per_dimension[0].object[0].geometry[0] > 100:
        s = 0.001

    points.Allocate(num_objects0d, 0)
    coordinate_type = grid_ggd.space[space_idx].coordinates_type

    # When length of coordinate type is greater than 2, the third coordinate is in
    # geometry[2], else it is 0.
    if len(coordinate_type) > 2:
        third_dim: Callable[[Any], int] = lambda e: getattr(e, "geometry")[2]
    else:
        third_dim: Callable[[Any], int] = lambda e: 0

    for obj in grid_ggd.space[space_idx].objects_per_dimension[0].object:
        if ids_name == "wall":
            points.InsertNextPoint(
                (obj.geometry[0] * s, obj.geometry[1] * s, third_dim(obj) * s)
            )
        else:
            points.InsertNextPoint(
                (obj.geometry[0] * s, third_dim(obj) * s, obj.geometry[1] * s)
            )
            # Use for changing orientation in paraview.


def _fill_vtk_cell_array_from_gs2(
    grid_ggd, subset_idx: int, ugrid: vtkUnstructuredGrid
) -> None:
    """
    _fill_vtk_cell_array_from_gs() for wall IDS.
    :param grid_ggd: a grid_ggd ids node.
    :param subset_idx: -1
    :param ugrid: the vtk unstructured grid instance.
    :return: None
    """
    grid = grid_ggd.space[0].objects_per_dimension
    num_cell = len(grid[2].object[:])  # 4- 0D, 1D, 2D, 3D objects
    ugrid.AllocateEstimate(num_cell, 10)

    # Uses only 2d cells.
    for j in range(len(grid[2].object[:])):
        obj = grid[2].object[j]
        obj_nodes = obj.nodes
        obj_dimension = 2
        # try:
        #     obj_boundary = obj.boundary
        # except:
        #     obj_boundary = []

        pt_ids = list(map(lambda val: val - 1, obj_nodes))
        npts = len(pt_ids)
        cell_type = _get_vtk_cell_type(obj_dimension, npts)
        ugrid.InsertNextCell(cell_type, npts, pt_ids)


def _fill_vtk_cell_array_from_gs(
    grid_ggd, subset_idx: int, ugrid: vtkUnstructuredGrid
) -> None:
    """
    Populate the cells in the vtk unstructured grid instance with elements from the
    grid_ggd/grid_subset IDS node.
    :param grid_ggd: a grid_ggd ids node.
    :param subset_idx: an index into the grid_ggd/grid_subset AoS.
    :param ugrid: the vtk unstructured grid instance.
    :return: None
    """
    grid_subset = grid_ggd.grid_subset[subset_idx]
    num_gs_el = len(grid_subset.element)

    if hasattr(grid_subset, "identifier"):
        pvlog.info(f"Reading {num_gs_el} elements from {grid_subset.identifier.name}")
    else:
        pvlog.info(
            f"Reading {num_gs_el} elements from grid_ggd/grid_subset[{subset_idx}]"
        )

    ugrid.AllocateEstimate(num_gs_el, 10)
    object_3d_pt_ids = vtkIdList()

    for element in grid_subset.element:
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

            pt_ids = list(
                map(lambda val: val - 1, obj_nodes)
            )  # offset by -1 as fortran indexing used in IMAS( 1,...n)
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
                    object_2d_pt_ids = list(
                        map(lambda val: val - 1, object_2d.nodes)
                    )  # offset by -1
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
    """
    Determines a suitable VTK cell type from given cell dimensionality and number of
    points for that cell.
    :param dimension: the number of dimensions for a cell.
    :param npts: the number of points for a cell.
    :return: VTK Cell Type
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
        return VTK_POLYHEDRON

    else:
        return VTK_EMPTY_CELL
