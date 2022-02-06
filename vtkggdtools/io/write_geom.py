"""
A collection of methods to write a VTK dataset into a grid_ggd IDS node.
These methods populate the grid_ggd/space and grid_ggd/grid_subset children.
"""
from collections import defaultdict, OrderedDict
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

from identifiers.coordinate_identifier import coordinate_identifier
from identifiers.ggd_geometry_content_identifier import ggd_geometry_content_identifier
from identifiers.ggd_subset_identifier import ggd_subset_identifier
from vtkmodules.numpy_interface import algorithms as algs
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonCore import vtkMath, vtkIdList, vtkDataArray
from vtkmodules.vtkCommonDataModel import vtkCell, vtkCellIterator, vtkEdgeTable, vtkGenericCell, vtkPointSet, \
    vtkUnstructuredGrid, vtkPolygon, vtkCellData
from vtkmodules.vtkFiltersCore import vtkAppendDataSets


@dataclass
class GridGGDRepresentable:
    """
    A representation of the index transforms required to go from point and cell IDS
    in the VTK unstructured grid to the object indices in a grid_ggd IDS node.
    The indices are stored in C/python index convention (0, ... n)
    Prior to assigning these indices to a grid_ggd/**/attribute, increment by 1.
    """
    ugrid: vtkUnstructuredGrid = None
    edges: OrderedDict = field(default_factory=OrderedDict)  #: Go from edge (p1,p2) to 1d grid_ggd object index.
    faces: OrderedDict = field(default_factory=OrderedDict)  #: Go from face (p1,..) to 2d grid_ggd object index.
    volumes: OrderedDict = field(default_factory=OrderedDict)  #: Go from a 3D cell's id to 3d grid_ggd object index.


def fill_grid_ggd_basic_geometry(dataset: vtkPointSet, space_idx: int, grid_ggd) -> GridGGDRepresentable:
    """
    Copy the basic geometry objects (0d, 1d, 2d, 3d) and space coordinates into a grid_ggd.
    The points for the dataset go into the grid_ggd/space[space_idx] AoS element.
    :param dataset: can be any derived instance of a vtkPointSet. converted into a vtkUnstructuredGrid for convenience.
    :param space_idx: an index into the grid_ggd/space AoS. This space will be populated with the dataset's points.
    :param grid_ggd: the grid_ggd IDS node.
    :return: an instance of GridGGDRepresentable that helps go from vtk point ids to indices into grid_ggd/**/object.
    """
    # Convert input dataset into a vtkUnstructuredGrid
    append_filter = vtkAppendDataSets()
    append_filter.SetMergePoints(False)
    append_filter.AddInputData(dataset)
    append_filter.Update(0)
    output = append_filter.GetOutput()
    grid_ggd.space[space_idx].objects_per_dimension.resize(4)

    print('Writing nodes..')
    _fill_0d_objects(output, space_idx, grid_ggd)
    print('Writing edges..')
    edges = _fill_1d_objects(output, space_idx, grid_ggd)
    print('Writing faces..')
    faces = _fill_2d_objects(output, space_idx, edges, grid_ggd)
    print('Writing volumes..')
    volumes = _fill_3d_objects(output, space_idx, faces, grid_ggd)
    print('Writing implicit subsets..')
    _fill_implicit_grid_subsets(space_idx, grid_ggd)
    print('Finished')

    return GridGGDRepresentable(output, edges, faces, volumes)


def convert_vtk_dataset_to_grid_subset_geometry(representable: GridGGDRepresentable, space_idx: int,
                                                subset_idx: int, name: str,
                                                grid_ggd) -> None:
    """
    Populate the grid_ggd/grid_subset IDS node with the elements which have a >0 value for the cell data 'name'd array.
    :param representable: A GridGGDRepresentable.
    :param space_idx: an index into the grid_ggd/space AoS.
    :param subset_idx: an index into the grid_ggd/grid_subset AoS.
    :param name: the name of the cell data array. All cells with array value >0 will go into the grid_subset with name.
    :param grid_ggd: the grid_ggd IDS node.
    :return: None
    """
    print(f'Writing subset: {name} into grid_ggd/grid_subset[{subset_idx}]')

    is_custom_subset = ggd_subset_identifier.get(name) is None
    if not is_custom_subset:
        grid_ggd.grid_subset[subset_idx].identifier.name = name
        grid_ggd.grid_subset[subset_idx].identifier.index = ggd_subset_identifier.get(name).get('index')
        grid_ggd.grid_subset[subset_idx].identifier.description = ggd_subset_identifier.get(name).get('description')
    else:
        grid_ggd.grid_subset[subset_idx].identifier.name = name
        grid_ggd.grid_subset[subset_idx].identifier.index = -subset_idx
        grid_ggd.grid_subset[subset_idx].identifier.description = ''

    cell_data: vtkCellData = representable.ugrid.GetCellData()
    data_array: vtkDataArray = cell_data.GetArray(name)
    as_np_array = dsa.vtkDataArrayToVTKArray(data_array)
    cell_types = dsa.vtkDataArrayToVTKArray(representable.ugrid.GetCellTypesArray())
    if len(algs.where(as_np_array > 0)):
        subset_cell_ids = algs.where(as_np_array > 0)[0]
    else:
        return

    num_elements = len(subset_cell_ids)
    subset_cell_types = cell_types[subset_cell_ids]
    grid_ggd.grid_subset[subset_idx].element.resize(num_elements)
    gen_cell = vtkGenericCell()
    max_ndims = 0

    for i, cell_id, cell_type in zip(range(num_elements), subset_cell_ids, subset_cell_types):
        gen_cell.SetCellType(cell_type)
        dim = gen_cell.GetCellDimension()

        pt_ids = vtkIdList()
        representable.ugrid.GetCellPoints(cell_id, pt_ids)
        num_pts = pt_ids.GetNumberOfIds()
        pts = [pt_ids.GetId(j) for j in range(num_pts)]
        max_ndims = max(dim, max_ndims)

        if dim == 1:
            grid_ggd.grid_subset[subset_idx].element[i].object.resize(num_pts - 1)
            edge_counter = 0
            for p1, p2 in zip(pts[:-1], pts[1:]):
                edge_tup = (p1, p2)
                if edge_tup not in representable.edges and tuple(reversed(edge_tup)) in representable.edges:
                    edge_tup = (p2, p1)
                if edge_tup in representable.edges:
                    grid_ggd.grid_subset[subset_idx].element[i].object[edge_counter].index = representable.edges.get(
                        edge_tup) + 1
                    grid_ggd.grid_subset[subset_idx].element[i].object[edge_counter].space = space_idx + 1
                    grid_ggd.grid_subset[subset_idx].element[i].object[edge_counter].dimension = dim + 1
                edge_counter += 1

        elif dim == 2:
            grid_ggd.grid_subset[subset_idx].element[i].object.resize(1)
            face_tup = pts
            if face_tup not in representable.faces and tuple(reversed(face_tup)) in representable.faces:
                face_tup = tuple(reversed(face_tup))
            if face_tup in representable.faces:
                grid_ggd.grid_subset[subset_idx].element[i].object[0].index = representable.faces.get(face_tup) + 1
                grid_ggd.grid_subset[subset_idx].element[i].object[0].space = space_idx + 1
                grid_ggd.grid_subset[subset_idx].element[i].object[0].dimension = dim + 1

        elif dim == 3:
            grid_ggd.grid_subset[subset_idx].element[i].object.resize(1)
            if cell_id in representable.volumes:
                grid_ggd.grid_subset[subset_idx].element[i].object[0].index = representable.volumes.get(cell_id) + 1
                grid_ggd.grid_subset[subset_idx].element[i].object[0].space = space_idx + 1
                grid_ggd.grid_subset[subset_idx].element[i].object[0].dimension = dim + 1

    grid_ggd.grid_subset[subset_idx].dimension = max_ndims
    print('Finished')

def __add_unique_edge(cell: vtkCell, etbl: vtkEdgeTable, edges: OrderedDict, pts: List[int], object1d_idx: int) -> int:
    """
    A convenient function to add unique edge.
    (p1, p2) is treated the same as (p2, p1).
    :param cell: a vtkCell whose edges we shall add into edges map.
    :param etbl: a vtkEdgeTable to keep track of inserted edges.
    :param edges: an edge map (go from edge point ids to 1d object index in grid_ggd).
    :param pts: a working list of points.
    :param object1d_idx: the current 1d object index.
    :return: the no. of edges that were inserted.
    """
    num_pts = cell.GetNumberOfPoints()
    edge_added = 0
    for i in range(num_pts):
        pts[1] = cell.GetPointId(i)
        if (i > 0) and (etbl.IsEdge(pts[0], pts[1]) == -1):
            etbl.InsertEdge(pts[0], pts[1])
            edges.update({(pts[0], pts[1]): object1d_idx + edge_added})
            edge_added += 1
        pts[0] = pts[1]

    return edge_added


def __add_unique_face(cell: vtkCell, faces: OrderedDict, object2d_idx: int) -> bool:
    """
    A convenient function to add unique face.
    (p1, p2, p3, p4) is treated the same as (p4, p3, p2, p1).
    :param cell: a vtkCell representing a face which we shall add into faces map.
    :param faces: a face map (go from face point ids to 2d object index in grid_ggd).
    :param object2d_idx: the current 2d object index.
    :return: True if the face was added, else False.
    """
    num_pts = cell.GetNumberOfPoints()
    face = []

    for i in range(num_pts):
        pt = cell.GetPointId(i)
        face.append(pt)

    if tuple(reversed(face)) in faces or tuple(face) in faces:
        return False

    faces.update({tuple(face): object2d_idx})
    return True


def _fill_0d_objects(dataset: vtkUnstructuredGrid, space_idx: int, grid_ggd) -> None:
    """
    Fill all the 0D objects in the grid_ggd IDS node.
    :param dataset: a vtkUnstructuredGrid instance.
    :param space_idx: an index into the grid_ggd/space AoS.
    :param grid_ggd: the grid_ggd IDS node.
    :return: None
    """
    # Determine the size of space/objects_per_dim[0]/object/geometry node.
    geom_size = 3
    points = dataset.GetPoints()
    bds = points.GetBounds()
    if (bds[4] == bds[5]) and (bds[4] == 0):
        geom_size = 2

    # Populate geometry metadata for 0-d objects.
    space = grid_ggd.space[space_idx]
    space.coordinates_type.resize(geom_size)
    if geom_size == 2:
        space.coordinates_type[0] = coordinate_identifier.get('r').get('index')
        space.coordinates_type[1] = coordinate_identifier.get('z').get('index')
    else:
        space.coordinates_type[0] = coordinate_identifier.get('x').get('index')
        space.coordinates_type[1] = coordinate_identifier.get('y').get('index')
        space.coordinates_type[2] = coordinate_identifier.get('z').get('index')

    objects_per_dim = space.objects_per_dimension
    geometry_content = objects_per_dim[0].geometry_content
    geometry_content.name = 'node_coordinates'
    geometry_content.index = ggd_geometry_content_identifier.get('node_coordinates').get('index')
    geometry_content.description = ggd_geometry_content_identifier.get('node_coordinates').get('description')

    # Fill up 0-dim objects.
    objects_0d = objects_per_dim[0].object
    coords = dsa.vtkDataArrayToVTKArray(points.GetData())
    objects_0d.resize(len(coords))
    for i in range(len(objects_0d)):
        object_ = objects_0d[i]
        object_.nodes.resize(1)
        object_.nodes[0] = i + 1
        object_.geometry.resize(geom_size)
        for j in range(geom_size):
            object_.geometry[j] = coords[i][j]


def _fill_1d_objects(dataset: vtkUnstructuredGrid, space_idx: int, grid_ggd):
    """
    Fill all the 1D objects in the grid_ggd IDS node.
    :param dataset: a vtkUnstructuredGrid instance.
    :param space_idx: an index into the grid_ggd/space AoS.
    :param grid_ggd: the grid_ggd IDS node.
    :return: None
    """
    points = dataset.GetPoints()
    space = grid_ggd.space[space_idx]
    objects_per_dim = space.objects_per_dimension
    geometry_content = objects_per_dim[1].geometry_content
    geometry_content.name = 'edge_areas'
    geometry_content.index = ggd_geometry_content_identifier.get('edge_areas').get('index')
    geometry_content.description = ggd_geometry_content_identifier.get('edge_areas').get('description')

    coords = dsa.vtkDataArrayToVTKArray(points.GetData())
    cells_iter: vtkCellIterator = dataset.NewCellIterator()

    # Get number of cells per dimension.
    cells_iter.InitTraversal()
    cell = vtkGenericCell()
    edge_table = vtkEdgeTable()
    edges = OrderedDict()
    num_objects_per_dim = 0
    while not cells_iter.IsDoneWithTraversal():
        cells_iter.GetCell(cell)
        if cell.GetCellDimension() > 0:
            num_edges = cell.GetNumberOfEdges()
            pts = [0, 0]
            for e in range(num_edges):
                edge = cell.GetEdge(e)
                num_insertions = __add_unique_edge(edge, edge_table, edges, pts, num_objects_per_dim)
                num_objects_per_dim += num_insertions
            else:
                if num_edges == 0:
                    # cell_type == polyline or line
                    num_insertions = __add_unique_edge(cell, edge_table, edges, pts, num_objects_per_dim)
                    num_objects_per_dim += num_insertions
        cells_iter.GoToNextCell()

    # Build the neighbor edge map.
    nei_edge_ids = defaultdict(set)
    for edge, object_1d_idx in edges.items():
        nei_edge_ids[edge[0]].add(object_1d_idx)
        nei_edge_ids[edge[1]].add(object_1d_idx)

    # Populate 1D object
    objects_1d = objects_per_dim[1].object
    objects_1d.resize(num_objects_per_dim)
    for edge, object_1d_idx in edges.items():
        object_ = objects_1d[object_1d_idx]
        # object_.geometry.resize(1)
        object_.measure = vtkMath.Distance2BetweenPoints(coords[edge[0]], coords[edge[1]]) ** 0.5
        object_.nodes.resize(2)
        object_.boundary.resize(2)
        neighbours = []
        for i in range(2):
            object_.nodes[i] = edge[i] + 1
            object_.boundary[i].index = edge[i] + 1
            for nei_edge_id in nei_edge_ids[edge[i]]:
                if nei_edge_id != object_1d_idx:
                    neighbours.append(nei_edge_id)

        for i in range(2):
            object_.boundary[i].neighbours.resize(len(neighbours))
            for j, neighbour in enumerate(neighbours):
                object_.boundary[i].neighbours[j] = neighbour + 1

    return edges


def _fill_2d_objects(dataset: vtkUnstructuredGrid, space_idx: int, edges: Dict[Tuple[int, int], int],
                     grid_ggd) -> OrderedDict:
    """
    Fill all the 2D objects in the grid_ggd IDS node.
    :param dataset: a vtkUnstructuredGrid instance.
    :param space_idx: an index into the grid_ggd/space AoS.
    :param edges: an edge map to go from edge point ids to the 1d object index under grid_ggd.
    :param grid_ggd: the grid_ggd IDS node.
    :return: face map
    """
    points = dataset.GetPoints()
    space = grid_ggd.space[space_idx]
    objects_per_dim = space.objects_per_dimension
    geometry_content = objects_per_dim[2].geometry_content
    geometry_content.name = 'face_indices_volume'
    geometry_content.index = ggd_geometry_content_identifier.get('face_indices_volume').get('index')
    geometry_content.description = ggd_geometry_content_identifier.get('face_indices_volume').get('description')

    cells_iter: vtkCellIterator = dataset.NewCellIterator()

    # Get number of cells per dimension.
    cells_iter.InitTraversal()
    cell = vtkGenericCell()
    faces = OrderedDict()
    num_objects_per_dim = 0
    while not cells_iter.IsDoneWithTraversal():
        cells_iter.GetCell(cell)
        if cell.GetCellDimension() > 1:
            num_faces = cell.GetNumberOfFaces()
            for f in range(num_faces):
                face_cell = cell.GetFace(f)
                inserted = __add_unique_face(face_cell, faces, num_objects_per_dim)
                num_objects_per_dim += 1 if inserted else 0
            else:
                if num_faces == 0:
                    # cell_type == polygon or triangle or quad
                    inserted = __add_unique_face(cell, faces, num_objects_per_dim)
                    num_objects_per_dim += 1 if inserted else 0
        cells_iter.GoToNextCell()

    # Build the neighbor face map.
    nei_face_ids = defaultdict(set)
    for face, object_2d_idx in faces.items():
        for e in range(len(face)):
            p1 = face[e]
            p2 = face[(e + 1) % len(face)]
            nei_face_ids[(p1, p2)].add(object_2d_idx)
            nei_face_ids[(p2, p1)].add(object_2d_idx)

    # Populate 2D object
    objects_2d = objects_per_dim[2].object
    objects_2d.resize(num_objects_per_dim)
    for face, object_2d_idx in faces.items():
        object_ = objects_2d[object_2d_idx]
        num_face_points = len(face)
        num_face_edges = len(face)
        # object_.geometry.resize(1)
        poly = vtkPolygon()
        poly.Initialize(num_face_points, face, points)
        n = [0, 0, 1]
        poly.ComputeNormal(points, num_face_points, face, n)
        object_.measure = poly.ComputeArea(points, num_face_points, face, n)
        object_.nodes.resize(num_face_points)
        object_.boundary.resize(num_face_edges)
        for i in range(num_face_points):
            object_.nodes[i] = face[i] + 1

        neighbours = []
        for e in range(num_face_edges):
            p1 = face[e]
            p2 = face[(e + 1) % num_face_edges]
            edge_pts = [p1, p2]
            edge_tup = (p1, p2)

            if edge_tup not in edges and tuple(reversed(edge_pts)) in edges:
                edge_tup = tuple(reversed(edge_pts))
            if edge_tup in edges:
                object_.boundary[e].index = edges.get(edge_tup) + 1
            else:
                object_.boundary[e].index = -1

            for nei_face_id in nei_face_ids.get(edge_tup):
                if nei_face_id != object_2d_idx:
                    neighbours.append(nei_face_id)

        for e in range(num_face_edges):
            object_.boundary[e].neighbours.resize(len(neighbours))
            for j, neighbour in enumerate(neighbours):
                object_.boundary[e].neighbours[j] = neighbour + 1

    return faces


def _fill_3d_objects(dataset: vtkUnstructuredGrid, space_idx: int, faces: OrderedDict, grid_ggd) -> OrderedDict:
    """
    Fill all the 3D objects in the grid_ggd IDS node.
    :param dataset: a vtkUnstructuredGrid instance.
    :param space_idx: an index into the grid_ggd/space AoS.
    :param faces: an face map to go from face point ids to the 2d object index under grid_ggd.
    :param grid_ggd: the grid_ggd IDS node.
    :return: volume map
    """
    space = grid_ggd.space[space_idx]
    objects_per_dim = space.objects_per_dimension
    geometry_content = objects_per_dim[3].geometry_content
    geometry_content.name = 'unspecified'
    geometry_content.index = ggd_geometry_content_identifier.get('unspecified').get('index')
    geometry_content.description = ggd_geometry_content_identifier.get('unspecified').get('description')

    cells_iter: vtkCellIterator = dataset.NewCellIterator()

    cells_iter.InitTraversal()
    # Get number of cells per dimension.
    cell_ids = []
    cell_ndims = []
    cell_id_to_obj_id_map = OrderedDict()
    while not cells_iter.IsDoneWithTraversal():
        cell_dim = cells_iter.GetCellDimension()
        cell_id = cells_iter.GetCellId()
        cell_ndims.append(cell_dim)
        if cell_dim == 3:
            cell_id_to_obj_id_map.update({cell_id: len(cell_ids)})
            cell_ids.append(cell_id)
        cells_iter.GoToNextCell()

    # Populate 3D objects
    num_objects_per_dim = len(cell_ids)
    objects_3d = objects_per_dim[3].object
    objects_3d.resize(num_objects_per_dim)
    neighbours_3d_id_list = vtkIdList()

    for i, cell_id in enumerate(cell_ids):
        cell: vtkCell = dataset.GetCell(cell_id)
        object_ = objects_3d[i]
        num_faces = cell.GetNumberOfFaces()

        cell_pts = []
        object_.boundary.resize(num_faces)
        neighbours = []
        for f in range(num_faces):
            face_cell = cell.GetFace(f)
            num_face_pts = face_cell.GetNumberOfPoints()

            face_pts = []
            for j in range(num_face_pts):
                face_pts.append(face_cell.GetPointId(j))
                cell_pts.append(face_cell.GetPointId(j))

            face_tup = tuple(face_pts)
            if face_tup not in faces and tuple(reversed(face_pts)) in faces:
                face_tup = tuple(reversed(face_pts))
            if face_tup in faces:
                object_.boundary[f].index = faces.get(face_tup) + 1
            else:
                object_.boundary[f].index = -1

            dataset.GetCellNeighbors(cell_id, face_cell.GetPointIds(), neighbours_3d_id_list)

            for j in range(neighbours_3d_id_list.GetNumberOfIds()):
                nei_cell_id = neighbours_3d_id_list.GetId(j)
                if cell_ndims[nei_cell_id] == 3:
                    neighbours.append(nei_cell_id)

        for f in range(num_faces):
            object_.boundary[f].neighbours.resize(len(neighbours))
            for j, neighbour in enumerate(neighbours):
                object_.boundary[f].neighbours[j] = cell_id_to_obj_id_map[neighbour] + 1

        object_.nodes.resize(len(cell_pts))
        for j, pt in enumerate(cell_pts):
            object_.nodes[j] = pt + 1

    return cell_id_to_obj_id_map


def _fill_implicit_grid_subsets(space_idx: int, grid_ggd):
    """
    Fill up the implicit grid subsets for 0D, 1D, 2D, 3D objects.
    :param space_idx: an index into the grid_ggd/space AoS.
    :param grid_ggd: the grid_ggd IDS node.
    :return:
    """
    grid_subset_id_name = ['nodes', 'edges', 'cells', 'volumes']
    # Populate grid_subset and space.
    for i in range(4):
        subset_idx = i
        dim = i
        name = grid_subset_id_name[dim]
        grid_ggd.grid_subset[subset_idx].identifier.name = name
        grid_ggd.grid_subset[subset_idx].identifier.index = ggd_subset_identifier.get(name).get('index')
        grid_ggd.grid_subset[subset_idx].identifier.description = ggd_subset_identifier.get(name).get('description')
        grid_ggd.grid_subset[subset_idx].dimension = dim + 1
        num_elements = len(grid_ggd.space[space_idx].objects_per_dimension[dim].object)
        grid_ggd.grid_subset[subset_idx].element.resize(num_elements)

        # Populate grid_subset
        for j in range(num_elements):
            grid_ggd.grid_subset[subset_idx].element[j].object.resize(1)
            grid_ggd.grid_subset[subset_idx].element[j].object[0].space = space_idx + 1
            grid_ggd.grid_subset[subset_idx].element[j].object[0].dimension = dim + 1
            grid_ggd.grid_subset[subset_idx].element[j].object[0].index = j + 1
