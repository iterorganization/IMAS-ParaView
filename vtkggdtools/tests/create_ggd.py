# Adapted from https://sharepoint.iter.org/departments/POP/CM/IMDesign/Data%20Model/sphinx/3.41/ggd_guide/examples.html # noqa

import imaspy
import numpy as np

from ..util import int64_to_int32


def fill_with_2_by_3_grid(grid_ggd):
    """
    Fills the grid_ggd of IDS with a simple rectangular grid containing 6 vertices,
    7 edges and 2 faces, arranged in the following manner for a 2x3 grid:

             E5        E6
        P5---------P4---------P3
        |          |          |
    E4  |    F0    |    F1    |E3
        |          |          |
        P0---------P1---------P2
             E0        E1

    Expects the grid_ggd to have a nonzero size and to be passed the first grid
    """

    num_vertices = 6
    num_edges = 7
    num_faces = 2

    # Set grid
    grid_ggd.identifier.name = "grid_example_2"
    grid_ggd.identifier.index = 1
    grid_ggd.identifier.description = "Grid - example 2"

    # Set space
    grid_ggd.space.resize(1)
    space = grid_ggd.space[0]
    space.identifier.name = "space_example_2"
    space.identifier.index = 1
    space.identifier.description = "Space - example 2"
    space.geometry_type.index = 0
    space.coordinates_type = int64_to_int32([1, 2])

    space.objects_per_dimension.resize(3)
    vertices = space.objects_per_dimension[0].object
    edges = space.objects_per_dimension[1].object
    face = space.objects_per_dimension[2].object

    # Set Vertices
    vertices.resize(num_vertices)
    vertices[0].geometry = [0.0, 0.0]  # P0
    vertices[1].geometry = [1.0, 0.0]  # P1
    vertices[2].geometry = [2.0, 0.0]  # P2
    vertices[3].geometry = [2.0, 1.0]  # P3
    vertices[4].geometry = [1.0, 1.0]  # P4
    vertices[5].geometry = [0.0, 1.0]  # P5

    # Set edges
    edges.resize(num_edges)
    edges[0].nodes = int64_to_int32([1, 2])  # E0 (P0 -> P1)
    edges[1].nodes = int64_to_int32([2, 3])  # E1 (P1 -> P2)
    edges[2].nodes = int64_to_int32([3, 4])  # E2 (P2 -> P3)
    edges[3].nodes = int64_to_int32([4, 5])  # E3 (P3 -> P4)
    edges[4].nodes = int64_to_int32([5, 6])  # E4 (P4 -> P5)
    edges[5].nodes = int64_to_int32([6, 1])  # E5 (P5 -> P0)
    edges[6].nodes = int64_to_int32([2, 5])  # E6 (P1 -> P4)

    # Set faces
    face.resize(num_faces)
    face[0].nodes = int64_to_int32([1, 2, 5, 6])  # F0 (P0, P1, P4, P5)
    face[1].nodes = int64_to_int32([2, 3, 4, 5])  # F1 (P1, P2, P3, P4)

    # Set subset
    grid_ggd.grid_subset.resize(3)
    grid_subsets = grid_ggd.grid_subset
    grid_subsets[0].dimension = 1
    grid_subsets[0].identifier.name = "nodes"
    grid_subsets[0].identifier.index = 1
    grid_subsets[0].identifier.description = "All nodes in the domain"

    grid_subsets[1].dimension = 2
    grid_subsets[1].identifier.name = "edges"
    grid_subsets[1].identifier.index = 2
    grid_subsets[1].identifier.description = "All edges in the domain"

    grid_subsets[2].dimension = 3
    grid_subsets[2].identifier.name = "faces"
    grid_subsets[2].identifier.index = 5
    grid_subsets[2].identifier.description = "All 2D cells in the domain"

    # Set elements for vertices
    grid_subsets[0].element.resize(6)
    for i, element in enumerate(grid_subsets[0].element):
        element.object.resize(1)
        element.object[0].space = 1
        element.object[0].dimension = 1
        element.object[0].index = i + 1

    # Set elements for edges
    grid_subsets[1].element.resize(7)
    for i, element in enumerate(grid_subsets[1].element):
        element.object.resize(1)
        element.object[0].space = 1
        element.object[0].dimension = 2
        element.object[0].index = i + 1

    # Set elements for faces
    grid_subsets[2].element.resize(2)
    for i, element in enumerate(grid_subsets[2].element):
        element.object.resize(1)
        element.object[0].space = 1
        element.object[0].dimension = 3
        element.object[0].index = i + 1

    return num_vertices, num_edges, num_faces


def fill_physical_quantities_2_by_3(ggd):

    # Fill time
    ggd[0].time = 1.0

    # Fill electron temperature
    ggd[0].electrons.temperature.resize(3)

    # Values for the 6 vertices
    ggd[0].electrons.temperature[0].grid_index = 1
    ggd[0].electrons.temperature[0].grid_subset_index = 1  # Nodes subset
    ggd[0].electrons.temperature[0].values = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6]

    ggd[0].electrons.temperature[1].grid_index = 1
    ggd[0].electrons.temperature[1].grid_subset_index = 2
    ggd[0].electrons.temperature[1].values = [2, 3, 4, 5, 6, 7, 8]

    ggd[0].electrons.temperature[2].grid_index = 1
    ggd[0].electrons.temperature[2].grid_subset_index = 5
    ggd[0].electrons.temperature[2].values = [3, 4]


def fill_with_simple_grid(grid_ggd_input):
    """
    Fills the grid_ggd of IDS with a simple rectangular grid containing 4 vertices,
    4 edges and a single face, arranged in the following manner:

             E2
        P3---------P2
        |          |
    E3  |    F0    |E1
        |          |
        P0---------P1
             E0
    """

    num_vertices = 4
    num_edges = 4
    num_faces = 1

    # Set grid
    grid_ggd = grid_ggd_input[0]
    grid_ggd.identifier.name = "grid_example_1"
    grid_ggd.identifier.index = 1
    grid_ggd.identifier.description = "Grid - example 1"

    # Set space
    grid_ggd.space.resize(1)
    space = grid_ggd.space[0]
    space.identifier.name = "space_example_1"
    space.identifier.index = 1
    space.identifier.description = "Space - example 1"
    space.geometry_type.index = 0
    space.coordinates_type = [1, 2]

    space.objects_per_dimension.resize(3)
    vertices = space.objects_per_dimension[0].object
    edges = space.objects_per_dimension[1].object
    face = space.objects_per_dimension[2].object

    # Set Vertices
    vertices.resize(4)
    vertices[0].geometry = [0.0, 0.0]
    vertices[1].geometry = [1.0, 0.0]
    vertices[2].geometry = [1.0, 1.0]
    vertices[3].geometry = [0.0, 1.0]

    # Set edges
    edges.resize(4)
    edges[0].nodes = [1, 2]
    edges[1].nodes = [2, 3]
    edges[2].nodes = [3, 4]
    edges[3].nodes = [4, 1]

    # Set face
    face.resize(1)
    face[0].nodes = [1, 2, 3, 4]

    # Set subset
    grid_ggd.grid_subset.resize(3)
    grid_subsets = grid_ggd.grid_subset
    grid_subsets[0].dimension = 1
    grid_subsets[0].identifier.name = "nodes"
    grid_subsets[0].identifier.index = 1
    grid_subsets[0].identifier.description = "All nodes in the domain"

    grid_subsets[1].dimension = 2
    grid_subsets[1].identifier.name = "edges"
    grid_subsets[1].identifier.index = 2
    grid_subsets[1].identifier.description = "All lines in the domain"

    grid_subsets[2].dimension = 3
    grid_subsets[2].identifier.name = "faces"
    grid_subsets[2].identifier.index = 5
    grid_subsets[2].identifier.description = "All 2D cells in the domain"

    # Set elements for vertices
    grid_subsets[0].element.resize(4)
    for i, element in enumerate(grid_subsets[0].element):

        element.object.resize(1)
        element.object[0].space = 1
        element.object[0].dimension = 1
        element.object[0].index = i + 1

    # Set elements for edges
    grid_subsets[1].element.resize(4)
    for i, element in enumerate(grid_subsets[1].element):

        element.object.resize(1)
        element.object[0].space = 1
        element.object[0].dimension = 2
        element.object[0].index = i + 1

    # Set elements for face
    grid_subsets[2].element.resize(1)
    grid_subsets[2].element[0].object.resize(1)
    grid_subsets[2].element[0].object[0].space = 1
    grid_subsets[2].element[0].object[0].dimension = 3
    grid_subsets[2].element[0].object[0].index = 1

    return num_vertices, num_edges, num_faces


def fill_physical_quantities(ggd):

    # Fill time
    ggd[0].time = 1.0

    # Fill electron temperature
    ggd[0].electrons.temperature.resize(3)

    ggd[0].electrons.temperature[0].grid_index = 1
    ggd[0].electrons.temperature[0].grid_subset_index = 1
    ggd[0].electrons.temperature[0].values = [1.1, 1.2, 1.3, 1.4]

    ggd[0].electrons.temperature[0].grid_index = 1
    ggd[0].electrons.temperature[0].grid_subset_index = 2
    ggd[0].electrons.temperature[0].values = [2.1, 2.2, 2.3, 2.4]

    ggd[0].electrons.temperature[1].grid_index = 1
    ggd[0].electrons.temperature[1].grid_subset_index = 5
    ggd[0].electrons.temperature[1].values = [2.5]


def main():
    factory = imaspy.IDSFactory()

    # Create IDS containing ggd
    ep = factory.edge_profiles()
    ep.time.resize(1)
    ep.ids_properties.homogeneous_time = imaspy.ids_defs.IDS_TIME_MODE_HOMOGENEOUS

    # Create grid and fill with physical quantities
    grid_ggd = ep.grid_ggd
    ggd = ep.ggd

    fill_with_2_by_3_grid(grid_ggd)
    fill_physical_quantities_2_by_3(ggd)
    # fill_with_simple_grid(grid_ggd)
    # fill_physical_quantities(ggd)
    imaspy.util.print_tree(ep)

    # Write ids to disk
    pulse, run, database = 1, 1, "example_ggd"
    entry = imaspy.DBEntry(imaspy.ids_defs.ASCII_BACKEND, database, pulse, run)
    entry.create()
    entry.put(ep)


if __name__ == "__main__":
    main()
