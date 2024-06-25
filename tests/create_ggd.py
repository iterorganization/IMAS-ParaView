# Adapted from https://sharepoint.iter.org/departments/POP/CM/IMDesign/Data%20Model/sphinx/3.41/ggd_guide/examples.html # noqa

import imaspy
import IPython


def fill_with_simple_grid(grid_ggd_input):
    """
    Fills the grid_ggd of IDS with a simple rectangular grid containing 4 vertices,
    4 edges and a single face, arranged in the following manner:

             E2
        P3---------P2
        |          |
    E3  |    C0    |E1
        |          |
        P0---------P1
             E0
    """

    # Set grid
    grid_ggd_input.resize(1)
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
    grid_subsets[1].identifier.name = "faces"
    grid_subsets[1].identifier.index = 2
    grid_subsets[1].identifier.description = "All lines in the domain"

    grid_subsets[2].dimension = 3
    grid_subsets[2].identifier.name = "cells"
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


def fill_physical_quantities(ggd):
    ggd.resize(1)

    # Fill time
    ggd[0].time = 1.0

    # Fill electron temperature
    ggd[0].electrons.temperature.resize(2)

    ggd[0].electrons.temperature[0].grid_index = 1
    ggd[0].electrons.temperature[0].grid_subset_index = 1
    ggd[0].electrons.temperature[0].values = [1.1, 1.2, 1.3, 1.4]

    ggd[0].electrons.temperature[1].grid_index = 1
    ggd[0].electrons.temperature[1].grid_subset_index = 2
    ggd[0].electrons.temperature[1].values = [2.5]

    # Fill ion species
    ggd[0].ion.resize(1)
    ggd[0].ion[0].label = "Ne+2"
    ggd[0].ion[0].state.resize(1)

    # Fill ion density
    ggd[0].ion[0].state[0].density.resize(2)

    ggd[0].ion[0].state[0].density[0].grid_index = 1
    ggd[0].ion[0].state[0].density[0].grid_subset_index = 1
    ggd[0].ion[0].state[0].density[0].values = [0.1, 0.2, 0.3, 0.4]

    ggd[0].ion[0].state[0].density[1].grid_index = 1
    ggd[0].ion[0].state[0].density[1].grid_subset_index = 2
    ggd[0].ion[0].state[0].density[1].values = [0.3]


def main():
    factory = imaspy.IDSFactory()

    # Create IDS containing ggd
    ep = factory.edge_profiles()
    ep.time.resize(1)
    ep.ids_properties.homogeneous_time = imaspy.ids_defs.IDS_TIME_MODE_HOMOGENEOUS

    # Create grid and fill with physical quantities
    grid_ggd = ep.grid_ggd
    ggd = ep.ggd

    fill_with_simple_grid(grid_ggd)
    fill_physical_quantities(ggd)
    imaspy.util.print_tree(ep)

    # Write ids to disk
    pulse, run, database = 1, 1, "example_ggd"
    entry = imaspy.DBEntry(imaspy.ids_defs.ASCII_BACKEND, database, pulse, run)
    entry.create()
    entry.put(ep)
    # IPython.embed()


if __name__ == "__main__":
    main()
