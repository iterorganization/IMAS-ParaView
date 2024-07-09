import imaspy
import numpy as np

from ..util import create_first_ggd, create_first_grid, int64_to_int32


def fill_with_2_by_3_grid(grid_ggd):
    """
    Fills the grid_ggd of an IDS with a simple rectangular grid containing 6 vertices,
    7 edges and 2 faces, arranged in the following manner for a 2x3 grid:

             E5        E6
        P5---------P4---------P3
        |          |          |
    E4  |    F0    |    F1    |E3
        |          |          |
        P0---------P1---------P2
             E0        E1

    Adapted from https://sharepoint.iter.org/departments/POP/CM/IMDesign/Data%20Model/sphinx/3.41/ggd_guide/examples.html # noqa
    """
    num_vertices = 6
    num_edges = 7
    num_faces = 2

    # Set grid
    grid_ggd.identifier.name = "random_grid"
    grid_ggd.identifier.index = 1
    grid_ggd.identifier.description = "A simple 2x3 grid consisting of 6 vertices,\
    7 edges and 2 faces"

    # Set space
    grid_ggd.space.resize(1)
    space = grid_ggd.space[0]
    space.identifier.name = "random_space"
    space.identifier.index = 1
    space.identifier.description = "A simple 2x3 grid consisting of 6 vertices,\
    7 edges and 2 faces"
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

    # Create subset
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

    # Create elements for vertices
    grid_subsets[0].element.resize(num_vertices)
    for i, element in enumerate(grid_subsets[0].element):
        element.object.resize(1)
        element.object[0].space = 1
        element.object[0].dimension = 1
        element.object[0].index = i + 1

    # Create elements for edges
    grid_subsets[1].element.resize(num_edges)
    for i, element in enumerate(grid_subsets[1].element):
        element.object.resize(1)
        element.object[0].space = 1
        element.object[0].dimension = 2
        element.object[0].index = i + 1

    # Create elements for faces
    grid_subsets[2].element.resize(num_faces)
    for i, element in enumerate(grid_subsets[2].element):
        element.object.resize(1)
        element.object[0].space = 1
        element.object[0].dimension = 3
        element.object[0].index = i + 1

    return num_vertices, num_edges, num_faces


def fill_vector_quantity(quantity, num_vertices, num_edges, num_faces):
    """
    Fills vector quantity with with random data for each vertex, edge and face.
    Only sets poloidal and toroidal components.
    """
    # Allocate memory for 3 entries: vertices, edges and faces
    quantity.resize(3)

    # Fill values for vertices
    quantity[0].grid_index = 1
    quantity[0].grid_subset_index = 1
    quantity[0].poloidal = np.random.rand(num_vertices)
    quantity[0].toroidal = np.random.rand(num_vertices)

    # Fill values for edges
    quantity[1].grid_index = 1
    quantity[1].grid_subset_index = 2
    quantity[1].poloidal = np.random.rand(num_edges)
    quantity[1].toroidal = np.random.rand(num_edges)

    # Fill values for faces
    quantity[2].grid_index = 1
    quantity[2].grid_subset_index = 5
    quantity[2].poloidal = np.random.rand(num_faces)
    quantity[2].toroidal = np.random.rand(num_faces)


def fill_scalar_quantity(quantity, num_vertices, num_edges, num_faces):
    """
    Fills scalar quantity with random data for each vertex, edge and face.
    """
    # Allocate memory for 3 entries: vertices, edges and faces
    quantity.resize(3)

    # Set 6 vertices
    quantity[0].grid_index = 1
    quantity[0].grid_subset_index = 1
    quantity[0].values = np.random.rand(num_vertices)

    # Set 7 edges
    quantity[1].grid_index = 1
    quantity[1].grid_subset_index = 2
    quantity[1].values = np.random.rand(num_edges)

    # Set 2 faces
    quantity[2].grid_index = 1
    quantity[2].grid_subset_index = 5
    quantity[2].values = np.random.rand(num_faces)


def fill_structure(quantity, num_vertices, num_edges, num_faces):
    """
    Recursively fills an IDS structure with random values.
    """

    for subquantity in quantity:
        # If subquantity is a struct array
        if (
            subquantity.metadata.data_type
            == imaspy.ids_data_type.IDSDataType.STRUCT_ARRAY
        ):
            # Only fill struct array if it is empty
            if not subquantity.has_value:
                # Fill scalar quantity
                if subquantity.metadata.structure_reference == "generic_grid_scalar":
                    fill_scalar_quantity(
                        subquantity, num_vertices, num_edges, num_faces
                    )

                # Fill vector quantity
                elif (
                    subquantity.metadata.structure_reference
                    == "generic_grid_vector_components"
                ):
                    fill_vector_quantity(
                        subquantity, num_vertices, num_edges, num_faces
                    )

                # Recursively fill struct array
                else:
                    subquantity.resize(1)
                    fill_structure(subquantity[0], num_vertices, num_edges, num_faces)

        # If subquantity is a structure
        elif (
            subquantity.metadata.data_type == imaspy.ids_data_type.IDSDataType.STRUCTURE
        ):
            # Only fill structure if it is empty
            if not subquantity.has_value:

                # Recursively structure
                fill_structure(subquantity, num_vertices, num_edges, num_faces)


def fill_ggd(ggd, num_vertices, num_edges, num_faces):
    """
    Fills all generic grid scalar and generic grid vector components
    with random values.
    """

    # Fill IDS structure with random values
    fill_structure(ggd, num_vertices, num_edges, num_faces)


def main():

    factory = imaspy.IDSFactory()

    # Fetch list of all IDSs
    ids_list = factory.ids_names()

    filled_ids = []
    for ids_name in ids_list:
        print(f"----- Filling {ids_name} -----")
        grid_ggd_is_filled = False
        ggd_is_filled = False

        # Retrieve IDS method name for IDS creation
        method = getattr(factory, ids_name, None)
        if method is not None:
            ids = method()

            # Create an empty grid_ggd
            grid_ggd = create_first_grid(ids)

            # Skip filling grid_ggd if it does not exist
            if grid_ggd is None:
                print(f"{ids} has no grid_ggd")
            else:
                # Fill GGD grid with a simple 2x3 grid
                num_vertices, num_edges, num_faces = fill_with_2_by_3_grid(grid_ggd)
                grid_ggd_is_filled = True
                print(f"filled grid_ggd for {ids}")

            # Create an empty GGD
            ggd = create_first_ggd(ids)

            # Skip filling GGD if it does not exist
            if ggd is None:
                print(f"{ids} has no ggd")
            else:
                # Fill the GGD with random values
                fill_ggd(ggd, num_vertices, num_edges, num_faces)
                ggd_is_filled = True
                print(f"filled ggd for {ids}")

            if ggd_is_filled or grid_ggd_is_filled:
                # imaspy.util.print_tree(ids)
                filled_ids.append(ids)

    print(f"Succesfully filled the following IDSs:\n{filled_ids}")


if __name__ == "__main__":
    main()
