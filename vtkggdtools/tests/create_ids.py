import imaspy
import numpy as np

from ..util import get_ggd_grid_path
from .create_ggd import fill_with_2_by_3_grid


def set_vector_quantity(quantity, num_vertices, num_edges, num_faces):
    """
    Fills vector quantity with random data, only sets radial and toroidal value
    """
    # Allocate memory for 3 entries: nodes, edges and faces
    quantity.resize(3)

    # Set 6 vertices
    quantity[0].grid_index = 1
    quantity[0].grid_subset_index = 1
    quantity[0].radial = np.random.rand(num_vertices)
    quantity[0].toroidal = np.random.rand(num_vertices)

    # Set 7 edges
    quantity[1].grid_index = 1
    quantity[1].grid_subset_index = 2
    quantity[1].radial = np.random.rand(num_edges)
    quantity[1].toroidal = np.random.rand(num_edges)

    # Set 2 faces
    quantity[2].grid_index = 1
    quantity[2].grid_subset_index = 5
    quantity[2].radial = np.random.rand(num_faces)
    quantity[2].toroidal = np.random.rand(num_faces)


def set_scalar_quantity(quantity, num_vertices, num_edges, num_faces):
    """
    Fills scalar quantity with random data
    """
    # Allocate memory for 3 entries: nodes, edges and faces
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
    Recursively fills an IDS structure
    """

    for subquantity in quantity:
        print(f"{subquantity=}")

        # STRUCT ARRAY
        if (
            subquantity.metadata.data_type
            == imaspy.ids_data_type.IDSDataType.STRUCT_ARRAY
        ):
            print(f"Filling struct_array: {subquantity}")
            # Fill scalar quantity
            if subquantity.metadata.structure_reference == "generic_grid_scalar":
                set_scalar_quantity(subquantity, num_vertices, num_edges, num_faces)

            # Fill vector quantity
            elif (
                subquantity.metadata.structure_reference
                == "generic_grid_vector_components"
            ):
                set_vector_quantity(subquantity, num_vertices, num_edges, num_faces)
            # Recursively fill structure
            else:
                subquantity.resize(1)
                fill_structure(subquantity[0], num_vertices, num_edges, num_faces)

        # STRUCTURE
        elif (
            subquantity.metadata.data_type == imaspy.ids_data_type.IDSDataType.STRUCTURE
        ):
            print(f"Filling structure: {subquantity}")
            fill_structure(subquantity, num_vertices, num_edges, num_faces)


def fill_ids(ggd, num_vertices, num_edges, num_faces):
    """
    Resizes the GGD and fills all generic grid scalar and generic grid vector components
    with random values.
    """
    # Create GGD for single time step
    ggd.resize(1)
    ggd = ggd[0]
    ggd.time = 1.0

    # Fill IDS structure with random values
    fill_structure(ggd, num_vertices, num_edges, num_faces)

    # Show fully filled ggd
    imaspy.util.print_tree(ggd)


def main():

    # ids_list = ["edge_profiles", "edge_sources"]
    ids_list = ["edge_profiles"]
    factory = imaspy.IDSFactory()
    for ids_name in ids_list:
        # Retrieve IDS method name for IDS creation
        method = getattr(factory, ids_name, None)
        if method is not None:
            ids = method()
            imaspy.util.inspect(ids)
            num_vertices, num_edges, num_faces = fill_with_2_by_3_grid(ids.grid_ggd)
            fill_ids(ids.ggd, num_vertices, num_edges, num_faces)


if __name__ == "__main__":
    main()
