import imaspy
import numpy as np

from ..util import create_first_ggd, create_first_grid
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
        # If subquantity is a struct array
        if (
            subquantity.metadata.data_type
            == imaspy.ids_data_type.IDSDataType.STRUCT_ARRAY
        ):
            # Only fill struct array if it is empty
            if not subquantity.has_value:
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

        # If subquantity is a structure
        elif (
            subquantity.metadata.data_type == imaspy.ids_data_type.IDSDataType.STRUCTURE
        ):
            # Only fill structure if it is empty
            if not subquantity.has_value:
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
    # ids_list = factory.ids_names()
    ids_list = [
        "distribution_sources",
        "distributions",
        "edge_profiles",
        "edge_sources",
        "edge_transport",
        "equilibrium",
        "mhd",
        "radiation",
        "tf",
        "transport_solver_numerics",
        "wall",
        "waves",
    ]
    for ids_name in ids_list:
        # Retrieve IDS method name for IDS creation
        method = getattr(factory, ids_name, None)
        if method is not None:
            ids = method()

            # Create an empty GGD grid
            grid = create_first_grid(ids)

            # Skip this IDS if it does not contain a GGD grid or a GGD
            if grid is None:
                print(f"{ids} has no grid_ggd")
                continue

            # Fill GGD grid with a simple 2x3 grid
            num_vertices, num_edges, num_faces = fill_with_2_by_3_grid(grid)
            print(f"filled grid_ggd for {ids}")

            # Fill GGD with random values
            ggd = create_first_ggd(ids)

            # Skip this IDS if it does not contain a GGD
            if ggd is None:
                print(f"{ids} has no ggd")
                continue

            # Fill the GGD with random values
            fill_ggd(ggd, num_vertices, num_edges, num_faces)
            print(f"filled ggd for {ids}")

            imaspy.util.print_tree(ids)


if __name__ == "__main__":
    main()
