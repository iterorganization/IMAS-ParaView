import random

import imaspy
import numpy as np
from create_ggd import fill_with_2_by_3_grid


def set_vector_quantity(quantity):
    """
    Fills vector quantity with random data, only sets radial and toroidal value
    """
    # Allocate memory for 3 entries: nodes, edges and faces
    quantity.resize(3)

    # Set 6 vertices
    quantity[0].grid_index = 1
    quantity[0].grid_subset_index = 1
    quantity[0].radial = np.random.rand(6)
    quantity[0].toroidal = np.random.rand(6)

    # Set 7 edges
    quantity[1].grid_index = 1
    quantity[1].grid_subset_index = 2
    quantity[1].radial = np.random.rand(7)
    quantity[1].toroidal = np.random.rand(7)

    # Set 2 faces
    quantity[2].grid_index = 1
    quantity[2].grid_subset_index = 5
    quantity[2].radial = np.random.rand(2)
    quantity[2].toroidal = np.random.rand(2)


def set_scalar_quantity(quantity):
    """
    Fills scalar quantity with random data
    """
    # Allocate memory for 3 entries: nodes, edges and faces
    quantity.resize(3)

    # Set 6 vertices
    quantity[0].grid_index = 1
    quantity[0].grid_subset_index = 1
    quantity[0].values = np.random.rand(6)

    # Set 7 edges
    quantity[1].grid_index = 1
    quantity[1].grid_subset_index = 2
    quantity[1].values = np.random.rand(7)

    # Set 2 faces
    quantity[2].grid_index = 1
    quantity[2].grid_subset_index = 5
    quantity[2].values = np.random.rand(2)


def fill_structure(quantity):
    """
    Recursively fills an IDS structure
    """
    if quantity.metadata.data_type == imaspy.ids_data_type.IDSDataType.STRUCTURE:
        print(f"Filling structure: {quantity}")

    # Handle IDS struct arrays
    elif quantity.metadata.data_type == imaspy.ids_data_type.IDSDataType.STRUCT_ARRAY:
        print(f"Filling struct_array: {quantity}")
        # Fill scalar or vector quantity
        if quantity.metadata.structure_reference == "generic_grid_scalar":
            set_scalar_quantity(quantity)
            return
        elif quantity.metadata.structure_reference == "generic_grid_vector_components":
            set_vector_quantity(quantity)
            return
        quantity.resize(1)
        quantity = quantity[0]

    # skip time
    elif quantity.metadata.name == "time":
        print(f"Skipped {quantity} quantity")
        return
    else:
        print("ERROR NOT RECOGNISED")
        print(quantity)

    for subquantity in quantity:
        print(f"{subquantity=}")

        # FLT
        if (
            subquantity.metadata.data_type == imaspy.ids_data_type.IDSDataType.FLT
            or subquantity.metadata.data_type == imaspy.ids_data_type.IDSDataType.INT
        ):
            print(f"{subquantity} IS A INT/FLT")
            # FIXME: make random float
            # imaspy.util.inspect(subquantity)
            # subquantity = np.random.rand(1)
            # imaspy.util.inspect(subquantity)

        elif subquantity.metadata.data_type == imaspy.ids_data_type.IDSDataType.STR:
            print(f"{subquantity} IS A STR")
            # FIXME: make random str?
            # imaspy.util.inspect(subquantity)
            # subquantity = np.random.rand(1)
            # imaspy.util.inspect(subquantity)

        # STRUCT ARRAY
        elif (
            subquantity.metadata.data_type
            == imaspy.ids_data_type.IDSDataType.STRUCT_ARRAY
        ):
            print(f"{subquantity} IS A STRUCT_ARRAY")

            # Fill scalar or vector quantity
            if subquantity.metadata.structure_reference == "generic_grid_scalar":
                set_scalar_quantity(subquantity)
            elif (
                subquantity.metadata.structure_reference
                == "generic_grid_vector_components"
            ):
                set_vector_quantity(subquantity)
            # Recursively fill structure
            else:
                fill_structure(subquantity)

        # STRUCTURE
        elif (
            subquantity.metadata.data_type == imaspy.ids_data_type.IDSDataType.STRUCTURE
        ):
            print(f"{subquantity} IS A STRUCTURE")
            fill_structure(subquantity)

        else:
            print(
                f"Error, did not recognise {subquantity.metadata.name} with metadata data type "
                f"{subquantity.metadata.data_type}"
            )
            exit()


def fill_ids(ggd):
    # Create GGD for single time step
    ggd.resize(1)
    ggd = ggd[0]
    ggd.time = 1.0

    print("-------------------------")
    for quantity in ggd:
        fill_structure(quantity)
        print(f"FINISHED FILLING {quantity}")
        if quantity.metadata.name != "time":
            # Show fully filled quantity
            imaspy.util.print_tree(quantity)

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
            fill_with_2_by_3_grid(ids.grid_ggd)
            fill_ids(ids.ggd)


if __name__ == "__main__":
    main()
