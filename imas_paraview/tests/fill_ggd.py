import logging

import imas
import numpy as np

from imas_paraview.ids_util import get_arrays_from_ids
from imas_paraview.util import create_first_ggd, create_first_grid, int32array

logger = logging.getLogger("imas_paraview")


def fill_NxN_grid(grid_ggd, N):
    """Fills the grid_ggd of an IDS with a uniform rectangular grid of size N x N,
    containing vertices, edges and faces.

    Adapted from https://sharepoint.iter.org/departments/POP/CM/IMDesign/Data%20Model/sphinx/3.41/ggd_guide/examples.html # noqa

    Args:
        grid_ggd: The GGD grid that will be filled with the N x N grid.
        N: The size of the N x N grid

    Returns:
        num_vertices: The number of vertices in the generated grid_ggd
        num_edges: The number of edges in the generated grid_ggd
        num_faces: The number of faces in the generated grid_ggd
    """

    # Set grid
    grid_ggd.identifier.name = "linear"
    grid_ggd.identifier.index = 1
    grid_ggd.identifier.description = f"A simple {N} x {N} rectangular grid"

    # Set space
    grid_ggd.space.resize(1)
    space = grid_ggd.space[0]
    space.identifier.name = "primary_standard"
    space.identifier.index = 1
    space.identifier.description = "Primary space defining the standard grid"
    space.geometry_type.index = 0
    space.coordinates_type.resize(2)
    space.coordinates_type[0] = 1
    space.coordinates_type[1] = 2

    space.objects_per_dimension.resize(3)
    num_vertices, num_edges, num_faces = build_grid(N, space)

    # Create subset
    grid_ggd.grid_subset.resize(3)
    grid_subsets = grid_ggd.grid_subset
    grid_subsets[0].dimension = 1
    grid_subsets[0].identifier.name = "vertices"
    grid_subsets[0].identifier.index = 1
    grid_subsets[0].identifier.description = "All vertices in the domain"

    grid_subsets[1].dimension = 2
    grid_subsets[1].identifier.name = "edges"
    grid_subsets[1].identifier.index = 2
    grid_subsets[1].identifier.description = "All edges in the domain"

    grid_subsets[2].dimension = 3
    grid_subsets[2].identifier.name = "faces"
    grid_subsets[2].identifier.index = 5
    grid_subsets[2].identifier.description = "All faces in the domain"

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


def build_grid(N, space):
    """Builds the vertices, edges and faces for an N x N grid.

    Args:
        N: Size of the grid
        space: Space AoS of the GGD grid

    Returns:
        num_vertices: The total number of vertices
        num_edges: The total number of edges
        num_faces: The total number of faces
    """
    num_vertices = set_vertices(N, space)
    num_edges = set_edges(N, space)
    num_faces = set_faces(N, space)
    return num_vertices, num_edges, num_faces


def set_vertices(N, space):
    """Sets the vertices for an N x N grid.

    Args:
        N: Size of the grid
        space: Space AoS of the GGD grid

    Returns:
        num_vertices: The total number of vertices
    """
    vertices = space.objects_per_dimension[0].object

    num_vertices = N * N
    vertices.resize(num_vertices)

    for i in range(N):
        for j in range(N):
            idx = i * N + j
            vertices[idx].geometry = [float(j), float(i)]
    return num_vertices


def set_edges(N, space):
    """Sets the edges for an N x N grid.

    Args:
        N: Size of the grid
        space: Space AoS of the GGD grid

    Returns:
        num_edges: The total number of edges
    """
    edges = space.objects_per_dimension[1].object

    num_edges = 2 * (N - 1) * N  # Total number of edges (horizontal + vertical)
    edges.resize(num_edges)

    edge_idx = 0

    # Set horizontal edges (left to right connections)
    for i in range(N):
        for j in range(N - 1):
            edges[edge_idx].nodes = int32array([i * N + j + 1, i * N + j + 2])
            edge_idx += 1

    # Set vertical edges (top to bottom connections)
    for i in range(N - 1):
        for j in range(N):
            edges[edge_idx].nodes = int32array([i * N + j + 1, (i + 1) * N + j + 1])
            edge_idx += 1
    return num_edges


def set_faces(N, space):
    """Sets the faces for an N x N grid.

    Args:
        N: Size of the grid
        space: Space AoS of the GGD grid

    Returns:
        num_faces: The total number of faces
    """
    face = space.objects_per_dimension[2].object

    num_faces = (N - 1) * (N - 1)  # Total number of faces (squares)
    face.resize(num_faces)

    face_idx = 0

    # Set faces for each square in the grid
    for i in range(N - 1):
        for j in range(N - 1):
            top_left = i * N + j + 1
            top_right = top_left + 1
            bottom_left = top_left + N
            bottom_right = bottom_left + 1
            face[face_idx].nodes = int32array(
                [top_left, top_right, bottom_right, bottom_left]
            )
            face_idx += 1
    return num_faces


def fill_vector_quantity(vector_quantity, num_vertices, num_edges, num_faces):
    """Fills vector quantity with with random data for each vertex, edge and face.
    Only the radial, poloidal and toroidal components of the vector quantity are filled.

    Args:
        vector_quantity: The vector quantity to be filled
        num_vertices: The number of vertices in the grid_ggd
        num_edges: The number of edges in the grid_ggd
        num_faces: The number of faces in the grid_ggd
    """
    # Allocate memory for 3 entries: vertices, edges and faces
    vector_quantity.resize(3)

    # Fill values for vertices
    vector_quantity[0].grid_index = 1
    vector_quantity[0].grid_subset_index = 1

    vector_quantity[0].radial = np.random.rand(num_vertices)
    vector_quantity[0].z = np.random.rand(num_vertices)

    # Fill values for edges
    vector_quantity[1].grid_index = 1
    vector_quantity[1].grid_subset_index = 2
    vector_quantity[1].radial = np.random.rand(num_edges)
    vector_quantity[1].z = np.random.rand(num_edges)

    # Fill values for faces
    vector_quantity[2].grid_index = 1
    vector_quantity[2].grid_subset_index = 5
    vector_quantity[2].radial = np.random.rand(num_faces)
    vector_quantity[2].z = np.random.rand(num_faces)


def fill_vector_rzphi_quantity(vector_quantity, num_vertices, num_edges, num_faces):
    """Fills vector rzphi quantity with with random data for each vertex, edge and face.
    Only the radial, toroidal and z components of the vector quantity are filled.

    Args:
        vector_quantity: The vector quantity to be filled
        num_vertices: The number of vertices in the grid_ggd
        num_edges: The number of edges in the grid_ggd
        num_faces: The number of faces in the grid_ggd
    """
    # Allocate memory for 3 entries: vertices, edges and faces
    vector_quantity.resize(3)

    # Fill values for vertices
    vector_quantity[0].grid_index = 1
    vector_quantity[0].grid_subset_index = 1
    vector_quantity[0].r = np.random.rand(num_vertices)
    vector_quantity[0].z = np.random.rand(num_vertices)

    # Fill values for edges
    vector_quantity[1].grid_index = 1
    vector_quantity[1].grid_subset_index = 2
    vector_quantity[1].r = np.random.rand(num_edges)
    vector_quantity[1].z = np.random.rand(num_edges)

    # Fill values for faces
    vector_quantity[2].grid_index = 1
    vector_quantity[2].grid_subset_index = 5
    vector_quantity[2].r = np.random.rand(num_faces)
    vector_quantity[2].z = np.random.rand(num_faces)


def fill_scalar_quantity(scalar_quantity, num_vertices, num_edges, num_faces):
    """Fills scalar quantity with random data for each vertex, edge and face.

    Args:
        scalar_quantity: The scalar quantity to be filled
        num_vertices: The number of vertices in the grid_ggd
        num_edges: The number of edges in the grid_ggd
        num_faces: The number of faces in the grid_ggd
    """
    # Allocate memory for 3 entries: vertices, edges and faces
    scalar_quantity.resize(3)

    # Set 6 vertices
    scalar_quantity[0].grid_index = 1
    scalar_quantity[0].grid_subset_index = 1
    scalar_quantity[0].values = np.random.rand(num_vertices)

    # Set 7 edges
    scalar_quantity[1].grid_index = 1
    scalar_quantity[1].grid_subset_index = 2
    scalar_quantity[1].values = np.random.rand(num_edges)

    # Set 2 faces
    scalar_quantity[2].grid_index = 1
    scalar_quantity[2].grid_subset_index = 5
    scalar_quantity[2].values = np.random.rand(num_faces)


def _generate_random_complex(num_entries):
    """Generates a list of random complex numbers.

    Args:
        num_entries: The length of the returned list

    Returns:
        a list of random complex numbers
    """
    return [np.random.rand() + 1j * np.random.rand() for _ in range(num_entries)]


def fill_complex_scalar_quantity(
    complex_scalar_quantity, num_vertices, num_edges, num_faces
):
    """Fills complex scalar quantity with random data for each vertex, edge and face.

    Args:
        complex_scalar_quantity: The scalar quantity to be filled
        num_vertices: The number of vertices in the grid_ggd
        num_edges: The number of edges in the grid_ggd
        num_faces: The number of faces in the grid_ggd
    """
    # Allocate memory for 3 entries: vertices, edges and faces
    complex_scalar_quantity.resize(3)

    # Set vertices
    complex_scalar_quantity[0].grid_index = 1
    complex_scalar_quantity[0].grid_subset_index = 1
    complex_scalar_quantity[0].values = _generate_random_complex(num_vertices)

    # Set edges
    complex_scalar_quantity[1].grid_index = 1
    complex_scalar_quantity[1].grid_subset_index = 2
    complex_scalar_quantity[1].values = _generate_random_complex(num_edges)

    # Set faces
    complex_scalar_quantity[2].grid_index = 1
    complex_scalar_quantity[2].grid_subset_index = 5
    complex_scalar_quantity[2].values = _generate_random_complex(num_faces)


def fill_ggd_data(ids, num_vertices, num_edges, num_faces):
    """Fills all generic grid scalar and generic grid vector components of a GGD
    with random values.

    Args:
        ids: The IDS for which the GGD will be filled
        num_vertices: The number of vertices in the grid_ggd
        num_edges: The number of edges in the grid_ggd
        num_faces: The number of faces in the grid_ggd
    """

    # Fill IDS structure with random values
    scalar_array_list, vector_array_list = get_arrays_from_ids(
        ids, get_empty_arrays=True, create_empty_structs=True
    )

    # Read scalar arrays
    for scalar_array in scalar_array_list:
        if scalar_array.metadata.structure_reference == "generic_grid_scalar":
            fill_scalar_quantity(scalar_array, num_vertices, num_edges, num_faces)
        elif scalar_array.metadata.structure_reference == "generic_grid_scalar_complex":
            fill_complex_scalar_quantity(
                scalar_array, num_vertices, num_edges, num_faces
            )

    # Read vector arrays
    for vector_array in vector_array_list:
        structure_reference = vector_array.metadata.structure_reference
        if structure_reference == "generic_grid_vector_components":
            fill_vector_quantity(vector_array, num_vertices, num_edges, num_faces)
        # From DDv4 onward `generic_grid_vector_components_rzphi` will be
        # replaced by `generic_grid_vector_components_rphiz`
        elif structure_reference in [
            "generic_grid_vector_components_rzphi",
            "generic_grid_vector_components_rphiz",
        ]:
            fill_vector_rzphi_quantity(vector_array, num_vertices, num_edges, num_faces)


def fill_ids(ids, time_steps=1, grid_size=2):
    """Fills the IDS with an N x N uniform GGD grid and fills all GGD arrays on this
    grid with random values.

    Args:
        ids: IDS to be filled.
        time_steps: Number of time steps to create in the IDS.
        grid_size: Size of the N x N grid. Defaults to 2, meaning a 2 x 2 grid.
    """

    # Create an empty grid_ggd
    grid_ggd = create_first_grid(ids)

    # Skip filling grid_ggd if it does not exist
    if grid_ggd is None:
        logger.warning(f"{ids.metadata.name} has no grid_ggd")
    else:
        # Create time steps
        ids.time = [float(t) for t in range(time_steps)]
        ids.ids_properties.homogeneous_time = imas.ids_defs.IDS_TIME_MODE_HOMOGENEOUS

        # Create grid and GGD AoS
        grid_ggd_aos = imas.util.get_parent(grid_ggd)

        grid_ggd_aos.resize(time_steps)
        ggd = create_first_ggd(ids)
        ggd_aos = imas.util.get_parent(ggd)
        ggd_aos.resize(time_steps)

        # Create uniform grids and fill them with random GGD data
        for i in range(time_steps):
            num_vertices, num_edges, num_faces = fill_NxN_grid(
                grid_ggd_aos[i], grid_size
            )
            logger.debug(f"filled grid_ggd at index {i}.")
        fill_ggd_data(ids, num_vertices, num_edges, num_faces)

    fill_ids_specific(ids)


def fill_ids_specific(ids):
    """Fill IDS-specific nodes such that these pass the IDS validation check.

    Args:
        ids: IDS to be filled.
    """
    ids_name = ids.metadata.name
    if ids_name == "wall":
        ids.description_ggd[0].thickness.resize(len(ids.time))
    elif ids_name == "runaway_electrons":
        ids.ggd_fluid.resize(len(ids.time))
