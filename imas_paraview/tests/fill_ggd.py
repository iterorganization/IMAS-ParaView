import logging

import imas
import numpy as np
from imas import identifiers
from imas.ids_structure import IDSStructure

from imas_paraview.ids_util import get_arrays_from_ids
from imas_paraview.util import create_first_ggd, create_first_grid, int32array

logger = logging.getLogger("imas_paraview")


def fill_NxN_grid(grid_ggd, N, create_3d_grid=False, create_volumes=False):
    """Fills the grid_ggd of an IDS with a uniform rectangular grid of size N x N,
    containing vertices, edges and faces, and optionally volumes.

    Adapted from https://imas-data-dictionary.readthedocs.io/en/latest/ggd_guide/examples.html

    Args:
        grid_ggd: The GGD grid that will be filled with the N x N grid.
        N: The size of the N x N grid
        create_3d_grid: If True, create a 2D grid in 3D space (in the X-Z plane at
            Y = 0). If False, create a 2D grid in the X-Y plane.
        create_volumes: If True, extrude the 2D grid by in the Z direction to
            create a layer of cube volume cells. Requires create_3d_grid=True.

    Returns:
        num_vertices: The number of vertices in the generated grid_ggd
        num_edges: The number of edges in the generated grid_ggd
        num_faces: The number of faces in the generated grid_ggd
        num_volumes: The number of volume cells, or 0 if create_volumes is False
    """
    if create_volumes and not create_3d_grid:
        raise RuntimeError("Cannot create volumes if grid is not 3D")

    # Set grid
    grid_ggd.identifier = identifiers.ggd_identifier.linear

    # Set space
    grid_ggd.space.resize(1)
    space = grid_ggd.space[0]
    space.identifier = identifiers.ggd_space_identifier.primary_standard
    space.geometry_type.index = 0  # standard, non-Fourier geometry

    num_dimension = 3 if create_3d_grid else 2
    space.coordinates_type.resize(num_dimension)
    coordinate_identifier = identifiers.coordinate_identifier

    # coordinates_type changed from INT_1D to AoS of identifiers in DD4.0.0
    if isinstance(space.coordinates_type, IDSStructure):
        space.coordinates_type[0] = coordinate_identifier.x
        space.coordinates_type[1] = coordinate_identifier.y
        if create_3d_grid:
            space.coordinates_type[2] = coordinate_identifier.z
    else:
        space.coordinates_type[0] = coordinate_identifier.x.index
        space.coordinates_type[1] = coordinate_identifier.y.index
        if create_3d_grid:
            space.coordinates_type[2] = coordinate_identifier.z.index

    space.objects_per_dimension.resize(4 if create_volumes else 3)
    num_vertices, num_edges, num_faces, num_volumes = build_grid(
        N, space, create_3d_grid, create_volumes
    )

    # Create subsets
    grid_ggd.grid_subset.resize(4 if create_volumes else 3)
    ggd_subset_identifier = identifiers.ggd_subset_identifier
    grid_subsets = grid_ggd.grid_subset
    grid_subsets[0].dimension = 1
    grid_subsets[0].identifier = ggd_subset_identifier.nodes
    grid_subsets[1].dimension = 2
    grid_subsets[1].identifier = ggd_subset_identifier.edges
    grid_subsets[2].dimension = 3
    grid_subsets[2].identifier = ggd_subset_identifier.cells

    # Create elements for nodes
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

    if create_volumes:
        grid_subsets[3].dimension = 4
        grid_subsets[3].identifier = ggd_subset_identifier.volumes
        grid_subsets[3].element.resize(num_volumes)
        for i, element in enumerate(grid_subsets[3].element):
            element.object.resize(1)
            element.object[0].space = 1
            element.object[0].dimension = 4
            element.object[0].index = i + 1

    return num_vertices, num_edges, num_faces, num_volumes


def build_grid(N, space, create_3d_grid, create_volumes):
    """Builds the vertices, edges and faces for an N x N grid.

    Args:
        N: Size of the grid
        space: Space AoS of the GGD grid
        create_3d_grid: If True, create a 2D grid in 3D space (in the X-Z plane at
            Y = 0). If False, create a 2D grid in the X-Y plane.
        create_volumes: If True, extrude the 2D grid by in the Z direction to
            create a layer of cube volume cells. Requires create_3d_grid=True.

    Returns:
        num_vertices: The total number of vertices
        num_edges: The total number of edges
        num_faces: The total number of faces
    """
    if create_volumes and not create_3d_grid:
        raise RuntimeError("Cannot create volumes if grid is not 3D")
    num_vertices = set_vertices(N, space, create_3d_grid, create_volumes)
    num_edges = set_edges(N, space)
    num_faces = set_faces(N, space)

    num_volumes = 0
    if create_volumes:
        num_volumes = set_volumes(N, space)

    return num_vertices, num_edges, num_faces, num_volumes


def set_vertices(N, space, create_3d_grid, create_volumes):
    """Sets the vertices for an N x N grid.

    Args:
        N: Size of the grid
        space: Space AoS of the GGD grid
        create_3d_grid: If True, create a 2D grid in 3D space (in the X-Z plane at
            Y = 0). If False, create a 2D grid in the X-Y plane.
        create_volumes: If True, extrude the 2D grid by in the Z direction.

    Returns:
        num_vertices: The total number of vertices
    """
    if create_volumes and not create_3d_grid:
        raise RuntimeError("Cannot create volumes if grid is not 3D")
    vertices = space.objects_per_dimension[0].object
    base_count = N * N
    total_vertices = base_count * 2 if create_volumes else base_count
    vertices.resize(total_vertices)

    for layer in range(2 if create_volumes else 1):
        z_val = float(layer)
        for i in range(N):
            for j in range(N):
                idx = (layer * base_count) + (i * N + j)
                if create_3d_grid:
                    # Using X-Y plane for the base, Z for extrusion
                    vertices[idx].geometry = [0.5 * float(j), float(i), z_val]
                else:
                    vertices[idx].geometry = [0.5 * float(j), float(i)]
    return total_vertices


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


def set_volumes(N, space):
    """Sets the volumes faces for an N x N x 1 grid.

    Args:
        N: Size of the grid
        space: Space AoS of the GGD grid

    Returns:
        num_volumes: The total number of faces
    """
    volumes = space.objects_per_dimension[3].object
    num_volumes = (N - 1) * (N - 1)
    volumes.resize(num_volumes)

    base_offset = N * N
    vol_idx = 0
    for i in range(N - 1):
        for j in range(N - 1):
            # Bottom layer
            b_tl = i * N + j + 1
            b_tr = b_tl + 1
            b_bl = b_tl + N
            b_br = b_bl + 1

            # Top layer
            t_tl = b_tl + base_offset
            t_tr = b_tr + base_offset
            t_bl = b_bl + base_offset
            t_br = b_br + base_offset

            volumes[vol_idx].nodes = int32array(
                [b_tl, b_tr, b_br, b_bl, t_tl, t_tr, t_br, t_bl]
            )
            vol_idx += 1
    return num_volumes


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


def fill_ids(
    ids,
    time_steps=1,
    grid_size=2,
    fill_ggd=True,
    create_3d_grid=False,
    dynamic_grid_size=False,
):
    """Fills the IDS with an N x N uniform GGD grid and optionally fills all GGD arrays
    on this grid with random values.

    Args:
        ids: IDS to be filled.
        time_steps: Number of time steps to create in the IDS.
        grid_size: Size of the N x N grid. Defaults to 2, meaning a 2 x 2 grid.
        fill_ggd: Whether to fill the GGD arrays on the grid. If set to False, only
            the grid itself will be created.
        create_3d_grid: If True, create a 2D grid in 3D space (in the X-Z plane at
            Y = 0). If False, create a 2D grid in the X-Y plane.
        dynamic_grid_size: If True, increase the grid size from `grid_size` by 1 for
            each subsequent time step. E.g. if `grid_size = 2` and `time_step=3`, then
            grids of sizes 2, 3 x 3, and 4 x 4 will be created.
    """

    # Create an empty grid_ggd
    grid_ggd = create_first_grid(ids)

    # Skip filling grid_ggd if it does not exist
    if grid_ggd is None:
        logger.warning("%s has no grid_ggd", ids.metadata.name)
        return

    # Create time steps
    ids.time = [float(t) for t in range(time_steps)]
    ids.ids_properties.homogeneous_time = imas.ids_defs.IDS_TIME_MODE_HOMOGENEOUS

    # Create grid and GGD AoS
    grid_ggd_aos = imas.util.get_parent(grid_ggd)

    grid_ggd_aos.resize(time_steps)

    # Create uniform grids for each time step
    for i in range(time_steps):
        num_vertices, num_edges, num_faces, _ = fill_NxN_grid(
            grid_ggd_aos[i],
            grid_size,
            create_3d_grid=create_3d_grid,
            create_volumes=False,
        )
        if dynamic_grid_size:
            grid_size += 1
        logger.debug("filled grid_ggd at index %d.", i)

    if fill_ggd:
        ggd = create_first_ggd(ids)
        ggd_aos = imas.util.get_parent(ggd)
        ggd_aos.resize(time_steps)
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
