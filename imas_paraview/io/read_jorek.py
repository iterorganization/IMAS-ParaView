"""
A collection of methods to read a grid_ggd IDS node into a VTK dataset.
These methods copy contents from the grid_ggd/space and grid_ggd/grid_subset
children into distinct vtkUnstructuredGrid objects for Bezier elements.
"""

import logging

import numpy as np
import vtk
from vtkmodules.util import numpy_support as npvtk
from vtkmodules.util.vtkConstants import VTK_DOUBLE
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

logger = logging.getLogger("imas_paraview")

prec = np.float64
vtk_prec = VTK_DOUBLE


def read_plasma_state(grid_ggd, ps_reader, plane_config, ugrid):
    """Reads plasma state data arrays from the ggd node in the IDS. These arrays are
        added as point data to the unstructured grid.

    Args:
        grid_ggd: The GGD grid.
        ps_reader: PlasmaStateReader which contains which array to load
        plane_config: Configuration for the bezier planes
        ugrid: The ugrid to load the GGD arrays for
    """
    n_plane = plane_config.n_plane

    phi = [plane_config.phi_start, plane_config.phi_end]
    val_tor1 = np.array([])
    nam = list()

    N_vertex = len(grid_ggd.space[0].objects_per_dimension[0].object)
    array_list = ps_reader.scalar_array_list + ps_reader.vector_array_list

    # Load selected arrays
    for attribute_array in array_list:
        name = ps_reader._create_name_with_units(attribute_array)
        if hasattr(attribute_array[0], "coefficients"):
            scalar_data = attribute_array[0].coefficients
        elif hasattr(attribute_array[0], "phi_coefficients"):
            scalar_data = attribute_array[0].phi_coefficients
        else:
            logger.warning(f"Could not load the coefficients for {name}, it is ignored")
            continue
        nam.append(name)
        if np.size(val_tor1) == 0:
            val_tor1 = np.array([scalar_data])
        else:
            val_tor1 = np.concatenate((val_tor1, np.array([scalar_data])), axis=0)

    n_val = len(nam)
    a = np.shape(val_tor1)
    n_tor = a[1] // N_vertex if len(a) > 1 else 0

    valu = np.reshape(val_tor1, (n_val, n_tor, N_vertex, 4))
    values = np.swapaxes(valu, 2, 3)
    values = np.swapaxes(values, 1, 2)

    n_period = grid_ggd.space[1].geometry_type.index
    n_sub = 4
    without_n0_mode = False

    n_plane = 1 + (n_plane - 1) * 2
    phis = calc_phis(n_plane, phi)

    HZ = toroidal_basis(n_tor, n_period, phis, without_n0_mode)

    vertex = get_vertex(grid_ggd)
    size = calc_size(grid_ggd)
    val = interp_scalars_3D(values, vertex, size, n_sub, HZ)
    if n_val > 0:
        val = val.reshape((n_val, -1))
    a = np.shape(val)
    val = val.reshape((a[0], a[1] // 16, 16))
    val[:, :, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]] = val[
        :, :, [0, 12, 15, 3, 4, 8, 11, 7, 1, 13, 14, 2, 5, 9, 10, 6]
    ]
    val = val.reshape((a[0], a[1]))

    for i in range(len(nam)):
        tmp = npvtk.numpy_to_vtk(val[i, :], deep=True, array_type=vtk_prec)
        tmp.SetName(nam[i])
        ugrid.GetPointData().AddArray(tmp)


def convert_grid_subset_to_unstructured_grid(
    grid_ggd, plane_config
) -> vtkUnstructuredGrid:
    """
    Copy the elements found in given grid_ggd/grid_subset IDS node into a
    vtkUnstructuredGrid instance. This method uses the supplied point coordinates in
    the form of a vtkPoints instance.

    Args:
        grid_ggd: The GGD grid.
        plane_config: Configuration for the bezier planes.

    Returns:
        The created unstructured grid.
    """
    output = vtkUnstructuredGrid()

    n_plane = plane_config.n_plane
    phi = [plane_config.phi_start, plane_config.phi_end]
    gr2d = grid_ggd.space[0]
    N_vertex = len(gr2d.objects_per_dimension[0].object)

    x = np.zeros((2, 4, N_vertex))
    for j, obj in enumerate(gr2d.objects_per_dimension[0].object):
        x[:, :, j] = obj.geometry_2d.value

    size = calc_size(grid_ggd)
    vertex = get_vertex(grid_ggd)

    # Everything we need to visualise data is now excavated from IDS file
    n_plane = 1 + (n_plane - 1) * 2

    phis = calc_phis(n_plane, phi)
    ien = None

    tmp = np.zeros((x.shape[0], x.shape[1], vertex.shape[0], vertex.shape[1]))
    for i in range(vertex.shape[0]):  # small loop over vertices (hardcode 4 here?)
        tmp[:, :, i, :] = x[:, :, vertex[i, :] - 1]
    # multiply by size[order, vertex, element]
    tmp[0, :, :, :] *= size
    tmp[1, :, :, :] *= size
    # Create output array
    xy = np.zeros((np.shape(tmp)[3] * 16, 2))
    for i in range(np.shape(tmp)[3]):
        x1 = tmp[0, :, :, i]
        y1 = tmp[1, :, :, i]
        x1[3, :] = x1[3, :] + x1[1, :] + x1[2, :]
        y1[3, :] = y1[3, :] + y1[1, :] + y1[2, :]
        x1[1:, :] += np.tile(x1[0, :], (3, 1))
        y1[1:, :] += np.tile(y1[0, :], (3, 1))
        xy[16 * i : 16 * (i + 1), 0] = np.ravel(x1)
        xy[16 * i : 16 * (i + 1), 1] = np.ravel(y1)
    RZ = xy

    n_xy = np.shape(RZ)[0]
    xyz = np.zeros((n_xy * n_plane, 3))
    for i in range(n_plane):
        # x = R * cos(phi); y = R * sin(phi), z = Z
        xyz[i * n_xy : (i + 1) * n_xy, 0] = np.ravel(RZ[:, 0] * np.cos(phis[i]))
        xyz[i * n_xy : (i + 1) * n_xy, 1] = np.ravel(RZ[:, 0] * np.sin(phis[i]))
        xyz[i * n_xy : (i + 1) * n_xy, 2] = np.ravel(RZ[:, 1])

    if n_plane == 1:
        index = np.array([0, 1, 2, 3, 4, 5, 9, 10, 7, 6, 8, 11, 12, 13, 15, 14])
        step = np.array([16 for i in range(16)])
        ien2 = np.array([index])
        for i in range(1, np.shape(xyz)[0] // 16):
            ien2 = np.concatenate((ien2, np.array([index + i * step])), axis=0)
        ien = np.insert(ien2, 0, 16, axis=1)
        etype = vtk.VTK_BEZIER_QUADRILATERAL

    else:
        alpha = (phi[1] - phi[0]) / (n_plane - 1)
        w = np.cos(np.deg2rad(alpha))
        w1 = np.ones((np.shape(xyz)[0]))
        ien = None
        # Connectivity list of one bezier cell:
        index = np.array(
            [
                0,
                1,
                2,
                3,
                0 + 2 * n_xy,
                1 + 2 * n_xy,
                2 + 2 * n_xy,
                3 + 2 * n_xy,
                4,
                5,
                9,
                10,
                7,
                6,
                8,
                11,
                4 + 2 * n_xy,
                5 + 2 * n_xy,
                9 + 2 * n_xy,
                10 + 2 * n_xy,
                7 + 2 * n_xy,
                6 + 2 * n_xy,
                8 + 2 * n_xy,
                11 + 2 * n_xy,
                # According to VTK docs, the sequence should be ([0, 1, 3, 2] + n_xy).
                # However that doesn't render correctly and we need to swap the 3 and
                # the 2... See also
                # https://www.kitware.com/main/wp-content/uploads/2020/03/Implementation-of-rational-Be%CC%81zier-cells-into-VTK-Report.pdf
                0 + n_xy,
                1 + n_xy,
                2 + n_xy,
                3 + n_xy,
                8 + n_xy,
                11 + n_xy,
                9 + n_xy,
                10 + n_xy,
                4 + n_xy,
                5 + n_xy,
                7 + n_xy,
                6 + n_xy,
                12,
                13,
                15,
                14,
                12 + 2 * n_xy,
                13 + 2 * n_xy,
                15 + 2 * n_xy,
                14 + 2 * n_xy,
                12 + n_xy,
                13 + n_xy,
                15 + n_xy,
                14 + n_xy,
            ]
        )
        # Connectivity for all cells between two planes
        repeat_in_plane = 16 * np.arange(n_xy // 16)
        repeat_planes = np.arange((n_plane - 1) // 2) * 2 * n_xy
        repeat_index = (repeat_in_plane + repeat_planes[:, np.newaxis]).ravel()
        ien = index + repeat_index[:, np.newaxis]

        # set rational weights of midplane points and adjust x/y values accordingly
        for i in range((n_plane - 1) // 2):
            w1[n_xy * (1 + 2 * i) : n_xy * (2 + 2 * i)] = w
            xyz[n_xy * (1 + 2 * i) : n_xy * (2 + 2 * i), 0:2] /= w

        # VTK cell list requires number of points as first element:
        ien = np.insert(ien, 0, 48, axis=1)
        etype = vtk.VTK_BEZIER_HEXAHEDRON

        weights = npvtk.numpy_to_vtk(w1, deep=True, array_type=vtk_prec)
        weights.SetName("RationalWeights")
        output.GetPointData().SetRationalWeights(weights)

        # Bezier degrees of the cells
        degrees = npvtk.numpy_to_vtk(
            np.repeat([[3, 3, 2]], len(ien), axis=0),
            deep=True,
            array_type=vtk.VTK_ID_TYPE,
        )
        degrees.SetName("HigherOrderDegrees")
        output.GetCellData().SetHigherOrderDegrees(degrees)

    pcoords = npvtk.numpy_to_vtk(xyz, deep=True, array_type=vtk_prec)
    points = vtk.vtkPoints()
    points.SetData(pcoords)

    cells = vtk.vtkCellArray()
    cells.SetCells(
        ien.shape[0], npvtk.numpy_to_vtk(ien, deep=True, array_type=vtk.VTK_ID_TYPE)
    )

    output.SetPoints(points)
    output.SetCells(etype, cells)

    return output


def calc_size(grid_ggd):
    """Calculate the size matrix for each face in the grid.
    Uses the geometry values from the 2D grid representation.

    Args:
        grid_ggd: The GGD grid.

    Returns:
        The calculated size matrix.
    """
    # size 1, d_{uk}, d_{vk}, d{uv}d{vk} as in Daan Van Vugt thesis
    gr2d = grid_ggd.space[0]
    n_face = len(gr2d.objects_per_dimension[2].object)
    size = np.empty((4, 4, n_face))

    for i, obj in enumerate(gr2d.objects_per_dimension[2].object):
        size[:, :, i] = obj.geometry_2d.value
    return size


def get_vertex(grid_ggd):
    """Retrieve vertex numbers in the R-Z plane.
    Extracts node indices for each face in the grid.

    Args:
        grid_ggd: The GGD grid.

    Returns:
        The vertex numbers.
    """
    # Get vertex numbers in the R,Z plane
    object_aos = grid_ggd.space[0].objects_per_dimension[2].object
    vertex = np.zeros((len(object_aos[0].nodes), len(object_aos)), dtype=np.int32)
    for i, obj in enumerate(object_aos):
        vertex[:, i] = obj.nodes.value
    return vertex


def calc_phis(n_plane, phi):
    """Compute an array of phi values for given planes.
    Determines if the values are periodic and converts them to radians.

    Args:
        n_plane: The number of planes.
        phi: Array containing the start and end phi values.

    Returns:
        the array of phi values for the planes.
    """
    periodic = False
    if n_plane == 1:
        phis = np.asarray([phi[0]])
    else:
        periodic = np.mod(phi[0] - phi[1], 360) < 1e-9
        phis = np.linspace(phi[0], phi[1], num=n_plane, endpoint=not periodic)

    phis = phis * np.pi / 180
    return phis


def toroidal_basis(n_tor, n_period, phis, without_n0_mode):
    # Setup toroidal coefficients for each plane and toroidal harmonic
    HZ = np.zeros((n_tor, len(phis)))
    for i in range(n_tor):
        mode = np.floor((i + 1) / 2) * n_period
        if i == 0:
            if not without_n0_mode:
                HZ[i, :] = 1
        elif i % 2 == 0:
            HZ[i, :] = np.sin(mode * phis)
        elif i % 2 == 1:
            HZ[i, :] = np.cos(mode * phis)
    return HZ


def interp_scalars(values, vertex, size, n_sub):
    """
    Interpolate scalars on 2D poloidal plane

    returns:
        values: interpolated values, values[var, harmonic, element, is, it]
    """
    # Multiply values[var,order,harm,vertex,element] with
    # size[order, vertex, element] and bf[order, vertex, s, t]
    return np.einsum(
        "lihjk,ijk,ijmn->lhkmn", values[:, :, :, vertex - 1], size, bf(n_sub)
    )


def interp_scalars_3D(values, vertex, size, n_sub, HZ):
    """
    Interpolate scalars on 2D planes * n_planes
    """
    vals = interp_scalars(values, vertex, size, n_sub)
    return np.einsum("lhkmn,hp->lpkmn", vals, HZ)


def basis_functions(s, t):
    """
    Calculate values of the basis functions at positions s and t
    Optionally put many values of s and t at once as numpy arrays.
    Dimension 0: order
    Dimension 1: vertex
    optional dimension 2, 3: position s, t
    """
    return np.asarray(
        [
            [
                (-1 + s) ** 2 * (1 + 2 * s) * (-1 + t) ** 2 * (1 + 2 * t),
                -(s**2 * (-3 + 2 * s) * (-1 + t) ** 2 * (1 + 2 * t)),
                s**2 * (-3 + 2 * s) * t**2 * (-3 + 2 * t),
                -((-1 + s) ** 2) * (1 + 2 * s) * t**2 * (-3 + 2 * t),
            ],
            [
                3 * (-1 + s) ** 2 * s * (-1 + t) ** 2 * (1 + 2 * t),
                -3 * (-1 + s) * s**2 * (-1 + t) ** 2 * (1 + 2 * t),
                3 * (-1 + s) * s**2 * t**2 * (-3 + 2 * t),
                -3 * (-1 + s) ** 2 * s * t**2 * (-3 + 2 * t),
            ],
            [
                3 * (-1 + s) ** 2 * (1 + 2 * s) * (-1 + t) ** 2 * t,
                -3 * s**2 * (-3 + 2 * s) * (-1 + t) ** 2 * t,
                3 * s**2 * (-3 + 2 * s) * (-1 + t) * t**2,
                -3 * (-1 + s) ** 2 * (1 + 2 * s) * (-1 + t) * t**2,
            ],
            [
                9 * (-1 + s) ** 2 * s * (-1 + t) ** 2 * t,
                -9 * (-1 + s) * s**2 * (-1 + t) ** 2 * t,
                9 * (-1 + s) * s**2 * (-1 + t) * t**2,
                -9 * (-1 + s) ** 2 * s * (-1 + t) * t**2,
            ],
        ]
    )


def basis_functions_s(s, t):
    """
    Calculate values of the basis functions derived to s at positions s and t
    Optionally put many values of s and t at once as numpy arrays.
    Dimension 0: order
    Dimension 1: vertex
    optional dimension 2, 3: position s, t
    """
    return np.asarray(
        [
            [
                6 * (-1 + s) * s * (-1 + t) ** 2 * (1 + 2 * t),
                -6 * (-1 + s) * s * (-1 + t) ** 2 * (1 + 2 * t),
                6 * (-1 + s) * s * t**2 * (-3 + 2 * t),
                -6 * (-1 + s) ** 2 * t**2 * (-3 + 2 * t),
            ],
            [
                3 * (-1 + s) * (-1 + 3 * s) * (-1 + t) ** 2 * (1 + 2 * t),
                -3 * s * (-2 + 3 * s) * (-1 + t) ** 2 * (1 + 2 * t),
                3 * s * (-2 + 3 * s) * t**2 * (-3 + 2 * t),
                -3 * (-1 + 3 * s) * (-1 + s) * t**2 * (-3 + 2 * t),
            ],
            [
                18 * (-1 + s) * s * (-1 + t) ** 2 * t,
                -18 * (-1 + s) * s * (-1 + t) ** 2 * t,
                18 * (-1 + s) * s * (-1 + t) * t**2,
                -18 * (-1 + s) * s * (-1 + t) * t**2,
            ],
            [
                9 * (-1 + s) * (-1 + 3 * s) * (-1 + t) ** 2 * t,
                -9 * s * (-2 + 3 * s) * (-1 + t) ** 2 * t,
                9 * s * (-2 + 3 * s) * (-1 + t) * t**2,
                -9 * (-1 + 3 * s) * (-1 + s) * (-1 + t) * t**2,
            ],
        ]
    )


def basis_functions_t(s, t):
    """
    Calculate values of the basis functions derived to t at positions s and t
    Optionally put many values of s and t at once as numpy arrays.
    Dimension 0: order
    Dimension 1: vertex
    optional dimension 2, 3: position s, t
    """
    return np.asarray(
        [
            [
                6 * (-1 + s) ** 2 * (1 + 2 * s) * (-1 + t) * t,
                -6 * s**2 * (-3 + 2 * s) * (-1 + t) * t,
                6 * s**2 * (-3 + 2 * s) * (-1 + t) * t,
                -6 * (-1 + s) ** 2 * (1 + 2 * s) * (-1 + t) * t,
            ],
            [
                18 * (-1 + s) ** 2 * s * (-1 + t) * t,
                -18 * (-1 + s) * s**2 * (-1 + t) * t,
                18 * (-1 + s) * s**2 * (-1 + t) * t,
                -18 * (-1 + s) ** 2 * s * (-1 + t) * t,
            ],
            [
                3 * (-1 + s) ** 2 * (1 + 2 * s) * (-1 + t) * (-1 + 3 * t),
                -3 * s**2 * (-3 + 2 * s) * (1 - 3 * t) * (-1 + t),
                3 * s**2 * (-3 + 2 * s) * t * (-2 + 3 * t),
                -3 * (-1 + s) ** 2 * (1 + 2 * s) * t * (-2 + 3 * t),
            ],
            [
                9 * (-1 + s) ** 2 * s * (-1 + t) * (-1 + 3 * t),
                -9 * (-1 + s) * s**2 * (-1 + t) * (-1 + 3 * t),
                9 * (-1 + s) * s**2 * t * (-2 + 3 * t),
                -9 * (-1 + s) ** 2 * s * t * (-2 + 3 * t),
            ],
        ]
    )


"""
Calculate basis functions at n_sub**2 points
"""


def bf(n_sub):
    # Get the basis functions at each of the points
    lin = np.linspace(0.0, 1.0, n_sub)
    s = np.tensordot(lin, [1] * n_sub, axes=0)
    t = s.transpose()
    return basis_functions(s, t)


def bf_s(n_sub):
    # Get the basis functions at each of the points
    lin = np.linspace(0.0, 1.0, n_sub)
    s = np.tensordot(lin, [1] * n_sub, axes=0)
    t = s.transpose()
    return basis_functions_s(s, t)


def bf_t(n_sub):
    # Get the basis functions at each of the points
    lin = np.linspace(0.0, 1.0, n_sub)
    s = np.tensordot(lin, [1] * n_sub, axes=0)
    t = s.transpose()
    return basis_functions_t(s, t)
