"""
A collection of methods to read a grid_ggd IDS node into a VTK dataset.
These methods copy contents from the grid_ggd/space and grid_ggd/grid_subset
children into distinct vtkUnstructuredGrid objects for Bezier elements.
"""
from typing import Callable, Any

from vtkmodules.vtkCommonCore import vtkPoints,  vtkIdList
from vtkmodules.vtkCommonDataModel import ( vtkUnstructuredGrid, vtkCellArray,
    VTK_BEZIER_QUADRILATERAL, VTK_BEZIER_HEXAHEDRON )
from vtkmodules.util.vtkConstants import VTK_FLOAT, VTK_ID_TYPE
from vtkmodules.util import numpy_support as npvtk

from paraview import logger as pvlog

import numpy as np
import vtk
import operator

# For the units: 
from vtkggdtools.imashelper import get_units

prec = np.float32
vtk_prec = VTK_FLOAT

def convert_grid_subset_to_unstructured_grid(ids_name: str, ids, aos_index_values: dict, n_plane: int, phi_start, phi_end: float) -> vtkUnstructuredGrid:
    """
    Copy the elements found in given grid_ggd/grid_subset IDS node into a vtkUnstructuredGrid instance.
    This method uses the supplied point coordinates in the form of a vtkPoints instance.
    :param grid_ggd: a grid_ggd ids node
    :param subset_idx: an index into grid_ggd/grid_subset
    :param n_plane: number of toroidal planes to be generated if 3D axysimetric 
    :param phi_start: start phi plane
    :param phi_end: end plane at phi in degrees
    :param vtk_grid_points: the point coordinates corresponding to 1d objects in the subset elements.
    :return:
    """
    output = vtkUnstructuredGrid()

    time_idx = aos_index_values.get('TimeIdx')
    try:
        if ids_name == 'mhd':
            ggd = ids.ggd[time_idx]
            mhdval = True
        elif ids_name == 'radiation':
            ggd = ids.grid_ggd[time_idx]
            radval = True
        else:
            raise IndexError
    except IndexError:
        return output

    phi = [phi_start, phi_end]
    gr2d = ids.grid_ggd[0].space[0]
    xyz0 = gr2d.objects_per_dimension[0].object.array
    ien0 = np.array(gr2d.objects_per_dimension[2].object.array)

    n_period = ids.grid_ggd[0].space[1].geometry_type.index

    x = np.zeros((2, 4, len(xyz0)))
    for j in range(np.shape(x)[2]):
        x[:, :, j] = gr2d.objects_per_dimension.array[0].object.array[j].geometry_2d

    # size 1, d_{uk}, d_{vk}, d{uv}d{vk} as in Daan Van Vugt thesis
    size = np.empty((4, 4, np.shape(ien0)[0]))
    for i in range(np.shape(size)[2]):
        size[:, :, i] = gr2d.objects_per_dimension.array[2].object.array[i].geometry_2d

    val_tor1 = np.array([])
    nam = list()
    
    if mhdval:
        data = ids.ggd.array[time_idx]
        ggd_path = 'ggd'

        quantities = {
            't_i_average': 'Ion Temperature (average)',
            'n_i_total': 'Ion Density (total)',
            'zeff': 'Z effective',
            'b_field_r': 'Magnetic Field Br',
            'b_field_z': 'Magnetic Field Bz',
            'b_field_tor': 'Magnetic Field Btor',
            'a_field_r': 'Magnetic Potential Ar',
            'a_field_z': 'Magnetic Potential Az',
            'a_field_tor': 'Magnetic Potential Ator',
            'psi': 'Poloidal Flux',
            'velocity_r': 'Plasma Velocity Vr',
            'velocity_z': 'Plasma Velocity Vz',
            'velocity_tor': 'Plasma Velocity Vtor',
            'velocity_parallel': 'Plasma Velocity Vparallel',
            'phi_potential': 'Electric Potential',
            'vorticity': 'Vorticity',
            'j_r': 'Current Density Jr',
            'j_z': 'Current Density Jz',
            'j_tor': 'Current Density Jtor',
            'mass_density': 'Mass Density'
        }
        for q_field, q_name in quantities.items():
            val_tor1, nam = value_in_IDS(ids_name, ggd_path, data, q_field, q_name, val_tor1, nam)
        

    elif radval:
        ggd_path = 'process/ggd'
        data = ids.process[0].ggd.array[time_idx]
        
        try:
            if np.size(val_tor1) == 0:
                val_tor1 = np.array([data_r.ion[0].emissivity[0].coefficients])
                nam.append("Radiation (W/m^3)")
            else:
                val_tor1 = np.concatenate((val_tor1,
                    np.array([data_r.ion[0].emissivity[0].coefficients])), axis=0)
                nam.append("Radiation (W/m^3)")
        except:
            print('No radiation values')

    else:
        print('No mhd or radiation values found')
        return output

    n_val = len(nam)
    a = np.shape(val_tor1)
    n_tor = int(a[1]/len(xyz0))
    valu = np.reshape(val_tor1, (n_val, n_tor, len(xyz0), 4))
    values = np.swapaxes(valu, 2, 3)
    values = np.swapaxes(values, 1, 2)

            #vertex
    ien0 = np.array(ids.grid_ggd.array[0].space.array[0].objects_per_dimension.array[2].object.array)
    ver = np.empty((np.shape(ien0)[0], np.shape(ien0[0].nodes)[0]))
    for i in range(np.shape(ien0)[0]):
        ver[i, :] = np.array(ien0[i].nodes)
    ver = ver.astype(int)
    vertex = np.swapaxes(ver, 1, 0)

    # Everything we need to visualise data is now excavated from IDS file
    n_plane = 1 + (n_plane - 1) * 2
    n_sub = 4
    without_n0_mode = False
    periodic = False
            
    if (n_plane == 1):
        phis = np.asarray([phi[0]])
    else:
        periodic = (np.mod(phi[0] - phi[1], 360) < 1e-9)
        phis = np.linspace(phi[0], phi[1], num=n_plane, endpoint=not periodic)

    phis = phis * np.pi / 180
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
        xy[16 * i:16 * (i + 1), 0] = np.ravel(x1)
        xy[16 * i:16 * (i + 1), 1] = np.ravel(y1)
    RZ = xy

    n_xy = np.shape(RZ)[0]
    xyz = np.zeros((n_xy * n_plane, 3))
    for i in range(n_plane):
        xyz[i * n_xy:(i + 1) * n_xy, 0] = np.ravel(RZ[:, 0] * np.cos(phis[i]))
        xyz[i * n_xy:(i + 1) * n_xy, 1] = np.ravel(RZ[:, 1])
        xyz[i * n_xy:(i + 1) * n_xy, 2] = np.ravel(RZ[:, 0] * np.sin(phis[i]))

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
        s = np.shape(xyz)[0] // n_plane
        w = np.cos(np.deg2rad(alpha))
        w1 = np.ones((np.shape(xyz)[0]))
        ien = None
        index = np.array([0, 1, 2, 3, 0 + 2 * s, 1 + 2 * s, 2 + 2 * s, 3 + 2 * s, 4, 5, 9, 10, 7, 6, 8, 11,
                          4 + 2 * s, 5 + 2 * s, 9 + 2 * s, 10 + 2 * s, 7 + 2 * s, 6 + 2 * s, 8 + 2 * s,
                          11 + 2 * s, s, 1 + s, 3 + s, 2 + s, 8 + s, 11 + s, 9 + s, 10 + s, 4 + s, 5 + s,
                          7 + s, 6 + s, 12, 13, 15, 14, 12 + 2 * s, 13 + 2 * s, 15 + 2 * s, 14 + 2 * s,
                          12 + s, 13 + s, 15 + s, 14 + s])
        for i in range((n_plane - 1) // 2):
            index2 = index + i * np.array([s * 2 for k in range(48)])
            w1[s + i * 2 * s: 2 * s + i * 2 * s] = np.array([w for i in range(s)])
            xyz[s + i * 2 * s: 2 * s + i * 2 * s, 0] = 1 / w * xyz[s + i * 2 * s: 2 * s + i * 2 * s, 0]
            xyz[s + i * 2 * s: 2 * s + i * 2 * s, 2] = 1 / w * xyz[s + i * 2 * s: 2 * s + i * 2 * s, 2]
            step = np.array([16 for i in range(48)])
            ien2 = np.zeros((np.shape(xyz)[0] // n_plane // 16, 48))
            for j in range(0, np.shape(xyz)[0] // n_plane // 16):
                ien2[j, :] = np.array([index2 + j * step])

            if (np.any(ien)):
                ien = np.concatenate((ien, ien2), axis=0)
            else:
                ien = ien2

        ien = np.insert(ien, 0, 48, axis=1)
        etype = vtk.VTK_BEZIER_HEXAHEDRON

        weights = npvtk.numpy_to_vtk(w1, deep=True, array_type=vtk_prec)
        weights.SetName("RationalWeights")
        output.GetPointData().SetRationalWeights(weights)
        n_c = int(np.shape(xyz)[0] / n_plane / 16 * (n_plane - 1) / 2)
        degrees = npvtk.numpy_to_vtk(np.array([[3, 3, 2] for k in range(n_c)]), deep=True,
                                     array_type=vtk.VTK_ID_TYPE)
        degrees.SetName("HigherOrderDegrees")

        output.GetCellData().SetHigherOrderDegrees(degrees)

    #Use for changing grid orientation in paraview
    z = np.copy(xyz[:,1])
    y = np.copy(xyz[:,2])
    xyz[:,1] = y
    xyz[:,2] = z

    pcoords = npvtk.numpy_to_vtk(xyz, deep=True, array_type=vtk_prec)
    points = vtk.vtkPoints()
    points.SetData(pcoords)

    cells = vtk.vtkCellArray()
    cells.SetCells(ien.shape[0], npvtk.numpy_to_vtk(ien, deep=True, array_type=vtk.VTK_ID_TYPE))

    output.SetPoints(points)
    output.SetCells(etype, cells)
            
    HZ = toroidal_basis(n_tor, n_period, phis, without_n0_mode)

    val = interp_scalars_3D(values, vertex, size, n_sub, HZ).reshape((n_val, -1))

    a = np.shape(val)
    val = val.reshape((a[0], a[1] // 16, 16))
    val[:, :, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]] = val[:, :, [0, 12, 15, 3, 4, 8, 11, 7, 1,
                                                                                           13, 14, 2, 5, 9, 10, 6]]
    val = val.reshape((a[0], a[1]))

    for i in range(len(nam)):
        tmp = npvtk.numpy_to_vtk(val[i, :], deep=True, array_type=vtk_prec)
        tmp.SetName(nam[i])
        output.GetPointData().AddArray(tmp)
    time2 = ids.time[time_idx]
    stime = npvtk.numpy_to_vtk(np.array([time2]), deep=True, array_type=vtk_prec)
    stime.SetName("TimeValue")
    output.GetFieldData().AddArray(stime)

    return output

def toroidal_basis(n_tor, n_period, phis, without_n0_mode):
    # Setup toroidal coefficients for each plane and toroidal harmonic
    HZ = np.zeros((n_tor,len(phis)))
    for i in range(n_tor):
        mode = np.floor((i+1)/2)*n_period
        if (i == 0):
            if (not without_n0_mode):
                HZ[i,:] = 1
        elif (i % 2 == 0):
            HZ[i,:] = np.sin(mode*phis)
        elif (i % 2 == 1):
            HZ[i,:] = np.cos(mode*phis)
    return HZ
"""
Interpolate scalars on 2D poloidal plane

returns:
    values: interpolated values, values[var, harmonic, element, is, it]
"""
def interp_scalars(values, vertex, size, n_sub):
    # Multiply values[var,order,harm,vertex,element] with
    # size[order, vertex, element] and bf[order, vertex, s, t]
    return np.einsum('lihjk,ijk,ijmn->lhkmn',
                    values[:,:,:,vertex-1],
                    size, bf(n_sub))
"""
Interpolate scalars on 2D planes * n_planes
"""
def interp_scalars_3D(values, vertex, size, n_sub, HZ):
    vals = interp_scalars(values, vertex, size, n_sub)
    return np.einsum('lhkmn,hp->lpkmn', vals, HZ)

"""
Calculate values of the basis functions at positions s and t
Optionally put many values of s and t at once as numpy arrays.
Dimension 0: order
Dimension 1: vertex
optional dimension 2, 3: position s, t
"""
def basis_functions(s,t):
    return np.asarray([
        [ (-1 + s)**2*(1 + 2*s)*(-1 + t)**2*(1 + 2*t),
         -(s**2*(-3 + 2*s)*(-1 + t)**2*(1 + 2*t)),
          s**2*(-3 + 2*s)*t**2*(-3 + 2*t),
         -(-1 + s)**2*(1 + 2*s)*t**2*(-3 + 2*t)],
        [ 3*(-1 + s)**2*s*(-1 + t)**2*(1 + 2*t),
         -3*(-1 + s)*s**2*(-1 + t)**2*(1 + 2*t),
         3*(-1 + s)*s**2*t**2*(-3 + 2*t),
         -3*(-1 + s)**2*s*t**2*(-3 + 2*t)],
        [ 3*(-1 + s)**2*(1 + 2*s)*(-1 + t)**2*t,
         -3*s**2*(-3 + 2*s)*(-1 + t)**2*t,
          3*s**2*(-3 + 2*s)*(-1 + t)*t**2,
         -3*(-1 + s)**2*(1 + 2*s)*(-1 + t)*t**2],
        [ 9*(-1 + s)**2*s*(-1 + t)**2*t,
         -9*(-1 + s)*s**2*(-1 + t)**2*t,
          9*(-1 + s)*s**2*(-1 + t)*t**2,
         -9*(-1 + s)**2*s*(-1 + t)*t**2]])

"""
Calculate values of the basis functions derived to s at positions s and t
Optionally put many values of s and t at once as numpy arrays.
Dimension 0: order
Dimension 1: vertex
optional dimension 2, 3: position s, t
"""
def basis_functions_s(s,t):
    return np.asarray([
        [ 6*(-1 + s)*s*(-1 + t)**2*(1 + 2*t),
         -6*(-1 + s)*s*(-1 + t)**2*(1 + 2*t),
          6*(-1 + s)*s*t**2*(-3 + 2*t),
         -6*(-1 + s)**2*t**2*(-3 + 2*t)],
        [ 3*(-1 + s)*(-1+3*s)*(-1 + t)**2*(1 + 2*t),
         -3*s*(-2 + 3*s)*(-1 + t)**2*(1 + 2*t),
          3*s*(-2 + 3*s)*t**2*(-3 + 2*t),
         -3*(-1 + 3*s)*(-1 + s)*t**2*(-3 + 2*t)],
        [ 18*(-1 + s)*s*(-1 + t)**2*t,
         -18*(-1 + s)*s*(-1 + t)**2*t,
          18*(-1 + s)*s*(-1 + t)*t**2,
         -18*(-1 + s)*s*(-1 + t)*t**2],
        [ 9*(-1 + s)*(-1+3*s)*(-1 + t)**2*t,
         -9*s*(-2 + 3*s)*(-1 + t)**2*t,
          9*s*(-2 + 3*s)*(-1 + t)*t**2,
         -9*(-1 + 3*s)*(-1 + s)*(-1 + t)*t**2]])

"""
Calculate values of the basis functions derived to t at positions s and t
Optionally put many values of s and t at once as numpy arrays.
Dimension 0: order
Dimension 1: vertex
optional dimension 2, 3: position s, t
"""
def basis_functions_t(s,t):
    return np.asarray([
        [ 6*(-1 + s)**2*(1 + 2*s)*(-1 + t)*t,
         -6*s**2*(-3 + 2*s)*(-1 + t)*t,
          6*s**2*(-3 + 2*s)*(-1 + t)*t,
         -6*(-1 + s)**2*(1 + 2*s)*(-1 + t)*t],
        [ 18*(-1 + s)**2*s*(-1 + t)*t,
         -18*(-1 + s)*s**2*(-1 + t)*t,
          18*(-1 + s)*s**2*(-1 + t)*t,
         -18*(-1 + s)**2*s*(-1 + t)*t],
        [ 3*(-1 + s)**2*(1 + 2*s)*(-1 + t)*(-1 + 3*t),
         -3*s**2*(-3 + 2*s)*(1 - 3*t)*(-1 + t),
          3*s**2*(-3 + 2*s)*t*(-2 + 3*t),
         -3*(-1 + s)**2*(1 + 2*s)*t*(-2 + 3*t)],
        [ 9*(-1 + s)**2*s*(-1 + t)*(-1 + 3*t),
         -9*(-1 + s)*s**2*(-1 + t)*(-1 + 3*t),
          9*(-1 + s)*s**2*t*(-2 + 3*t),
         -9*(-1 + s)**2*s*t*(-2 + 3*t)]])


"""
Calculate basis functions at n_sub**2 points
"""
def bf(n_sub):
    # Get the basis functions at each of the points
    lin = np.linspace(0.0, 1.0, n_sub)
    s  = np.tensordot(lin, [1]*n_sub, axes=0)
    t  = s.transpose()
    return basis_functions(s, t)
def bf_s(n_sub):
    # Get the basis functions at each of the points
    lin = np.linspace(0.0, 1.0, n_sub)
    s  = np.tensordot(lin, [1]*n_sub, axes=0)
    t  = s.transpose()
    return basis_functions_s(s, t)
def bf_t(n_sub):
    # Get the basis functions at each of the points
    lin = np.linspace(0.0, 1.0, n_sub)
    s  = np.tensordot(lin, [1]*n_sub, axes=0)
    t  = s.transpose()
    return basis_functions_t(s, t)

def value_in_IDS(ids_name, ggd_path, ids_data, field:str, name:str, valu:np.ndarray, names:list):
    """
    Check if IDS contains chosen value and add it to array of all values
    """
    try:
        new_val = operator.attrgetter(field)(ids_data).array[0].coefficients
        #pvlog.info(f'Getting units for {ids_name}/{ggd_path}/{field}')
        units = get_units(ids_name, f'{ggd_path}/{field}')
        name = f"{name} {units}"
        if not(names):
            names = list()
            names.append(name)
        else:
            names.append(name)
    
        if np.size(valu) == 0:
            valu = np.array([new_val])
            return valu, names

        valu = np.concatenate((valu, np.array([new_val])), axis=0)
        return valu, names

    except Exception as e:
        #pvlog.exception(e)
        return valu, names
