import numpy as np
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonCore import vtkDoubleArray
from vtkmodules.vtkCommonDataModel import vtkPointData, vtkCellData, vtkUnstructuredGrid


def read_plasma_state(ids_name: str, ids_obj, aos_index_values: dict,
                      subset_idx: int, grid_ggd, ugrid: vtkUnstructuredGrid) -> None:
    """
    Reads plasma state data arrays from the ggd node. These arrays are added as point data or cell data to the
    unstructured grid.
    :param ids_name: name of the top level IDS.
    :param ids_obj: an ids_obj of type ids_name
    :param aos_index_values: the values that shall be used to navigate the AoS and reach the scalar arrays.
    :param subset_idx: an index into grid_ggd/grid_subset AoS
    :param grid_ggd: the grid_ggd (unused)
    :param ugrid: the unstructured grid instance.
    :return: None
    """
    if ids_name == 'distribution_sources':
        # ggd is at /distribution_sources/source(i1)/ggd(itime)
        source_idx = aos_index_values.get('SourceIdx')
        time_idx = aos_index_values.get('TimeIdx')
        try:
            ggd = ids_obj.source[source_idx].ggd[time_idx]
        except IndexError:
            return

        name = 'Particle Density'
        _add_aos_scalar_array_to_vtk_field_data(ggd.particles, subset_idx, name, ugrid)

    elif ids_name == 'distributions':
        # ggd is at /distributions/distribution(i1)/ggd(itime)
        distribution_idx = aos_index_values.get('DistributionIdx')
        time_idx = aos_index_values.get('TimeIdx')
        try:
            ggd = ids_obj.distribution[distribution_idx].ggd[time_idx]
        except IndexError:
            return

        num_expansions = len(ggd.expansion)
        expansion_components = []

        for i in range(num_expansions):
            expansion = ggd.expansion[i].grid_subset
            if subset_idx >= len(expansion):
                continue
            expansion_components.append(expansion[subset_idx].values)
            name = f'Distribution Function Order: {i}'
            _add_aos_scalar_array_to_vtk_field_data(expansion, subset_idx, name, ugrid)

        if len(expansion_components):
            total_df = np.sum(expansion_components, axis=0)
            name = 'Distribution Function'
            _add_scalar_array_to_vtk_field_data(total_df, name, ugrid)

    elif ids_name == 'edge_profiles':
        # ggd is at /edge_profiles/ggd(itime)
        time_idx = aos_index_values.get('TimeIdx')
        try:
            ggd = ids_obj.ggd[time_idx]
        except IndexError:
            return

        # electrons
        #  - temperature
        name = 'Electron Temperature (eV)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.electrons.temperature, subset_idx, name, ugrid)
        #  - density
        name = 'Electron Density (m^-3)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.electrons.density, subset_idx, name, ugrid)
        #  - density_fast
        name = 'Electron Density Fast (m^-3)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.electrons.density_fast, subset_idx, name, ugrid)
        #  - pressure
        name = 'Electron Pressure (Pa)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.electrons.pressure, subset_idx, name, ugrid)
        #  - pressure_fast_perpendicular
        name = 'Electron Pressure Fast Perpendicular (Pa)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.electrons.pressure_fast_perpendicular, subset_idx, name, ugrid)
        #  - pressure_fast_parallel
        name = 'Electron Pressure Fast Parallel (Pa)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.electrons.pressure_fast_parallel, subset_idx, name, ugrid)
        #  - velocity
        name = 'Electron Velocity (m.s^-1)'
        _add_aos_vector_array_to_vtk_field_data(ggd.electrons.velocity, subset_idx, name, ugrid)
        #  - distribution_function
        name = 'Electron Distribution Function'
        _add_aos_scalar_array_to_vtk_field_data(ggd.electrons.distribution_function, subset_idx, name, ugrid)

        # TODO: the remaining arrays. ion(). neutral(), ..

    elif ids_name == 'edge_sources':
        # ggd is at /edge_sources/ggd(itime)
        time_idx = aos_index_values.get('TimeIdx')
        try:
            ggd = ids_obj.ggd[time_idx]
        except IndexError:
            return

        # electrons
        #  - particles
        name = 'Electron Particle Density (m^-3.s^-1)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.electrons.particles, subset_idx, name, ugrid)
        #  - density
        name = 'Electron Energy (W.m^-3)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.electrons.energy, subset_idx, name, ugrid)

        # TODO: the remaining arrays. ion(). neutral(), ..

    elif ids_name == 'edge_transport':
        # ggd is at /edge_transport/model(i1)/ggd(itime)
        model_idx = aos_index_values.get('ModelIdx')
        time_idx = aos_index_values.get('TimeIdx')
        try:
            ggd = ids_obj.model[model_idx].model.ggd[time_idx]
        except IndexError:
            return

        # TODO: all arrays

    elif ids_name == 'equilibrium':
        # ggd is at /equilibrium/time_slice(itime)/ggd(i1)
        time_idx = aos_index_values.get('TimeIdx')
        ggd_idx = aos_index_values.get('GridGGDIdx')
        try:
            ggd = ids_obj.time_slice[time_idx].ggd[ggd_idx]
        except IndexError:
            return

        # - r
        name = 'Major Radius (m)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.r, subset_idx, name, ugrid)
        # - z
        name = 'Height (m)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.z, subset_idx, name, ugrid)
        # - psi
        name = 'Poloidal Flux (Wb)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.psi, subset_idx, name, ugrid)
        # - phi
        name = 'Toroidal Flux (Wb)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.phi, subset_idx, name, ugrid)
        # - theta
        name = 'Poloidal Angle (rad)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.theta, subset_idx, name, ugrid)
        # - j_tor
        name = 'Toroidal Plasma Current Density (A.m^-2)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.j_tor, subset_idx, name, ugrid)
        # - j_parallel
        name = 'Parallel Plasma Current Density (A.m^-2)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.j_parallel, subset_idx, name, ugrid)
        # - b_field_r
        name = 'Magnetic Field Br (T)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.b_field_r, subset_idx, name, ugrid)
        # - b_field_z
        name = 'Magnetic Field Bz (T)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.b_field_z, subset_idx, name, ugrid)
        # - b_field_tor
        name = 'Magnetic Field Btor (T)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.b_field_tor, subset_idx, name, ugrid)

    elif ids_name == 'mhd':
        # TODO: all arrays
        pass
    elif ids_name == 'radiation':
        # TODO: all arrays
        pass
    elif ids_name == 'tf':
        # TODO: all arrays
        pass
    elif ids_name == 'transport_solver_numerics':
        # TODO: all arrays
        pass

    elif ids_name == 'wall':
        # ggd is at /wall/description_ggd(i1)/ggd(itime)
        description_ggd_idx = aos_index_values.get('DescriptionGGDIdx')
        time_idx = aos_index_values.get('TimeIdx')
        try:
            ggd = ids_obj.description_ggd[description_ggd_idx].ggd[time_idx]
        except IndexError:
            return

        # - power_density
        name = 'Power Density (W.m^-2)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.power_density, subset_idx, name, ugrid)
        # - temperature
        name = 'Temperature  (K)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.temperature, subset_idx, name, ugrid)

    elif ids_name == 'waves':
        # TODO: all arrays
        pass


def _add_scalar_array_to_vtk_field_data(array: np.ndarray, name: str, ugrid: vtkUnstructuredGrid) -> None:
    """
    Add a named array as scalars to a vtkUnstructuredGrid instance
    :param array: the numpy array.
    :param name: the name string.
    :param ugrid: an instance of vtkUnstructuredGrid
    :return: None
    """
    point_data: vtkPointData = ugrid.GetPointData()
    num_points = ugrid.GetNumberOfPoints()
    cell_data: vtkCellData = ugrid.GetCellData()
    num_cells = ugrid.GetNumberOfCells()

    vtk_arr = dsa.numpyTovtkDataArray(array, name)
    if len(array) == num_points:
        point_data.AddArray(vtk_arr)
    if len(array) == num_cells:
        cell_data.AddArray(vtk_arr)


def _add_aos_scalar_array_to_vtk_field_data(aos_scalar_node, subset_idx: int, name: str, ugrid: vtkUnstructuredGrid):
    """
    Add the array under the aos_scalar_node to the unstructured grid.
    :param aos_scalar_node: A node with scalar values for each grid subset.
    :param subset_idx: an index into aos_scalar_node
    :param name: this becomes the array name in VTK
    :param ugrid: an unstructured grid instance
    :return: None
    """
    if subset_idx >= len(aos_scalar_node):
        return
    if hasattr(aos_scalar_node[subset_idx], 'values') and len(aos_scalar_node[subset_idx].values):
        _add_scalar_array_to_vtk_field_data(aos_scalar_node[subset_idx].values, name, ugrid)


def _add_aos_vector_array_to_vtk_field_data(aos_vector_node, subset_idx: int, name: str, ugrid: vtkUnstructuredGrid):
    """
    Add the array under the aos_vector_node to the unstructured grid.
    :param aos_vector_node: A node with component vectors for each grid subset.
    :param subset_idx: an index into aos_scalar_node
    :param name: this becomes the array name in VTK
    :param ugrid: an unstructured grid instance
    :return: None
    """
    if subset_idx >= len(aos_vector_node):
        return

    point_data: vtkPointData = ugrid.GetPointData()
    num_points = ugrid.GetNumberOfPoints()
    cell_data: vtkCellData = ugrid.GetCellData()
    num_cells = ugrid.GetNumberOfCells()

    vtk_arr = vtkDoubleArray()
    vtk_arr.SetName(name)
    vtk_arr.SetNumberOfComponents(5)
    num_tuples = 0

    for i, component_name in enumerate(['radial', 'diamagnetic', 'parallel', 'poloidal', 'toroidal']):
        values = getattr(aos_vector_node[subset_idx], component_name)
        vtk_arr.SetComponentName(i, component_name.capitalize())

        if len(values):
            scalar_arr = dsa.numpyTovtkDataArray(values, name + '-' + component_name.capitalize())

            if num_tuples == 0:
                num_tuples = scalar_arr.GetNumberOfTuples()
                vtk_arr.SetNumberOfTuples(num_tuples)

            vtk_arr.CopyComponent(i, scalar_arr, 0)

    if num_tuples == num_points:
        point_data.AddArray(vtk_arr)
    if num_tuples == num_cells:
        cell_data.AddArray(vtk_arr)
