import numpy as np
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonCore import vtkDoubleArray
from vtkmodules.vtkCommonDataModel import vtkPointData, vtkCellData, vtkUnstructuredGrid
import paraview


def read_plasma_state(ids_name: str, ids_obj, aos_index_values: dict,
                      subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:
    """
    Reads plasma state data arrays from the ggd node. These arrays are added as point data or cell data to the
    unstructured grid.
    This is just a dispatching routine, the real work is done in the called
    funtions read_<IDS_name>(ids_obj, aos_index_values, subset_idx, ugrid)
    :param ids_name: name of the top level IDS.
    :param ids_obj: an ids_obj of type ids_name
    :param aos_index_values: the values that shall be used to navigate the AoS and reach the scalar arrays.
    :param subset_idx: an index into grid_ggd/grid_subset AoS
    :param ugrid: the unstructured grid instance.
    :return: None
    """
    if ids_name == 'distribution_sources':
        read_distribution_sources(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == 'distributions':
        read_distributions(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == 'edge_profiles':
        read_edge_profiles(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == 'edge_sources':
        read_edge_sources(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == 'edge_transport':
        read_edge_transport(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == 'equilibrium':
        read_equilibrium(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == 'mhd':
        read_mhd(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == 'radiation':
        read_radiation(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == 'tf':
        read_tf(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == 'transport_solver_numerics':
        read_transport_solver_numerics(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == 'wall':
        read_wall(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == 'waves':
        read_waves(ids_obj, aos_index_values, subset_idx, ugrid)

    else:
        paraview.logger.info(f"Reading plasma state from IDS {ids_name} not implemented.")
        

def read_distribution_sources(ids_obj, aos_index_values: dict,
                              subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:
    
    # ggd is at /distribution_sources/source(i1)/ggd(itime)
    source_idx = aos_index_values.get('SourceIdx')
    time_idx = aos_index_values.get('TimeIdx')
    try:
        ggd = ids_obj.source[source_idx].ggd[time_idx]
    except IndexError:
        return

    name = 'Particle Density'
    _add_aos_scalar_array_to_vtk_field_data(ggd.particles, subset_idx, name, ugrid)


def read_distributions(ids_obj, aos_index_values: dict,
                       subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:
    
    # ggd is at /distributions/distribution(i1)/ggd(itime)
    distribution_idx = aos_index_values.get('DistributionIdx')
    time_idx = aos_index_values.get('TimeIdx')
    try:
        ggd = ids_obj.distribution[distribution_idx].ggd[time_idx]
    except IndexError:
        return

    # - expansion(i2)
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


def read_edge_profiles(ids_obj, aos_index_values: dict,
                       subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:
    
    # ggd is at /edge_profiles/ggd(itime)
    time_idx = aos_index_values.get('TimeIdx')
    try:
        ggd = ids_obj.ggd[time_idx]
    except IndexError:
        return

    # electrons
    e = ggd.electrons # shortcut
    e_description = {
        'Electron Temperature (eV)': e.temperature,
        'Electron Density (m^-3)': e.density,
        'Electron Density Fast (m^-3)': e.density_fast,
        'Electron Pressure (Pa)': e.pressure,
        'Electron Pressure Fast Perpendicular (Pa)': e.pressure_fast_perpendicular,
        'Electron Pressure Fast Parallel (Pa)': e.pressure_fast_parallel,
        'Electron Distribution Function': e.distribution_function
    }
    _multi_add_aos_scalar_to_vtk(e_description, subset_idx, ugrid)
    name = 'Electron Velocity (m.s^-1)'
    _add_aos_vector_array_to_vtk_field_data(e.velocity, subset_idx, name, ugrid)

    # ions
    for i in ggd.ion:
        i_description = {}
        i_description[i.label+' Temperature (eV)'] = i.temperature 
        i_description[i.label+' Density (m^-3)'] = i.density 
        i_description[i.label+' Density Fast (m^-3)'] = i.density_fast 
        i_description[i.label+' Pressure (Pa)'] = i.pressure 
        i_description[i.label+' Pressure Fast Perpendicular (Pa)'] = i.pressure_fast_perpendicular 
        i_description[i.label+' Pressure Fast Parallel (Pa)'] = i.pressure_fast_parallel 
        i_description[i.label+' Kinetic Energy Density'] = i.energy_density_kinetic 
        _multi_add_aos_scalar_to_vtk(i_description, subset_idx, ugrid)
        name = i.label+' Velocity (m.s^-1)'
        _add_aos_vector_array_to_vtk_field_data(i.velocity, subset_idx, name, ugrid)

    # neutrals
    for n in ggd.neutral:
        n_description = {}
        n_description[n.label+' Temperature (eV)'] = n.temperature 
        n_description[n.label+' Density (m^-3)'] = n.density 
        n_description[n.label+' Density Fast (m^-3)'] = n.density_fast 
        n_description[n.label+' Pressure (Pa)'] = n.pressure 
        n_description[n.label+' Pressure Fast Perpendicular (Pa)'] = n.pressure_fast_perpendicular 
        n_description[n.label+' Pressure Fast Parallel (Pa)'] = n.pressure_fast_parallel 
        n_description[n.label+' Kinetic Energy Density'] = n.energy_density_kinetic 
        _multi_add_aos_scalar_to_vtk(n_description, subset_idx, ugrid)
        name = n.label+' Velocity (m.s^-1)'
        _add_aos_vector_array_to_vtk_field_data(n.velocity, subset_idx, name, ugrid)

    # TODO: the remaining arrays. ion(). neutral(), ..


def read_edge_sources(ids_obj, aos_index_values: dict,
                      subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:
    
    # ggd is at /edge_sources/ggd(itime)
    time_idx = aos_index_values.get('TimeIdx')
    try:
        ggd = ids_obj.ggd[time_idx]
    except IndexError:
        return

    # electrons
    #  - particles(i1)
    name = 'Electron Particle Density (m^-3.s^-1)'
    _add_aos_scalar_array_to_vtk_field_data(ggd.electrons.particles, subset_idx, name, ugrid)
    #  - density(i1)
    name = 'Electron Energy (W.m^-3)'
    _add_aos_scalar_array_to_vtk_field_data(ggd.electrons.energy, subset_idx, name, ugrid)

    # TODO: the remaining arrays. ion(). neutral(), ..

def read_edge_transport(ids_obj, aos_index_values: dict,
                        subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:

    # ggd is at /edge_transport/model(i1)/ggd(itime)
    model_idx = aos_index_values.get('ModelIdx')
    time_idx = aos_index_values.get('TimeIdx')
    try:
        ggd = ids_obj.model[model_idx].model.ggd[time_idx]
    except IndexError:
        return

    # TODO: all arrays
    
def read_equilibrium(ids_obj, aos_index_values: dict,
                     subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:

    # ggd is at /equilibrium/time_slice(itime)/ggd(i1)
    time_idx = aos_index_values.get('TimeIdx')
    ggd_idx = aos_index_values.get('GridIdx')
    try:
        ggd = ids_obj.time_slice[time_idx].ggd[ggd_idx]
    except IndexError:
        return

    # - r(i1)
    name = 'Major Radius (m)'
    _add_aos_scalar_array_to_vtk_field_data(ggd.r, subset_idx, name, ugrid)
    # - z(i1)
    name = 'Height (m)'
    _add_aos_scalar_array_to_vtk_field_data(ggd.z, subset_idx, name, ugrid)
    # - psi(i1)
    name = 'Poloidal Flux (Wb)'
    _add_aos_scalar_array_to_vtk_field_data(ggd.psi, subset_idx, name, ugrid)
    # - phi(i1)
    name = 'Toroidal Flux (Wb)'
    _add_aos_scalar_array_to_vtk_field_data(ggd.phi, subset_idx, name, ugrid)
    # - theta(i1)
    name = 'Poloidal Angle (rad)'
    _add_aos_scalar_array_to_vtk_field_data(ggd.theta, subset_idx, name, ugrid)
    # - j_tor(i1)
    name = 'Toroidal Plasma Current Density (A.m^-2)'
    _add_aos_scalar_array_to_vtk_field_data(ggd.j_tor, subset_idx, name, ugrid)
    # - j_parallel(i1)
    name = 'Parallel Plasma Current Density (A.m^-2)'
    _add_aos_scalar_array_to_vtk_field_data(ggd.j_parallel, subset_idx, name, ugrid)
    # - b_field_r(i1)
    name = 'Magnetic Field Br (T)'
    _add_aos_scalar_array_to_vtk_field_data(ggd.b_field_r, subset_idx, name, ugrid)
    # - b_field_z(i1)
    name = 'Magnetic Field Bz (T)'
    _add_aos_scalar_array_to_vtk_field_data(ggd.b_field_z, subset_idx, name, ugrid)
    # - b_field_tor(i1)
    name = 'Magnetic Field Btor (T)'
    _add_aos_scalar_array_to_vtk_field_data(ggd.b_field_tor, subset_idx, name, ugrid)


def read_mhd(ids_obj, aos_index_values: dict,
             subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:

    # ggd is at /mhd/ggd(itime)
    time_idx = aos_index_values.get('TimeIdx')
    try:
        ggd = ids_obj.ggd[time_idx]
    except IndexError:
        return

    # TODO: all arrays


def read_radiation(ids_obj, aos_index_values: dict,
                   subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:
    
    # ggd is at /radiation/process(i1)/ggd(itime)
    # TODO: ProcessIdx is not detected by generator.py and hence it is not present.
    process_idx = aos_index_values.get('ProcessIdx') or 0
    time_idx = aos_index_values.get('TimeIdx')
    try:
        ggd = ids_obj.process[process_idx].ggd[time_idx]
    except IndexError:
        return

    # TODO: all arrays


def read_tf(ids_obj, aos_index_values: dict,
            subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:
    # A 'ggd' node is absent,
    # but there are scalar arrays.
    # /tf/field_map(itime)/
    #  - b_field_r(i1)
    #  - b_field_z(i1)
    #  - b_field_tor(i1)
    #  - a_field_r(i1)
    #  - a_field_z(i1)
    #  - a_field_tor(i1)
    time_idx = aos_index_values.get('TimeIdx')
    for attr_name in ['b_field_r', 'b_field_z', 'b_field_tor', 'a_field_r', 'a_field_z', 'a_field_tor']:
        try:
            aos_scalar_node = getattr(ids_obj.field_map[time_idx], attr_name)
        except IndexError:
            continue
        component_name = attr_name.split('_')[-1]
        name = f'Magnetic Field B{component_name} (T)'
        _add_aos_scalar_array_to_vtk_field_data(aos_scalar_node, subset_idx, name, ugrid)


def read_transport_solver_numerics(ids_obj, aos_index_values: dict,
                                   subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:
    
    # ggd is at /transport_solver_numerics/boundary_conditions_ggd(itime)
    time_idx = aos_index_values.get('TimeIdx')
    try:
        ggd = ids_obj.boundary_conditions_ggd[time_idx]
    except IndexError:
        return

    # - current(i1)
    name = f'Current Boundary Condition - {ggd.current.identifier.name.capitalize()}'
    _add_aos_scalar_array_to_vtk_field_data(ggd.current, subset_idx, name, ugrid)
    # - electrons/particles(i1)
    try:
        name = f'Electron Density Boundary Condition - {ggd.electrons.particles[subset_idx].identifier.name} (m^-3.s^-1)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.electrons.particles, subset_idx, name, ugrid)
    except IndexError:
        # can raise from name = "...", when subset_idx >= ggd.electrons.particles
        pass
    # - electrons/energy(i1)
    try:
        name = f'Electron Energy Boundary Condition - {ggd.electrons.energy[subset_idx].identifier.name} (W.m^-3)'
        _add_aos_scalar_array_to_vtk_field_data(ggd.electrons.energy, subset_idx, name, ugrid)
    except IndexError:
        # can raise from name = "...", when subset_idx >= ggd.electrons.energy
        pass

    # TODO: ion arrays


def read_wall(ids_obj, aos_index_values: dict,
              subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:
    
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


def read_waves(ids_obj, aos_index_values: dict,
               subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:

    # A 'ggd' node is absent,
    # but there are scalar and vector arrays.
    # /waves/coherent_wave(i1)/full_wave(itime)/
    #  - e_field/plus(i2)
    #  - e_field/minus(i2)
    #  - e_field/parallel(i2)
    #  - e_field/normal(i2)
    #  - e_field/bi_normal(i2)
    #  - b_field/plus(i2)
    #  - b_field/minus(i2)
    #  - b_field/parallel(i2)
    #  - b_field/normal(i2)
    #  - b_field/bi_normal(i2)
    #  - k_perpendicular(i2)
    coherent_wave_idx = aos_index_values.get('CoherentWaveIdx')
    time_idx = aos_index_values.get('TimeIdx')

    try:
        e_field = ids_obj.coherent_wave[coherent_wave_idx].full_wave[time_idx].e_field
        name = 'Electric Field - LH polarized (V.m^-1)'
        _add_aos_scalar_array_to_vtk_field_data(e_field.plus, subset_idx, name, ugrid)
        name = 'Electric Field - RH polarized (V.m^-1)'
        _add_aos_scalar_array_to_vtk_field_data(e_field.minus, subset_idx, name, ugrid)
        name = 'Electric Field - Parallel to B (V.m^-1)'
        _add_aos_scalar_array_to_vtk_field_data(e_field.parallel, subset_idx, name, ugrid)
        name = 'Electric Field - Normal to flux surface (V.m^-1)'
        _add_aos_scalar_array_to_vtk_field_data(e_field.normal, subset_idx, name, ugrid)
        name = 'Electric Field - Tangential to flux surface (V.m^-1)'
        _add_aos_scalar_array_to_vtk_field_data(e_field.bi_normal, subset_idx, name, ugrid)
    except IndexError:
        pass

    try:
        # WARNING: Are these units correct?
        b_field = ids_obj.coherent_wave[coherent_wave_idx].full_wave[time_idx].b_field
        name = 'Magnetic Field - Parallel to B (V.m^-1)'
        _add_aos_scalar_array_to_vtk_field_data(b_field.parallel, subset_idx, name, ugrid)
        name = 'Magnetic Field - Normal to flux surface (V.m^-1)'
        _add_aos_scalar_array_to_vtk_field_data(b_field.normal, subset_idx, name, ugrid)
        name = 'Magnetic Field - Tangential to flux surface (V.m^-1)'
        _add_aos_scalar_array_to_vtk_field_data(b_field.bi_normal, subset_idx, name, ugrid)
    except IndexError:
        pass

    try:
        # WARNING: Are these units correct?
        name = 'Perpendicular Wave Vector (V.m^-1)'
        _add_aos_scalar_array_to_vtk_field_data(
            ids_obj.coherent_wave[coherent_wave_idx].full_wave[time_idx].k_perpendicular, subset_idx, name, ugrid)
    except IndexError:
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
    print(vtk_arr, array)
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


def _multi_add_aos_scalar_to_vtk(scalar_description: dict, subset_idx: int, ugrid: vtkUnstructuredGrid):
    """
    Add the dict of names and arrays in scalar_description to the unstructured grid.
    :param scalar_description: A dictionary of {'name, aos_scalar_node} 
    Each aos_scalar_node is a node with an array of scalar values for each grid subset.
    Each name becomes the respective array name in VTK.
    :param subset_idx: an index into aos_scalar_node
    :param ugrid: an unstructured grid instance
    :return: None
    """
    for name, node in scalar_description.items():
        _add_aos_scalar_array_to_vtk_field_data(node, subset_idx, name, ugrid)
    

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

