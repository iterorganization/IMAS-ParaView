import numpy as np
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonCore import vtkDoubleArray
from vtkmodules.vtkCommonDataModel import vtkPointData, vtkCellData, vtkUnstructuredGrid
import paraview

from imas import dd_units
dd_units = dd_units.DataDictionaryUnits()

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
        paraview.logger.info(f'Reading plasma state from IDS {ids_name} not implemented.')
        

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

    # for finding units:
    ids_name = 'edge_profiles'
    ggd_path = 'ggd'
    
    # electrons
    units_path = ggd_path + '/electrons'
    # scalar arrays:
    quantities = {
        'temperature': 'Electron Temperature',
        'density': 'Electron Density',
        'density_fast': 'Electron Density (fast)',
        'pressure': 'Electron Pressure',
        'pressure_fast_perpendicular': 'Electron Pressure (fast perpendicular)',
        'pressure_fast_parallel': 'Electron Pressure (fast parallel)',
        'distribution_function': 'Electron Distribution Function'
    }
    for q_name in quantities:
        units = dd_units.get_units(ids_name, f'{units_path}/{q_name}')
        _add_aos_scalar_array_to_vtk_field_data(getattr(ggd.electrons, q_name),
                                                subset_idx,
                                                f'{quantities[q_name]} [${units}$]',
                                                ugrid)
    # vector array:
    units = dd_units.get_units(ids_name, units_path+'/velocity')
    name = f'Electron Velocity [${units}$]'
    _add_aos_vector_array_to_vtk_field_data(ggd.electrons.velocity, subset_idx, name, ugrid)

    # ions, neutrals = heavy particles = hp
    # Add ' (all ions)' to the ion name to distinguish from state below:
    ion_names = [' (all ions)']*len(ggd.ion) + ['']*len(ggd.neutral)
    for hp, i_name in list(zip([*ggd.ion, *ggd.neutral], ion_names)):
        hp_name = hp.label + i_name # ' add (all ions)' to ions or '' to neutrals
        if i_name == '': units_hp_path = ggd_path+'/neutral'
        else: units_hp_path = ggd_path+'/ion'
        # scalar arrays:
        quantities = {
            'temperature': f'{hp_name} Temperature',
            'density': f'{hp_name} Density',
            'density_fast': f'{hp_name} Density (fast)',
            'pressure': f'{hp_name} Pressure',
            'pressure_fast_perpendicular': f'{hp_name} Pressure (fast perpendicular)',
            'pressure_fast_parallel': f'{hp_name} Pressure (fast parallel)',
            'energy_density_kinetic': f'{hp_name} Kinetic Energy Density'
        }
        for q_name in quantities:
            units = dd_units.get_units(ids_name, f'{units_hp_path}/{q_name}')
            _add_aos_scalar_array_to_vtk_field_data(getattr(hp, q_name),
                                                    subset_idx,
                                                    f'{quantities[q_name]} [${units}$]',
                                                    ugrid)            
        # vector array:
        units = dd_units.get_units(ids_name, units_hp_path+'/velocity')
        name = f'{hp_name} Velocity [${units}$]'
        _add_aos_vector_array_to_vtk_field_data(hp.velocity, subset_idx, name, ugrid)
        # heavy particles: state
        units_hp_path += '/state'
        for state in hp.state:
            # scalars:
            quantities = {
                'temperature': f'{state.label} Temperature',
                'density': f'{state.label} Density',
                'density_fast': f'{state.label} Density (fast)',
                'pressure': f'{state.label} Pressure',
                'pressure_fast_perpendicular': f'{state.label} Pressure (fast perpendicular)',
                'pressure_fast_parallel': f'{state.label} Pressure (fast parallel)',
                'energy_density_kinetic': f'{state.label} Kinetic Energy Density',
                'distribution_function': f'{state.label} Distribution Function'
            }
            if i_name: # only for ions
                quantities['z_average'] = f'{state.label} $<Z>$'
                quantities['z_square_average'] = f'{state.label} $<Z^2>$'
                quantities['ionisation_potential'] = f'{state.label} Ionisation Potential'
            for q_name in quantities:
                units = dd_units.get_units(ids_name, f'{units_hp_path}/{q_name}')
                _add_aos_scalar_array_to_vtk_field_data(getattr(state, q_name),
                                                        subset_idx,
                                                        f'{quantities[q_name]} [${units}$]',
                                                        ugrid)
            # vectors:
            quantities = {
                'velocity': f'{state.label} Velocity',
                'velocity_diamagnetic': f'{state.label} Velocity (diamagnetic)',
                'velocity_exb': f'{state.label} Velocity (ExB)'
            }
            for q_name in quantities:
                units = dd_units.get_units(ids_name, f'{units_hp_path}/{q_name}')
                _add_aos_vector_array_to_vtk_field_data(getattr(state, q_name),
                                                        subset_idx,
                                                        f'{quantities[q_name]} [${units}$]',
                                                        ugrid)

    # Other quantities;
    # The tuple has ('text description', scalar?, units?)
    units_path = ggd_path
    quantities = {
        't_i_average': 'Average t_i',
        'n_i_total_over_n_e': 'Total n_i over n_e',
        'zeff': 'Zeff',
        'pressure_thermal': 'Pressure (thermal)',
        'pressure_perpendicular': 'Pressure (perpendicular)',
        'pressure_parallel': 'Pressure (parallel)',
        'j_parallel': 'Current (parallel)',
        'phi_potential': 'Potential Phi'
    }
    for q_name in quantities:
        units = dd_units.get_units(ids_name, f'{units_path}/{q_name}')
        _add_aos_scalar_array_to_vtk_field_data(getattr(ggd, q_name),
                                                subset_idx,
                                                f'{quantities[q_name]} [${units}$]',
                                                ugrid)

    # other vector quantities
    quantities = {
        'j_total': 'Current (total)',
        'j_anomalous': 'Current (anomalous)',
        'j_inertial': 'Current (inertial)',
        'j_ion_neutral_friction': 'Current (ion+neutral friction)',
        'j_parallel_viscosity': 'Current (parallel viscosity)',
        'j_perpendicular_viscosity': 'Current (perpendicular viscosity)',
        'j_heat_viscosity': 'Current (heat viscosity)',
        'j_pfirsch_schlueter': 'Current (Pfirsch-Schlueter effects)',
        'j_diamagnetic': 'Current (diamagnetic)',
        'e_field': 'E'
    }
    for q_name in quantities:
        units = dd_units.get_units(ids_name, f'{units_path}/{q_name}')
        _add_aos_vector_array_to_vtk_field_data(getattr(ggd, q_name),
                                                subset_idx,
                                                f'{quantities[q_name]} [${units}$]',
                                                ugrid)
        
    # TODO: ggd_fast...

def read_edge_profiles_lecad(ids_obj, aos_index_values: dict,
                       subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:
        # ggd is at /edge_profiles/ggd(itime)
        time_idx = aos_index_values.get('TimeIdx')
        try:
            ggd = ids_obj.ggd[time_idx]
        except IndexError:
            return

        # quantities for electrons, ions and neutrals
        # tuple shape: (object attribute, if values are scalars or not)
        # velocity values are vectors, other are scalars
        quantityDict = {
            'Density [$m^{-3}$]' : ('density', True),
            'Temperature [eV]' : ('temperature', True),
            'Density Fast [$m^{-3}$]' : ('density_fast', True),
            'Pressure [Pa]' : ('pressure', True),
            'Pressure Fast Perpendicular [Pa]' : ('pressure_fast_perpendicular', True),
            'Pressure Fast Parallel [Pa]' : ('pressure_fast_parallel', True),
            'Velocity [$m/s$]' : ('velocity', False),
        }

        for i in quantityDict:
            (path, scalar) = quantityDict[i]
            # electrons
            if scalar:
                _add_aos_scalar_array_to_vtk_field_data(getattr(ggd.electrons, path), 
                    subset_idx, 'Electron ' + i, ugrid)
            else:
                _add_aos_vector_array_to_vtk_field_data(getattr(ggd.electrons, path), 
                    subset_idx, 'Electron ' + i, ugrid)

            # ions
            for j in range(len(ggd.ion)):
                # find ion label
                if ggd.ion[j].label != '':
                    ionLabel = ggd.ion[j].label
                else:
                    ionLabel = ggd.ion[j].state[0].label
                for k in range(len(getattr(ggd.ion[j], path))):
                    # Check if there is actually an array
                    if len(getattr(ggd.ion[j], path)[k].values) > 0:
                        if scalar:
                            _add_aos_scalar_array_to_vtk_field_data(getattr(ggd.ion[j], path), 
                                subset_idx, 'Ion ' + ionLabel + ' ' + i, ugrid)
                        else:
                            _add_aos_vector_array_to_vtk_field_data(getattr(ggd.ion[j], path), 
                                subset_idx, 'Ion ' + ionLabel + ' ' + i, ugrid)
                    else:
                        break
            # neutrals
            for j in range(len(ggd.neutral)):
                # find neutral label
                if ggd.neutral[j].label != '':
                    neutralLabel = ggd.neutral[j].label
                else:
                    neutralLabel = ggd.neutral[j].state[0].label
                try:
                    for k in range(len(getattr(ggd.neutral[j], path))):
                        # Check if there is actually an array
                        if len(getattr(ggd.neutral[j], path)[k].values) > 0:
                            if scalar:
                                _add_aos_scalar_array_to_vtk_field_data(getattr(ggd.neutral[j], path), 
                                    subset_idx, 'Neutral ' + neutralLabel + ' ' + i, ugrid)
                            else:
                                _add_aos_vector_array_to_vtk_field_data(getattr(ggd.neutral[j], path), 
                                    subset_idx, 'Neutral ' + neutralLabel + ' ' + i, ugrid)
                        else:
                            break
                except:
                    break
                    
        # TODO ggd_fast



def read_edge_sources(ids_obj, aos_index_values: dict,
                      subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:
    
    # ggd is at /edge_sources/source/ggd(itime)
    time_idx = aos_index_values.get('TimeIdx')
    model_idx = aos_index_values.get('ModelIdx')
    # For units reading:
    ids_name = 'edge_sources'
    ggd_path = 'source/ggd'
    
    for source in ids_obj.source:
        try:
            ggd = source.ggd[time_idx]
        except IndexError:
            return
        s_name = source.identifier.name
        
        # electrons
        units_path = ggd_path + '/electrons'
        units = dd_units.get_units(ids_name, units_path+'/particles')
        name = f'Electron Density ({s_name}) [${units}$]'
        _add_aos_scalar_array_to_vtk_field_data(ggd.electrons.particles, subset_idx, name, ugrid)
        units = dd_units.get_units(ids_name, units_path+'/energy')
        name = f'Electron Energy ({s_name}) [${units}$]'
        _add_aos_scalar_array_to_vtk_field_data(ggd.electrons.energy, subset_idx, name, ugrid)

        # ions, neutrals = heavy particles = hp
        # Add ' (all ions)' to the ion name to distinguish from state below:
        ion_names = [' (all ions)']*len(ggd.ion) + ['']*len(ggd.neutral)
        for hp, i_name in list(zip([*ggd.ion, *ggd.neutral], ion_names)):
            hp_name = hp.label + i_name # add ' (all ions)' to ions or '' to neutrals
            if i_name == '': units_hp_path = ggd_path+'/neutral'
            else: units_hp_path = ggd_path+'/ion'
            units = dd_units.get_units(ids_name, units_hp_path+'/particles')
            name = f'{hp_name} Density ({s_name}) [${units}$]'
            _add_aos_scalar_array_to_vtk_field_data(hp.particles, subset_idx, name, ugrid)
            units = dd_units.get_units(ids_name, units_hp_path+'/energy')
            name = f'{hp_name} Energy ({s_name}) [${units}$]'
            _add_aos_scalar_array_to_vtk_field_data(hp.energy, subset_idx, name, ugrid)
            units = dd_units.get_units(ids_name, units_hp_path+'/momentum')
            name = f'{hp_name} Momentum ({s_name}) [${units}$]'
            _add_aos_vector_array_to_vtk_field_data(hp.momentum, subset_idx, name, ugrid)
            # heavy particles: state
            for state in hp.state:
                # sometimes state.label is not filled:
                try: st_name = state.label
                except: st_name = '? ' + hp_name
                units_path = units_hp_path + '/state'
                units = dd_units.get_units(ids_name, units_path+'/particles')
                name = f'{st_name} Density ({s_name}) [${units}$]'
                _add_aos_scalar_array_to_vtk_field_data(state.particles, subset_idx, name, ugrid)
                units = dd_units.get_units(ids_name, units_path+'/energy')
                name = f'{st_name} Energy ({s_name}) [${units}$]'
                _add_aos_scalar_array_to_vtk_field_data(state.energy, subset_idx, name, ugrid)
                units = dd_units.get_units(ids_name, units_path+'/momentum')
                name = f'{st_name} Momentum ({s_name}) [${units}$]'
                _add_aos_vector_array_to_vtk_field_data(state.momentum, subset_idx, name, ugrid)

        # total_ion_energy
        units = dd_units.get_units(ids_name, ggd_path+'/total_ion_energy')
        name = f'Total Ion Energy ({s_name}) [${units}$]'
        _add_aos_scalar_array_to_vtk_field_data(ggd.total_ion_energy, subset_idx, name, ugrid)

        # momentum
        units = dd_units.get_units(ids_name, ggd_path+'/momentum')
        name = f'Momentum ({s_name}) [${units}$]'
        _add_aos_vector_array_to_vtk_field_data(ggd.momentum, subset_idx, name, ugrid)

        # current
        units = dd_units.get_units(ids_name, ggd_path+'/current')
        name = f'Current ({s_name}) [${units}$]'
        _add_aos_scalar_array_to_vtk_field_data(ggd.current, subset_idx, name, ugrid)

        # TODO: ggd_fast

def read_edge_transport(ids_obj, aos_index_values: dict,
                        subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:

    # ggd is at /edge_transport/model(i1)/ggd(itime)
    # TODO: ModelIdx is not detected by generator.py and hence it is not present.
    # model_idx = aos_index_values.get('ModelIdx') or 0
    # time_idx = aos_index_values.get('TimeIdx')
    # try:
    #     ggd = ids_obj.model[model_idx].ggd[time_idx]
    # except IndexError:
    #     return

    time_idx = aos_index_values.get('TimeIdx')
    ids_name = 'edge_transport'
    ggd_path = 'model/ggd'
    for model in ids_obj.model:
        try:
            ggd = model.ggd[time_idx]
        except IndexError:
            continue

        # Use model name but only if we have more than one model:
        if len(ids_obj.model) > 1: m_name = f'({model.identifier.name})'
        else: m_name = ''

        # conductivity
        units_path = ggd_path + '/conductivity'
        units = dd_units.get_units(ids_name, units_path)
        name = f'Conductivity {m_name} [${units}$]' 
        _add_aos_vector_array_to_vtk_field_data(ggd.conductivity, subset_idx, name, ugrid)

        # electrons
        electrons = ggd.electrons
        for q_name in ['particles', 'energy']:
            units_path = ggd_path + '/electrons/' + q_name
            quantity = getattr(electrons, q_name)
            units = dd_units.get_units(ids_name, units_path+'/d')
            name = f'Electron Diffusivity ({q_name}) [${units}$] {m_name}'
            _add_aos_scalar_array_to_vtk_field_data(quantity.d, subset_idx, name, ugrid)
            units = dd_units.get_units(ids_name, units_path+'/v')
            name = f'Electron Convection ({q_name}) [${units}$] {m_name}'
            _add_aos_scalar_array_to_vtk_field_data(quantity.v, subset_idx, name, ugrid)
            units = dd_units.get_units(ids_name, units_path+'/flux')
            name = f'Electron Flux ({q_name}) [{units}] {m_name}' 
            _add_aos_scalar_array_to_vtk_field_data(quantity.flux, subset_idx, name, ugrid)
            name = f'Electron Flux Limiter Coefficient ({q_name}) {m_name}'
            _add_aos_scalar_array_to_vtk_field_data(quantity.flux_limiter, subset_idx, name, ugrid)
    
        # total_ion_energy
        units_path = ggd_path + '/total_ion_energy'
        quantity = ggd.total_ion_energy
        q_name = 'Total Ion Energy'
        units = dd_units.get_units(ids_name, units_path+'/d')
        name = f'{q_name} Diffusivity [${units}$] {m_name}'
        _add_aos_scalar_array_to_vtk_field_data(quantity.d, subset_idx, name, ugrid)
        units = dd_units.get_units(ids_name, units_path+'/v')
        name = f'{q_name} Convection [${units}$] {m_name}'
        _add_aos_scalar_array_to_vtk_field_data(quantity.v, subset_idx, name, ugrid)
        units = dd_units.get_units(ids_name, units_path+'/flux')
        name = f'{q_name} Flux [${units}$] {m_name}'
        _add_aos_scalar_array_to_vtk_field_data(quantity.flux, subset_idx, name, ugrid)
        name = f'{q_name} Flux Limiter Coefficient {m_name}'
        _add_aos_scalar_array_to_vtk_field_data(quantity.flux_limiter, subset_idx, name, ugrid)

        # momentum
        units_path = ggd_path + '/momentum'
        quantity =ggd.momentum
        q_name = 'Momentum'
        units = dd_units.get_units(ids_name, units_path+'/d')
        name = f'{q_name} Diffusivity [${units}$] {m_name}'
        _add_aos_vector_array_to_vtk_field_data(quantity.d, subset_idx, name, ugrid)
        units = dd_units.get_units(ids_name, units_path+'/v')
        name = f'{q_name} Convection [${units}$] {m_name}'
        _add_aos_vector_array_to_vtk_field_data(quantity.v, subset_idx, name, ugrid)
        units = dd_units.get_units(ids_name, units_path+'/flux')
        name = f'{q_name} Flux [${units}$] {m_name}'
        _add_aos_vector_array_to_vtk_field_data(quantity.flux, subset_idx, name, ugrid)
        name = f'{q_name} Flux Limiter Coefficient {m_name}'
        _add_aos_vector_array_to_vtk_field_data(quantity.flux_limiter, subset_idx, name, ugrid)

        # ions, neutrals = heavy particles = hp
        # Add ' (all ions)' to the ion name to distinguish from state below:
        ion_names = [' (all ions)']*len(ggd.ion) + ['']*len(ggd.neutral)
        for hp, i_name in list(zip([*ggd.ion, *ggd.neutral], ion_names)):
            hp_name = hp.label + i_name # add ' (all ions)' to ions or '' to neutrals
            if i_name == '': units_hp_path = ggd_path+'/neutral'
            else: units_hp_path = ggd_path+'/ion'
            # heavy particles: particles, energy
            for q_name in ['particles', 'energy']:
                quantity = getattr(hp, q_name)
                units_path = units_hp_path + '/' + q_name
                units = dd_units.get_units(ids_name, units_path+'/d')
                name = f'{hp_name} Diffusivity ({q_name}) [${units}$] {m_name}'
                _add_aos_scalar_array_to_vtk_field_data(quantity.d, subset_idx, name, ugrid)
                units = dd_units.get_units(ids_name, units_path+'/v')
                name = f'{hp_name} Convection ({q_name}) [${units}$] {m_name}'
                _add_aos_scalar_array_to_vtk_field_data(quantity.v, subset_idx, name, ugrid)
                units = dd_units.get_units(ids_name, units_path+'/flux')
                name = f'{hp_name} Flux ({q_name}) [${units}$] {m_name}'
                _add_aos_scalar_array_to_vtk_field_data(quantity.flux, subset_idx, name, ugrid)
                name = f'{hp_name} Flux Limiter Coefficient ({q_name}) {m_name}'
                _add_aos_scalar_array_to_vtk_field_data(quantity.flux_limiter, subset_idx, name, ugrid)
            # heavy particles: momentum
            quantity =hp.momentum
            q_name = 'Momentum'
            units_path = units_hp_path + '/momentum'
            units = dd_units.get_units(ids_name, units_path+'/d')
            name = f'{hp_name} Diffusivity ({q_name}) [${units}$] {m_name}'
            _add_aos_vector_array_to_vtk_field_data(quantity.d, subset_idx, name, ugrid)
            units = dd_units.get_units(ids_name, units_path+'/v')
            name = f'{hp_name} Convection ({q_name}) [${units}$] {m_name}'
            _add_aos_vector_array_to_vtk_field_data(quantity.v, subset_idx, name, ugrid)
            units = dd_units.get_units(ids_name, units_path+'/flux')
            name = f'{hp_name} Flux ({q_name}) [${units}$] {m_name}'
            _add_aos_vector_array_to_vtk_field_data(quantity.flux, subset_idx, name, ugrid)
            name = f'{hp_name} Flux Limiter Coefficient ({q_name}) {m_name}'
            _add_aos_vector_array_to_vtk_field_data(quantity.flux_limiter, subset_idx, name, ugrid)
            # heavy particles: state
            for state in hp.state:
                for q_name in ['particles', 'energy']:
                    units_path = units_hp_path + '/state/' + q_name
                    quantity = getattr(state, q_name)
                    units = dd_units.get_units(ids_name, units_path+'/d')
                    name = f'{state.label} Diffusivity ({q_name}) [${units}$] {m_name}'
                    _add_aos_scalar_array_to_vtk_field_data(quantity.d, subset_idx, name, ugrid)
                    units = dd_units.get_units(ids_name, units_path+'/v')
                    name = f'{state.label} Convection ({q_name}) [${units}$] {m_name}'
                    _add_aos_scalar_array_to_vtk_field_data(quantity.v, subset_idx, name, ugrid)
                    units = dd_units.get_units(ids_name, units_path+'/flux')
                    name = f'{state.label} Flux ({q_name}) [${units}$] {m_name}'
                    _add_aos_scalar_array_to_vtk_field_data(quantity.flux, subset_idx, name, ugrid)
                    name = f'{state.label} Flux Limiter Coefficient ({q_name}) {m_name}'
                    _add_aos_scalar_array_to_vtk_field_data(quantity.flux_limiter, subset_idx, name, ugrid)
                
                q_name = 'momentum'
                units_path = units_hp_path + '/state/' + q_name
                quantity = state.momentum
                units = dd_units.get_units(ids_name, units_path+'/d')
                name = f'{state.label} Diffusivity ({q_name}) [${units}$] {m_name}'
                _add_aos_vector_array_to_vtk_field_data(quantity.d, subset_idx, name, ugrid)
                units = dd_units.get_units(ids_name, units_path+'/v')
                name = f'{state.label} Convection ({q_name}) [${units}$] {m_name}'
                _add_aos_vector_array_to_vtk_field_data(quantity.v, subset_idx, name, ugrid)
                units = dd_units.get_units(ids_name, units_path+'/flux')
                name = f'{state.label} Flux ({q_name}) [${units}$] {m_name}'
                _add_aos_vector_array_to_vtk_field_data(quantity.flux, subset_idx, name, ugrid)
                name = f'{state.label} Flux Limiter Coefficient ({q_name}) {m_name}'
                _add_aos_vector_array_to_vtk_field_data(quantity.flux_limiter, subset_idx, name, ugrid)
    # TODO: ggd_fast
    
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
    name = 'Major Radius [$m$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.r, subset_idx, name, ugrid)
    # - z(i1)
    name = 'Height [$m$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.z, subset_idx, name, ugrid)
    # - psi(i1)
    name = 'Poloidal Flux [$Wb$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.psi, subset_idx, name, ugrid)
    # - phi(i1)
    name = 'Toroidal Flux [$Wb$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.phi, subset_idx, name, ugrid)
    # - theta(i1)
    name = 'Poloidal Angle [$rad$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.theta, subset_idx, name, ugrid)
    # - j_tor(i1)
    name = 'Toroidal Plasma Current Density [$A m^{-2}$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.j_tor, subset_idx, name, ugrid)
    # - j_parallel(i1)
    name = 'Parallel Plasma Current Density [$A m^{-2}$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.j_parallel, subset_idx, name, ugrid)
    # - b_field_r(i1)
    name = 'Magnetic Field Br [$T$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.b_field_r, subset_idx, name, ugrid)
    # - b_field_z(i1)
    name = 'Magnetic Field Bz [$T$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.b_field_z, subset_idx, name, ugrid)
    # - b_field_tor(i1)
    name = 'Magnetic Field Btor [$T$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.b_field_tor, subset_idx, name, ugrid)


def read_mhd(ids_obj, aos_index_values: dict,
             subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:

    # ggd is at /mhd/ggd(itime)
    time_idx = aos_index_values.get('TimeIdx')
    try:
        ggd = ids_obj.ggd[time_idx]
    except IndexError:
        return

    # - electrons.temperature
    name = 'Electron Temperature [$eV$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.electron.temperature, subset_idx, name, ugrid)
    # - t_i_average
    name = 'Ion Temperature (average) [$eV$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.t_i_average, subset_idx, name, ugrid)
    # - n_i_total
    name = 'Ion Density (total) [$m^{-3}$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.n_i_total, subset_idx, name, ugrid)
    # - zeff
    name = 'Z effective'
    _add_aos_scalar_array_to_vtk_field_data(ggd.zeff, subset_idx, name, ugrid)
    # - b_field_r
    name = 'Magnetic Field Br [$T$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.b_field_r, subset_idx, name, ugrid)
    # - b_field_z
    name = 'Magnetic Field Bz [$T$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.b_field_z, subset_idx, name, ugrid)
    # - b_field_tor
    name = 'Magnetic Field Btor [$T$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.b_field_tor, subset_idx, name, ugrid)
    # - a_field_r
    name = 'Magnetic Potential Ar [$T m$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.a_field_r, subset_idx, name, ugrid)
    # - a_field_z
    name = 'Magnetic Potential Az [$T m$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.a_field_z, subset_idx, name, ugrid)
    # - a_field_tor
    name = 'Magnetic Potential Ator [$T m$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.a_field_tor, subset_idx, name, ugrid)
    # - PSI
    name = 'Poloidal Flux [$Wb$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.psi, subset_idx, name, ugrid)
    # - velocity_r
    name = 'Plasma Velocity Vr [$T m$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.velocity_r, subset_idx, name, ugrid)
    # - velocity_z
    name = 'Plasma Velocity Vz [$T m$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.velocity_z, subset_idx, name, ugrid)
    # - velocity_tor
    name = 'Plasma Velocity Vtor [$T m$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.velocity_tor, subset_idx, name, ugrid)
    # - velocity_parallel
    name = 'Plasma Velocity Vparallel [$T m$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.velocity_parallel, subset_idx, name, ugrid)
    # - phi_potential
    name = 'Electric Potential [$V$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.phi_potential, subset_idx, name, ugrid)
    # - vorticity
    name = 'Vorticity [$s^{-1}$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.vorticity, subset_idx, name, ugrid)
    # - J_r
    name = 'Current Density Jr [$A m^{-2}$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.j_r, subset_idx, name, ugrid)
    # - J_z
    name = 'Current Density Jz [$A m^{-2}$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.j_z, subset_idx, name, ugrid)
    # - J_phi
    name = 'Current Density Jtor [$A m^{-2}$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.j_tor, subset_idx, name, ugrid)
    # - mass_density
    name = 'Mass Density [$kg m^{-3}$]'
    _add_aos_scalar_array_to_vtk_field_data(ggd.mass_density, subset_idx, name, ugrid)

    # NOT TESTED! Can't find a data entry to test.

def read_radiation(ids_obj, aos_index_values: dict,
                   subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:
    
    # ggd is at /radiation/process(i1)/ggd(itime)
    # TODO: ProcessIdx is not detected by generator.py and hence it is not present.
    # process_idx = aos_index_values.get('ProcessIdx') or 0
    # time_idx = aos_index_values.get('TimeIdx')
    # try:
    #     ggd = ids_obj.process[process_idx].ggd[time_idx]
    # except IndexError:
    #     return

    time_idx = aos_index_values.get('TimeIdx')
    for process in ids_obj.process:
        try:
            ggd = process.ggd[time_idx]
        except IndexError:
            continue
        p_name = process.identifier.name

        # electrons.emissivity
        quantity = ggd.electrons.emissivity
        name = f'Electron Emissivity ({p_name}) ' + '[$W m^{-3}$]'
        _add_aos_scalar_array_to_vtk_field_data(quantity, subset_idx, name, ugrid)
        
        # ions, neutrals = heavy particles = hp
        for hp in [*ggd.ion, *ggd.neutral]:

            # heavy particles: emissivity
            quantity = hp.emissivity
            name = f'{hp.label} Emissivity ({p_name}) ' + '[$W m^{-3}$]'
            _add_aos_scalar_array_to_vtk_field_data(quantity, subset_idx, name, ugrid)
            
            # heavy particles: state
            for state in hp.state:
                quantity = state.emissivity
                name = f'{state.label} Emissivity ({p_name}) ' + '[$W m^{-3}$]'
                _add_aos_scalar_array_to_vtk_field_data(quantity, subset_idx, name, ugrid)


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
        name = f'Magnetic Field B{component_name} [$T$]'
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
        name = f'Electron Density Boundary Condition - {ggd.electrons.particles[subset_idx].identifier.name} ' + '[$m^{-3} s^{-1}$]'
        _add_aos_scalar_array_to_vtk_field_data(ggd.electrons.particles, subset_idx, name, ugrid)
    except IndexError:
        # can raise from name = "...", when subset_idx >= ggd.electrons.particles
        pass
    # - electrons/energy(i1)
    try:
        name = f'Electron Energy Boundary Condition - {ggd.electrons.energy[subset_idx].identifier.name} ' + '[$W m^{-3}$]'
        _add_aos_scalar_array_to_vtk_field_data(ggd.electrons.energy, subset_idx, name, ugrid)
    except IndexError:
        # can raise from name = "...", when subset_idx >= ggd.electrons.energy
        pass

    # ions
    for ion in ggd.ion:
        # - ion/particles
        quantity = ion.particles
        name = f'{ion.label} Density Boundary Condition ({quantity.identifier.name}) ' + '[$m^{-3} s^{-1}$]'
        _add_aos_scalar_array_to_vtk_field_data(quantity, subset_idx, name, ugrid)
        # - ion/energy
        quantity = ion.energy
        name = f'{ion.label} Energy Boundary Condition ({quantity.identifier.name}) ' + '[$W m^{-3}$]'
        _add_aos_scalar_array_to_vtk_field_data(quantity, subset_idx, name, ugrid)

        for state in ion.state:
            # - ion/state/particles    
            quantity = state.particles
            name = f'{state.label} Density Boundary Condition ({quantity.identifier.name}) ' + '[$m^{-3} s^{-1}$]'
            _add_aos_scalar_array_to_vtk_field_data(quantity, subset_idx, name, ugrid)
            # - ion/state/energy
            quantity = state.energy
            name = f'{state.label} Energy Boundary Condition ({quantity.identifier.name}) ' + '[$W m^{-3}$]'
            _add_aos_scalar_array_to_vtk_field_data(quantity, subset_idx, name, ugrid)

            
def read_wall(ids_obj, aos_index_values: dict,
              subset_idx: int, ugrid: vtkUnstructuredGrid) -> None:
    
    # ggd is at /wall/description_ggd(i1)/ggd(itime)
    description_ggd_idx = aos_index_values.get('DescriptionGgdIdx')
    time_idx = aos_index_values.get('TimeIdx')

    try:
        if description_ggd_idx != None:
            ggd = ids_obj.description_ggd[description_ggd_idx].ggd[time_idx]
        else:
            ggd = ids_obj.description_ggd[0].ggd[time_idx]
    except IndexError:
        return

    # Power density and temperature read if data is provided.
    # - power_density
    name = 'Power Density [$W m^{-2}$]'
    try:
        _add_aos_scalar_array_to_vtk_field_data(ggd.power_density, subset_idx, name, ugrid)
    except:
        print('No power density in IDS')
        pass

    # - temperature
    name = 'Temperature [$K$]'
    try:
        _add_aos_scalar_array_to_vtk_field_data(ggd.temperature, subset_idx, name, ugrid)
    except:
        print('No temperature in IDS.')
        pass

    
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
        name = 'Electric Field - LH polarized [$V/m$]'
        _add_aos_scalar_array_to_vtk_field_data(e_field.plus, subset_idx, name, ugrid)
        name = 'Electric Field - RH polarized [$V/m$]'
        _add_aos_scalar_array_to_vtk_field_data(e_field.minus, subset_idx, name, ugrid)
        name = 'Electric Field - Parallel to B [$V/m$]'
        _add_aos_scalar_array_to_vtk_field_data(e_field.parallel, subset_idx, name, ugrid)
        name = 'Electric Field - Normal to flux surface [$V/m$]'
        _add_aos_scalar_array_to_vtk_field_data(e_field.normal, subset_idx, name, ugrid)
        name = 'Electric Field - Tangential to flux surface [$V/m$]'
        _add_aos_scalar_array_to_vtk_field_data(e_field.bi_normal, subset_idx, name, ugrid)
    except IndexError:
        pass

    try:
        # WARNING: Are these units correct?
        b_field = ids_obj.coherent_wave[coherent_wave_idx].full_wave[time_idx].b_field
        name = 'Magnetic Field - Parallel to B [$V/m$]'
        _add_aos_scalar_array_to_vtk_field_data(b_field.parallel, subset_idx, name, ugrid)
        name = 'Magnetic Field - Normal to flux surface [$V/m$]'
        _add_aos_scalar_array_to_vtk_field_data(b_field.normal, subset_idx, name, ugrid)
        name = 'Magnetic Field - Tangential to flux surface [$V/m$]'
        _add_aos_scalar_array_to_vtk_field_data(b_field.bi_normal, subset_idx, name, ugrid)
    except IndexError:
        pass

    try:
        # WARNING: Are these units correct?
        name = 'Perpendicular Wave Vector [$V/m$]'
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

    # To see interpolation of point data on cells, just point data is necessary.
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

    # For wall IDS nodes, edges, cells, volumes in one partition.
    if subset_idx == -1:
        for i in range(4):
            if hasattr(aos_scalar_node[i], 'values') and len(aos_scalar_node[i].values):
                _add_scalar_array_to_vtk_field_data(aos_scalar_node[i].values, name, ugrid)
    else:
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
    :param subset_idx: an index into aos_vector_node
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

def _multi_add_aos_vector_to_vtk(vector_description: dict, subset_idx: int, ugrid: vtkUnstructuredGrid):
    """
    Add the dict of names and arrays in vector_description to the unstructured grid.
    :param vector_description: A dictionary of {'name, aos_vector_node} 
    Each aos_vector_node is a node with an array of vector values for each grid subset.
    Each name becomes the respective array name in VTK.
    :param subset_idx: an index into aos_vector_node
    :param ugrid: an unstructured grid instance
    :return: None
    """
    for name, node in vector_description.items():
        _add_aos_vector_array_to_vtk_field_data(node, subset_idx, name, ugrid)
    
