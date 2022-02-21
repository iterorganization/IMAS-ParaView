from typing import Sequence, Union

from vtkmodules.vtkCommonCore import vtkDataArray, vtkIdList
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid, VTK_LINE, VTK_POLYGON, VTK_POLYHEDRON, vtkPointData, \
    vtkCellData
from vtkmodules.vtkFiltersCore import vtkPointDataToCellData
from vtkmodules.numpy_interface import dataset_adapter as dsa

from vtkggdtools.io.representables import GridSubsetRepresentable, GridGGDRepresentable


def write_plasma_state(ids_name: str, ids_obj, aos_index_values: dict, space_idx: int,
                       rep: GridGGDRepresentable, subset_rep: GridSubsetRepresentable, grid_ggd) -> None:
    """
    Write the plasma state attached to the vtkUnstructuredGrid in the form of vtkPointData/vtkCellData
    to the ggd IDS node under the given top level ids node.

    :param ids_name: name of the top level IDS. (Each IDS has different set of plasma state arrays.)
    :param ids_obj: an ids_obj of type ids_name
    :param aos_index_values: the values that shall be used to navigate the AoS and reach the scalar arrays.
    :param space_idx: an index into grid_ggd/space AoS
    :param subset_rep: representation of ids_obj/grid_ggd/grid_subset node for VTK index transforms.
    :param rep: representation of ids_obj/grid_ggd/space/objects_per_dimension for VTK index transforms.
    :param grid_ggd: the accompanying grid_ggd node.
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
        # _write_aos_scalar_node_from_vtk_field_data(name, subset_rep, rep, space_idx, ggd.particles, grid_ggd)

    elif ids_name == 'distributions':
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

    elif ids_name == 'edge_profiles':
        # ggd is at /edge_profiles/ggd(itime)
        time_idx = aos_index_values.get('TimeIdx')
        try:
            ggd = ids_obj.ggd[time_idx]
        except IndexError:
            return

        # electrons
        #  - temperature(i1)
        name = 'Electron Temperature (eV)'
        #  - density(i1)
        name = 'Electron Density (m^-3)'
        #  - density_fast(i1)
        name = 'Electron Density Fast (m^-3)'
        #  - pressure(i1)
        name = 'Electron Pressure (Pa)'
        #  - pressure_fast_perpendicular(i1)
        name = 'Electron Pressure Fast Perpendicular (Pa)'
        #  - pressure_fast_parallel(i1)
        name = 'Electron Pressure Fast Parallel (Pa)'
        #  - velocity(i1)
        name = 'Electron Velocity (m.s^-1)'
        #  - distribution_function(i1)
        name = 'Electron Distribution Function'

        # TODO: the remaining arrays. ion(). neutral(), ..

    elif ids_name == 'edge_sources':
        # ggd is at /edge_sources/ggd(itime)
        time_idx = aos_index_values.get('TimeIdx')
        try:
            ggd = ids_obj.ggd[time_idx]
        except IndexError:
            return

        # electrons
        #  - particles(i1)
        name = 'Electron Particle Density (m^-3.s^-1)'
        #  - density(i1)
        name = 'Electron Energy (W.m^-3)'

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
        ggd_idx = aos_index_values.get('GridIdx')

        if time_idx >= len(ids_obj.time_slice):
            ids_obj.time_slice.resize(time_idx + 1, keep=True)
        if ggd_idx >= len(ids_obj.time_slice[time_idx].ggd):
            ids_obj.time_slice[time_idx].ggd.resize(ggd_idx + 1, keep=True)
        ggd = ids_obj.time_slice[time_idx].ggd[ggd_idx]

        # - r(i1)
        name = 'Major Radius (m)'
        _write_aos_scalar_node_from_vtk_field_data(name, subset_rep, rep, space_idx, ggd.r, grid_ggd)
        # - z(i1)
        name = 'Height (m)'
        _write_aos_scalar_node_from_vtk_field_data(name, subset_rep, rep, space_idx, ggd.z, grid_ggd)
        # - psi(i1)
        name = 'Poloidal Flux (Wb)'
        _write_aos_scalar_node_from_vtk_field_data(name, subset_rep, rep, space_idx, ggd.psi, grid_ggd)
        # - phi(i1)
        name = 'Toroidal Flux (Wb)'
        _write_aos_scalar_node_from_vtk_field_data(name, subset_rep, rep, space_idx, ggd.phi, grid_ggd)
        # - theta(i1)
        name = 'Poloidal Angle (rad)'
        _write_aos_scalar_node_from_vtk_field_data(name, subset_rep, rep, space_idx, ggd.theta, grid_ggd)
        # - j_tor(i1)
        name = 'Toroidal Plasma Current Density (A.m^-2)'
        _write_aos_scalar_node_from_vtk_field_data(name, subset_rep, rep, space_idx, ggd.j_tor, grid_ggd)
        # - j_parallel(i1)
        name = 'Parallel Plasma Current Density (A.m^-2)'
        _write_aos_scalar_node_from_vtk_field_data(name, subset_rep, rep, space_idx, ggd.j_parallel, grid_ggd)
        # - b_field_r(i1)
        name = 'Magnetic Field Br (T)'
        _write_aos_scalar_node_from_vtk_field_data(name, subset_rep, rep, space_idx, ggd.b_field_r, grid_ggd)
        # - b_field_z(i1)
        name = 'Magnetic Field Bz (T)'
        _write_aos_scalar_node_from_vtk_field_data(name, subset_rep, rep, space_idx, ggd.b_field_z, grid_ggd)
        # - b_field_tor(i1)
        name = 'Magnetic Field Btor (T)'
        _write_aos_scalar_node_from_vtk_field_data(name, subset_rep, rep, space_idx, ggd.b_field_tor, grid_ggd)

    elif ids_name == 'mhd':
        # ggd is at /mhd/ggd(itime)
        time_idx = aos_index_values.get('TimeIdx')
        try:
            ggd = ids_obj.ggd[time_idx]
        except IndexError:
            return

        # TODO: all arrays

    elif ids_name == 'radiation':
        # ggd is at /radiation/process(i1)/ggd(itime)
        # TODO: ProcessIdx is not detected by generator.py and hence it is not present.
        process_idx = aos_index_values.get('ProcessIdx') or 0
        time_idx = aos_index_values.get('TimeIdx')
        try:
            ggd = ids_obj.process[process_idx].ggd[time_idx]
        except IndexError:
            return

        # TODO: all arrays

    elif ids_name == 'tf':
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

    elif ids_name == 'transport_solver_numerics':
        # ggd is at /transport_solver_numerics/boundary_conditions_ggd(itime)
        time_idx = aos_index_values.get('TimeIdx')
        try:
            ggd = ids_obj.boundary_conditions_ggd[time_idx]
        except IndexError:
            return

        # - current(i1)
        name = f'Current Boundary Condition - {ggd.current.identifier.name.capitalize()}'
        # - electrons/particles(i1)
        # - electrons/energy(i1)
        # TODO: ion arrays

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
        # - temperature
        name = 'Temperature  (K)'

    elif ids_name == 'waves':
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
            name = 'Electric Field - RH polarized (V.m^-1)'
            name = 'Electric Field - Parallel to B (V.m^-1)'
            name = 'Electric Field - Normal to flux surface (V.m^-1)'
            name = 'Electric Field - Tangential to flux surface (V.m^-1)'
        except IndexError:
            pass

        try:
            # WARNING: Are these units correct?
            b_field = ids_obj.coherent_wave[coherent_wave_idx].full_wave[time_idx].b_field
            name = 'Magnetic Field - Parallel to B (V.m^-1)'
            name = 'Magnetic Field - Normal to flux surface (V.m^-1)'
            name = 'Magnetic Field - Tangential to flux surface (V.m^-1)'
        except IndexError:
            pass

        try:
            # WARNING: Are these units correct?
            name = 'Perpendicular Wave Vector (V.m^-1)'
        except IndexError:
            pass


def _write_scalar_array_from_vtk_field_data(vtk_array: vtkDataArray, subset_idx: int, aos_scalar_node):
    if vtk_array is None:
        return
    as_np = dsa.vtkDataArrayToVTKArray(vtk_array)
    aos_scalar_node[subset_idx].values.resize(len(as_np))
    for i in range(len(as_np)):
        aos_scalar_node[subset_idx].values[i] = as_np[i]


def _interpolate_point_data_to_cell_data(name: str, ugrid: vtkUnstructuredGrid, cells: Sequence,
                                         cell_types: Sequence) -> Union[vtkDataArray, None]:
    dataset = ugrid.NewInstance()
    dataset.SetPoints(ugrid.GetPoints())
    dataset.GetPointData().AddArray(ugrid.GetPointData().GetArray(name))
    dataset.AllocateEstimate(ugrid.GetNumberOfCells(), 10)

    face_stream = vtkIdList()
    for cell, cell_type in zip(cells, cell_types):
        if cell_type != VTK_POLYHEDRON:
            dataset.InsertNextCell(cell_type, len(cell), cell)
        else:
            dataset.GetFaceStream(cell, face_stream)
            dataset.InsertNextCell(VTK_POLYHEDRON, face_stream)

    converter = vtkPointDataToCellData()
    converter.SetInputData(dataset)
    converter.SetProcessAllArrays(True)
    converter.Update(0)
    output = converter.GetOutput()

    return output.GetCellData().GetArray(name)


def _write_aos_scalar_node_from_vtk_field_data(name: str, subset_rep: GridSubsetRepresentable,
                                               rep: GridGGDRepresentable, space_idx: int,
                                               aos_scalar_node, grid_ggd):
    aos_scalar_node.resize(subset_rep.num_subsets)
    point_data: vtkPointData = subset_rep.ugrid.GetPointData()
    cell_data: vtkCellData = subset_rep.ugrid.GetCellData()

    if point_data.HasArray(name):
        print(f'Writing plasma state {name} for {4 + subset_rep.num_subsets} subsets ..')
        subset_idx = 0
        # deal with array from point data
        # write array from point data for 'nodes'
        arr = subset_rep.ugrid.GetPointData().GetArray(name)
        _write_scalar_array_from_vtk_field_data(arr, subset_idx, aos_scalar_node)
        aos_scalar_node[subset_idx].grid_index = grid_ggd.space[space_idx].identifier.index
        aos_scalar_node[subset_idx].grid_subset_index = grid_ggd.grid_subset[subset_idx].identifier.index

        # interpolate array from points to edges.
        subset_idx += 1
        cell_types = [VTK_LINE] * len(rep.edges)
        arr = _interpolate_point_data_to_cell_data(name, subset_rep.ugrid, list(rep.edges.keys()), cell_types)
        # write array from point data for 'edges'
        _write_scalar_array_from_vtk_field_data(arr, subset_idx, aos_scalar_node)
        aos_scalar_node[subset_idx].grid_index = grid_ggd.space[space_idx].identifier.index
        aos_scalar_node[subset_idx].grid_subset_index = grid_ggd.grid_subset[subset_idx].identifier.index

        # interpolate array from points to faces.
        subset_idx += 1
        cell_types = [VTK_POLYGON] * len(rep.faces)
        arr = _interpolate_point_data_to_cell_data(name, subset_rep.ugrid, list(rep.faces.keys()), cell_types)
        # write array from point data for 'faces'
        _write_scalar_array_from_vtk_field_data(arr, subset_idx, aos_scalar_node)
        aos_scalar_node[subset_idx].grid_index = grid_ggd.space[space_idx].identifier.index
        aos_scalar_node[subset_idx].grid_subset_index = grid_ggd.grid_subset[subset_idx].identifier.index

        # interpolate array from points to volumes.
        subset_idx += 1
        cell_types = [VTK_POLYHEDRON] * len(rep.volumes)
        arr = _interpolate_point_data_to_cell_data(name, subset_rep.ugrid, list(rep.volumes.keys()), cell_types)
        # write array from point data for 'volumes'
        _write_scalar_array_from_vtk_field_data(arr, subset_idx, aos_scalar_node)
        aos_scalar_node[subset_idx].grid_index = grid_ggd.space[space_idx].identifier.index
        aos_scalar_node[subset_idx].grid_subset_index = grid_ggd.grid_subset[subset_idx].identifier.index

        # write array from point data for other subsets interpolating when necessary
        for subset_idx in range(subset_idx + 1, subset_rep.num_subsets):
            # interpolate array from points to elements
            arr = _interpolate_point_data_to_cell_data(name, subset_rep.ugrid, subset_rep.element_list.get(subset_idx),
                                                       subset_rep.subset_cell_types.get(subset_idx))
            _write_scalar_array_from_vtk_field_data(arr, subset_idx, aos_scalar_node)
        print('Finished')

    if cell_data.HasArray(name):
        pass
        # deal with array from cell data
        # interpolate array from cells to nodes.
        # write array from cell data for 'nodes'
        # interpolate array from cells to edges.
        # write array from cell data for 'edges'
        # interpolate array from cells to faces.
        # write array from cell data for 'faces'
        # interpolate array from cells to volumes.
        # write array from cell data for 'volumes'
        # write array from cell data for other subsets interpolating when necessary

