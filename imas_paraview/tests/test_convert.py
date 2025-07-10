import imas
import numpy as np
import pytest
import vtk
from imas import DBEntry
from imas.ids_defs import IDS_TIME_MODE_HETEROGENEOUS, IDS_TIME_MODE_HOMOGENEOUS
from imas.ids_path import IDSPath
from vtk.util.numpy_support import vtk_to_numpy

from imas_paraview.convert import Converter, InterpSettings
from imas_paraview.io.read_ps import PlasmaStateReader
from imas_paraview.tests.fill_ggd import fill_ids, fill_NxN_grid, fill_scalar_quantity
from imas_paraview.util import get_grid_ggd


def test_ggd_to_vtk(dummy_ids):
    """Test if ggd_to_vtk converts all GGD arrays to VTK arrays."""
    ps = PlasmaStateReader(dummy_ids)
    scalar_paths, vector_paths, _, _ = ps.load_paths_from_ids(return_empty=True)
    ggd_names = names_from_ids(dummy_ids, scalar_paths, vector_paths)

    # Check if names of VTK object match the GGD array names
    converter = Converter(dummy_ids)
    vtk_object = converter.ggd_to_vtk()

    vtk_array_names = names_from_vtk(vtk_object)
    assert vtk_array_names == ggd_names


def test_ggd_to_vtk_reference(tmp_path):
    """Test if ggd_to_vtk works with a reference grid."""

    with imas.DBEntry(f"{tmp_path}/test.nc", "w") as dbentry:
        ids = dbentry.factory.new("edge_profiles")
        fill_ids(ids)
        dbentry.put(ids)
        ids2 = dbentry.factory.new("edge_sources")
        ids2.ids_properties.homogeneous_time = IDS_TIME_MODE_HOMOGENEOUS
        ids2.grid_ggd.resize(1)
        ids2.time = [0]

        converter = Converter(ids2, dbentry=dbentry)
        assert converter.ggd_to_vtk() is None

        ids2.grid_ggd[0].path = "#edge_profiles/grid_ggd(1)"
        converter = Converter(ids2, dbentry=dbentry)
        assert converter.ggd_to_vtk() is not None


def test_ggd_to_vtk_index(dummy_ids_five_steps):
    """Test if ggd_to_vtk works with different time indices."""
    ps = PlasmaStateReader(dummy_ids_five_steps)
    scalar_paths, vector_paths, _, _ = ps.load_paths_from_ids(return_empty=True)
    ggd_names = names_from_ids(dummy_ids_five_steps, scalar_paths, vector_paths)

    converter = Converter(dummy_ids_five_steps)
    # Check if names of VTK object match the GGD array names
    for time_idx in range(5):
        vtk_object = converter.ggd_to_vtk(time_idx=time_idx)
        vtk_array_names = names_from_vtk(vtk_object)
        assert vtk_array_names == ggd_names


def test_ggd_to_vtk_time(dummy_ids_five_steps):
    """Test if ggd_to_vtk works with different times."""
    ps = PlasmaStateReader(dummy_ids_five_steps)
    scalar_paths, vector_paths, _, _ = ps.load_paths_from_ids(return_empty=True)
    ggd_names = names_from_ids(dummy_ids_five_steps, scalar_paths, vector_paths)

    converter = Converter(dummy_ids_five_steps)
    # Check if names of VTK object match the GGD array names
    for time in range(5):
        vtk_object = converter.ggd_to_vtk(time=time)
        vtk_array_names = names_from_vtk(vtk_object)
        assert vtk_array_names == ggd_names


def test_ggd_to_vtk_out_of_bounds(dummy_ids_five_steps):
    """Test if ggd_to_vtk fails when given an index which is not in the IDS."""
    time_idx = 6
    converter = Converter(dummy_ids_five_steps)
    vtk_object = converter.ggd_to_vtk(time_idx=time_idx)
    assert vtk_object is None


def test_ggd_to_vtk_heterogeneous():
    ids = imas.IDSFactory(version="3.41.0").new("edge_profiles")
    ids.ids_properties.homogeneous_time = IDS_TIME_MODE_HETEROGENEOUS
    # 5 data points, two grids:
    ids.time = [0.0, 0.1, 0.2, 0.3, 0.4]
    ids.ggd.resize(5)
    for ggd, time in zip(ids.ggd, ids.time):
        ggd.time = time

    ids.grid_ggd.resize(2)
    ids.grid_ggd[0].time = 0
    ids.grid_ggd[1].time = 0.2

    num_vertices, num_edges, num_faces = fill_NxN_grid(ids.grid_ggd[0], 2)
    for i in range(2):
        fill_scalar_quantity(ids.ggd[i].zeff, num_vertices, num_edges, num_faces)
    num_vertices, num_edges, num_faces = fill_NxN_grid(ids.grid_ggd[1], 3)
    for i in range(2, 5):
        fill_scalar_quantity(ids.ggd[i].zeff, num_vertices, num_edges, num_faces)

    converter = Converter(ids)
    for i in range(5):
        vtk_object = converter.ggd_to_vtk(time_idx=i)
        assert vtk_object is not None

        num_faces = vtk_object.GetPartition(2, 0).GetNumberOfCells()
        assert num_faces == (1 if i < 2 else 4)


def test_ggd_to_vtk_subset():
    """Test for conversion of a subset of GGD arrays."""
    ids = imas.IDSFactory(version="3.41.0").new("edge_profiles")
    fill_ids(ids)

    # Convert subset of filled paths
    es_scalar_paths = [
        IDSPath("ggd/ion/state/density"),
        IDSPath("ggd/electrons/temperature"),
        IDSPath("ggd/phi_potential"),
    ]
    es_vector_paths = [
        IDSPath("ggd/neutral/state/velocity"),
        IDSPath("ggd/ion/velocity"),
        IDSPath("ggd/e_field"),
    ]
    converter = Converter(ids)
    vtk_object = converter.ggd_to_vtk(
        scalar_paths=es_scalar_paths, vector_paths=es_vector_paths
    )
    vtk_array_names = names_from_vtk(vtk_object)
    ggd_names = names_from_ids(ids, es_scalar_paths, es_vector_paths)
    assert vtk_array_names == ggd_names


def test_ggd_to_vtk_subset_time_index(dummy_ids_five_steps):
    """Test if ggd_to_vtk returns None when given time and index values."""
    converter = Converter(dummy_ids_five_steps)
    vtk_object = converter.ggd_to_vtk(time=5, time_idx=6)
    assert vtk_object is None


def test_ggd_to_vtk_solps(test_data_dir):
    with DBEntry(test_data_dir / "iter_db-123364-1.nc", "r") as entry:
        ids = entry.get("edge_profiles", autoconvert=False)
        converter = Converter(ids)
        vtk_object = converter.ggd_to_vtk()
        num_pds = vtk_object.GetNumberOfPartitionedDataSets()
        grid_subsets = ids.grid_ggd[0].grid_subset
        assert num_pds == len(grid_subsets)
        # Check if names of subsets match partition names
        for i in range(num_pds):
            name = vtk_object.GetMetaData(i).Get(vtk.vtkCompositeDataSet.NAME())
            assert name == grid_subsets[i].identifier.name

        # Check points array
        pd = vtk_object.GetPartitionedDataSet(0)
        vtk_grid = pd.GetPartition(0)
        vtk_point_data = vtk_grid.GetPointData()
        elec_temp = ids.ggd[0].electrons.temperature[0].values
        vtk_elec = vtk_point_data.GetArray("Electrons Temperature [eV]")
        np_vtk_elec = vtk_to_numpy(vtk_elec)
        assert np.array_equal(elec_temp, np_vtk_elec)

        # Check cells array
        pd = vtk_object.GetPartitionedDataSet(4)
        vtk_grid = pd.GetPartition(0)
        vtk_face_data = vtk_grid.GetCellData()
        vtk_face_array = vtk_face_data.GetArray("Electrons Temperature [eV]")
        np_vtk_face = vtk_to_numpy(vtk_face_array)
        elec_temp_face = ids.ggd[0].electrons.temperature[4].values
        assert np.array_equal(elec_temp_face, np_vtk_face)


def test_ggd_to_vtk_jorek(test_data_dir):
    with DBEntry(test_data_dir / "iter_dis-113112-1.nc", "r") as entry:
        ids = entry.get("plasma_profiles", autoconvert=False)
        converter = Converter(ids)

        plane_config = InterpSettings(n_plane=3, phi_start=0, phi_end=180)
        for t in range(3):
            vtk_object = converter.ggd_to_vtk(time_idx=t, plane_config=plane_config)
            pd = vtk_object.GetPartitionedDataSet(0)
            vtk_grid = pd.GetPartition(0)
            vtk_point_data = vtk_grid.GetPointData()
            vtk_elec = vtk_point_data.GetArray("Electrons Temperature [eV]")
            assert vtk_elec is not None


def assert_cache(converter, hits, misses):
    """Assert that the cache of a converter has a certain number of hits and misses."""
    cache = converter.get_grids.cache_info()
    assert cache.hits == hits
    assert cache.misses == misses


def test_ggd_to_vtk_grid_caching(dummy_ids):
    """Test if ggd_to_vtk caches grids if the same time step is provided."""
    converter = Converter(dummy_ids)
    vtk_object1 = converter.ggd_to_vtk()
    assert_cache(converter, 0, 1)
    vtk_object2 = converter.ggd_to_vtk()
    assert_cache(converter, 1, 1)
    vtk_object3 = converter.ggd_to_vtk()
    assert_cache(converter, 2, 1)

    assert (
        names_from_vtk(vtk_object1)
        == names_from_vtk(vtk_object2)
        == names_from_vtk(vtk_object3)
    )


def test_ggd_to_vtk_grid_caching_time_dependent(dummy_ids_five_steps):
    """Test if ggd_to_vtk caches grids if the same time step is provided."""
    converter = Converter(dummy_ids_five_steps)
    for time_idx in range(5):
        _ = converter.ggd_to_vtk(time_idx=time_idx)
        assert_cache(converter, 0, time_idx + 1)

    for time_idx in range(5):
        _ = converter.ggd_to_vtk(time_idx=time_idx)
        assert_cache(converter, time_idx + 1, 5)


def test_convert_to_xml(dummy_ids_five_steps, tmp_path):
    """Test if convert_to_xml converts a single index."""
    output_file = tmp_path / "test"
    converter = Converter(dummy_ids_five_steps)
    converter.write_to_xml(output_file)
    assert_output_exists(dummy_ids_five_steps, 0, output_file)


def test_convert_to_xml_out_of_bounds(dummy_ids_five_steps, tmp_path):
    """Test if convert_to_xml fails when given an index which is not in the IDS."""
    output_file = tmp_path / "test_.vtpc"
    converter = Converter(dummy_ids_five_steps)
    with pytest.raises(RuntimeError):
        converter.write_to_xml(output_file, [6])


def test_convert_to_xml_index(dummy_ids_five_steps, tmp_path):
    """Test if convert_to_xml converts a single index."""
    converter = Converter(dummy_ids_five_steps)

    for time_idx in range(5):
        output_file = tmp_path / f"test_{time_idx}"
        converter.write_to_xml(output_file, [time_idx])
        assert_output_exists(dummy_ids_five_steps, time_idx, output_file)


def test_convert_to_xml_index_list(dummy_ids_five_steps, tmp_path):
    """Test if convert_to_xml converts an index list."""

    output_file = tmp_path / "test"
    time_idx = [0, 1, 2, 3, 4]
    converter = Converter(dummy_ids_five_steps)
    converter.write_to_xml(output_file, time_idx)
    for time_idx in range(5):
        assert_output_exists(dummy_ids_five_steps, time_idx, output_file)


def assert_output_exists(ids, time_idx, output_file):
    """Assert the VTK object is correctly written to disk."""
    output_dir = output_file.parent / output_file.stem

    # Check if vtpc file and the directory exists
    assert output_file.exists()
    assert output_dir.is_dir()

    # Check if the vtu files exist
    grid_ggd = get_grid_ggd(ids, time_idx)
    num_subsets = len(grid_ggd.grid_subset)
    for n in range(num_subsets):
        vtu_file = output_file.stem + f"_{n}_0.vtu"
        assert output_dir / vtu_file


def names_from_ids(ids, scalar_paths, vector_paths):
    """Convert the names of the GGD array paths to the names given in Paraview."""
    ps = PlasmaStateReader(ids)
    ps.load_arrays_from_path(0, scalar_paths, vector_paths)
    ggd_names = set()
    for array in ps.scalar_array_list + ps.vector_array_list:
        name = ps._create_name_with_units(array)
        ggd_names.add(name)
    return ggd_names


def names_from_vtk(vtk_object):
    """Extract the array names from the VTK object."""
    n_partds = vtk_object.GetNumberOfPartitionedDataSets()
    array_names = set()
    for i in range(n_partds):
        part_ds = vtk_object.GetPartitionedDataSet(i)
        n_part = part_ds.GetNumberOfPartitions()
        for j in range(n_part):
            part = part_ds.GetPartition(j)
            cell_data = part.GetCellData()
            num_arrays = cell_data.GetNumberOfArrays()
            for k in range(num_arrays):
                array_name = cell_data.GetArrayName(k)
                array_names.add(array_name)
    return array_names
