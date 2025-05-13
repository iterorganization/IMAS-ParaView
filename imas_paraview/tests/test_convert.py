import imas
import pytest
from imas.ids_defs import IDS_TIME_MODE_HOMOGENEOUS
from imas.ids_path import IDSPath

from imas_paraview.convert import Converter
from imas_paraview.io.read_ps import PlasmaStateReader
from imas_paraview.tests.fill_ggd import fill_ids
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
