import imaspy
import pytest
from imaspy.ids_path import IDSPath

from vtkggdtools.convert import _get_nearest_time_idx, convert_to_xml, ggd_to_vtk
from vtkggdtools.io.read_ps import PlasmaStateReader
from vtkggdtools.tests.fill_ggd import fill_ids
from vtkggdtools.util import get_grid_ggd


def test_ggd_to_vtk(dummy_ids):
    """Test if ggd_to_vtk converts all GGD arrays to VTK arrays."""
    ps = PlasmaStateReader(dummy_ids)
    scalar_paths, vector_paths, _, _ = ps.load_paths_from_ids(return_empty=True)
    ggd_names = names_from_ids(dummy_ids, scalar_paths, vector_paths)

    # Check if names of VTK object match the GGD array names
    vtk_object = ggd_to_vtk(dummy_ids)
    vtk_array_names = names_from_vtk(vtk_object)
    assert vtk_array_names == ggd_names


def test_ggd_to_vtk_index(dummy_ids_five_steps):
    """Test if ggd_to_vtk works with different time indices."""
    ps = PlasmaStateReader(dummy_ids_five_steps)
    scalar_paths, vector_paths, _, _ = ps.load_paths_from_ids(return_empty=True)
    ggd_names = names_from_ids(dummy_ids_five_steps, scalar_paths, vector_paths)

    # Check if names of VTK object match the GGD array names
    for time_idx in range(5):
        vtk_object = ggd_to_vtk(dummy_ids_five_steps, time_idx=time_idx)
        vtk_array_names = names_from_vtk(vtk_object)
        assert vtk_array_names == ggd_names


def test_ggd_to_vtk_time(dummy_ids_five_steps):
    """Test if ggd_to_vtk works with different times."""
    ps = PlasmaStateReader(dummy_ids_five_steps)
    scalar_paths, vector_paths, _, _ = ps.load_paths_from_ids(return_empty=True)
    ggd_names = names_from_ids(dummy_ids_five_steps, scalar_paths, vector_paths)

    # Check if names of VTK object match the GGD array names
    for time in range(5):
        vtk_object = ggd_to_vtk(dummy_ids_five_steps, time=time)
        vtk_array_names = names_from_vtk(vtk_object)
        assert vtk_array_names == ggd_names


def test_ggd_to_vtk_out_of_bounds(dummy_ids_five_steps):
    """Test if ggd_to_vtk fails when given an index which is not in the IDS."""
    time_idx = 6
    vtk_object = ggd_to_vtk(dummy_ids_five_steps, time_idx=time_idx)
    assert vtk_object is None


def test_ggd_to_vtk_subset():
    """Test for conversion of a subset of GGD arrays."""
    ids = imaspy.IDSFactory(version="3.41.0").new("edge_profiles")
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
    vtk_object = ggd_to_vtk(
        ids, scalar_paths=es_scalar_paths, vector_paths=es_vector_paths
    )
    vtk_array_names = names_from_vtk(vtk_object)
    ggd_names = names_from_ids(ids, es_scalar_paths, es_vector_paths)
    assert vtk_array_names == ggd_names


def test_ggd_to_vtk_subset_time_index(dummy_ids_five_steps):
    """Test if ggd_to_vtk returns None when given time and index values."""
    vtk_object = ggd_to_vtk(dummy_ids_five_steps, time=5, time_idx=6)
    assert vtk_object is None


def test_convert_to_xml_no_time(dummy_ids_five_steps, tmp_path):
    """Test if convert_to_xml raises an error if no time step is selected."""
    output_file = tmp_path / "test_.vtpc"
    with pytest.raises(RuntimeError):
        convert_to_xml(dummy_ids_five_steps, output_file)


def test_convert_to_xml_out_of_bounds(dummy_ids_five_steps, tmp_path):
    """Test if convert_to_xml fails when given an index which is not in the IDS."""
    output_file = tmp_path / "test_.vtpc"
    with pytest.raises(RuntimeError):
        convert_to_xml(dummy_ids_five_steps, output_file, index=6)


def test_convert_to_xml_index(dummy_ids_five_steps, tmp_path):
    """Test if convert_to_xml converts a single index."""
    for time_idx in range(5):
        output_file = tmp_path / f"test_{time_idx}"
        convert_to_xml(dummy_ids_five_steps, output_file, index=time_idx)
        assert_output_exists(dummy_ids_five_steps, time_idx, output_file)


def test_convert_to_xml_time(dummy_ids_five_steps, tmp_path):
    """Test if convert_to_xml converts a single time value."""
    for time in range(5):
        output_file = tmp_path / f"test_{time}"
        convert_to_xml(dummy_ids_five_steps, output_file, time=time)

        time_idx = _get_nearest_time_idx(dummy_ids_five_steps, time)
        assert_output_exists(dummy_ids_five_steps, time_idx, output_file)


def test_convert_to_xml_index_range(dummy_ids_five_steps, tmp_path):
    """Test if convert_to_xml converts a range of indices."""
    min_index = 1
    max_index = 3

    output_file = tmp_path / "test"
    convert_to_xml(
        dummy_ids_five_steps, output_file, index_range=[min_index, max_index]
    )

    for time_idx in range(min_index, max_index + 1):
        assert_output_exists(dummy_ids_five_steps, time_idx, output_file)


def test_convert_to_xml_time_range(dummy_ids_five_steps, tmp_path):
    """Test if convert_to_xml converts a range of time steps."""
    min_time = 1.1
    max_time = 3.2

    output_file = tmp_path / "test"
    convert_to_xml(dummy_ids_five_steps, output_file, time_range=[min_time, max_time])

    for time in range(int(min_time), int(max_time) + 1):
        time_idx = _get_nearest_time_idx(dummy_ids_five_steps, time)
        assert_output_exists(dummy_ids_five_steps, time_idx, output_file)


def test_convert_to_xml_all_times(dummy_ids_five_steps, tmp_path):
    """Test if convert_to_xml converts all times in the IDS."""

    output_file = tmp_path / "test"
    convert_to_xml(dummy_ids_five_steps, output_file, all_times=True)

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
