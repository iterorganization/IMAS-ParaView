from vtk import vtkXMLPartitionedDataSetCollectionWriter

from vtkggdtools.convert import ggd_to_vtk
from vtkggdtools.util import get_grid_ggd


def test_ggd_to_vtk(dummy_ids, tmp_path):
    time_idx = 0
    output_file = tmp_path / "test.vtpc"
    vtk_object = ggd_to_vtk(dummy_ids)
    writer = vtkXMLPartitionedDataSetCollectionWriter()
    writer.SetInputData(vtk_object)
    writer.SetFileName(output_file)
    writer.Write()
    output_dir = output_file.parent / output_file.stem

    # Check if vtpc file and the directory exists
    assert output_file.exists()
    assert output_dir.is_dir()

    # Check if the vtu files exist
    grid_ggd = get_grid_ggd(dummy_ids, time_idx)
    num_subsets = len(grid_ggd.grid_subset)
    for n in range(num_subsets):
        vtu_file = output_file.stem + f"_{n}_0.vtu"
        assert output_dir / vtu_file
