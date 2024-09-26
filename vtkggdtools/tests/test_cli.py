from click.testing import CliRunner

from vtkggdtools.cli import convert_ggd_to_vtk, print_version


def test_version():
    runner = CliRunner()
    result = runner.invoke(print_version)
    assert result.exit_code == 0


def test_ggd2vtk(tmp_path):
    uri = "imas:hdf5?path=vtkggdtools/tests/test_ids/"
    runner = CliRunner()
    output_path = tmp_path / "test"
    output_str = str(output_path)
    ids_name = "edge_profiles"
    occurrence = "0"
    args = [uri, ids_name, output_str, occurrence]
    result = runner.invoke(convert_ggd_to_vtk, args)
    assert result.exit_code == 0

    # Check if vtpc file and the directory containing vtu files exists
    assert output_path.with_suffix(".vtpc").exists()
    assert output_path.is_dir()
    for i in range(3):
        vtu_file = output_path / f"test_{i}_0.vtu"
        assert vtu_file.exists()
