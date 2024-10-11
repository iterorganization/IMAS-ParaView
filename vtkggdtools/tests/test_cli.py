import imaspy
from click.testing import CliRunner

from vtkggdtools.cli import cli, convert_ggd_to_vtk


def test_version():
    runner = CliRunner()

    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert "vtkggdtools version:" in result.output
    assert "IMASPy version:" in result.output
    assert "Default data dictionary version:" in result.output
    assert "Available data dictionary versions:" in result.output
    assert "Access Layer core version:" in result.output


def test_ggd2vtk(tmp_path, dummy_ids):
    uri = f"imas:hdf5?path={tmp_path}/testdb"
    with imaspy.DBEntry(uri, "w") as dbentry:
        dbentry.put(dummy_ids)

    runner = CliRunner()
    file_name = "test"
    output_path = tmp_path / file_name
    ids_name = dummy_ids.metadata.name
    output_str = str(output_path)
    occurrence = "0"
    args = [uri, ids_name, output_str, occurrence]
    result = runner.invoke(convert_ggd_to_vtk, args)
    assert result.exit_code == 0

    # Check if vtpc file and the directory containing vtu files exists
    output_dir = output_path / file_name
    assert output_dir.is_dir()
    assert output_dir.with_suffix(".vtpc").exists()

    if ids_name == "wall":
        vtu_file = output_dir / f"{file_name}_0_0.vtu"
        assert vtu_file.exists()
    else:
        for i in range(3):
            vtu_file = output_dir / f"{file_name}_{i}_0.vtu"
            assert vtu_file.exists()
