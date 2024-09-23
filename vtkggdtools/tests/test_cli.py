from click.testing import CliRunner

from vtkggdtools.cli import print_version


def test_vtkggdtools_version():
    runner = CliRunner()
    result = runner.invoke(print_version)
    assert result.exit_code == 0
