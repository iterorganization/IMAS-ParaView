import imas
import numpy as np
import pytest
from click import UsageError
from click.testing import CliRunner

from imas_paraview.cli import (
    cli,
    convert_ggd_to_vtk,
    parse_index,
    parse_time,
    parse_uri,
)


@pytest.mark.skip(reason="no IMAS-Core available")
def test_version():
    runner = CliRunner()

    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert "IMAS-ParaView version:" in result.output
    assert "IMAS-Python version:" in result.output
    assert "Default data dictionary version:" in result.output
    assert "Available data dictionary versions:" in result.output
    assert "Access Layer core version:" in result.output


def test_ggd2vtk(tmp_path, dummy_ids):
    uri = f"{tmp_path}/testdb.nc"
    with imas.DBEntry(uri, "w") as dbentry:
        dbentry.put(dummy_ids)

    runner = CliRunner()
    file_name = "test"
    output_path = tmp_path / file_name
    ids_name = dummy_ids.metadata.name
    output_str = str(output_path)
    occurrence = "0"
    uri_in = f"{uri}#{ids_name}:{occurrence}"
    args = [uri_in, output_str]
    result = runner.invoke(convert_ggd_to_vtk, args)
    assert result.exit_code == 0

    # The vtkXMLPartitionedDataSetCollectionWriter will output the following:
    # .
    # ├── test
    # │   ├── test_0_0.vtu
    # │   ├── test_1_0.vtu
    # │   ├── test_2_0.vtu
    # |   ├── ...
    # └── test.vtpc

    # Check if vtpc file and the directory containing vtu files exists
    assert output_path.is_dir()
    output_dir = output_path / f"{ids_name}_0"
    assert output_dir.with_suffix(".vtpc").exists()

    if ids_name == "wall":
        vtu_file = output_dir / f"{ids_name}_0_0_0.vtu"
        assert vtu_file.exists()
    else:
        for i in range(3):
            vtu_file = output_dir / f"{ids_name}_0_{i}_0.vtu"
            assert vtu_file.exists()


def test_parse_uri():
    uri_path, ids, occurrence = parse_uri("imas:hdf5?path=testdb#edge_profiles")
    assert uri_path == "imas:hdf5?path=testdb"
    assert ids == "edge_profiles"
    assert occurrence == 0
    uri_path, ids, occurrence = parse_uri("imas:hdf5?path=testdb2#edge_sources:1")
    assert uri_path == "imas:hdf5?path=testdb2"
    assert ids == "edge_sources"
    assert occurrence == 1
    with pytest.raises(UsageError):
        parse_uri("imas:hdf5?path=testdb")


def test_parse_index():
    # Assert correct usage
    valid_cases = [
        ("5", [5]),
        ("3:6", [3, 4, 5, 6]),
        ("2,4,6,7", [2, 4, 6, 7]),
        ("3:3", [3]),
        ("2,4,6,2", [2, 4, 6]),  # Duplicates should be removed
    ]
    for input_str, expected_output in valid_cases:
        assert parse_index(input_str) == expected_output

    # Assert invalid usage (raises UsageError)
    invalid_cases = [
        "",
        ",",
        ":",
        "a:b",
        "a,b,c",
        "3e",
        "5:3",
        "3::5",
        ":5",
        "5:",
        "-1",
        "-3:-1",
        "-5,2,-1",
        "1.2",
        "3:5.0",
        "2,3.0,4",
    ]
    for input_str in invalid_cases:
        with pytest.raises(UsageError):
            parse_index(input_str)


def test_parse_time():
    ids_times = np.array([1.1, 3.3, 6.6, 9.9, 10, 15])

    # Assert correct usage
    valid_cases = [
        ("3.3", [1]),
        ("3.3,6.6,10", [1, 2, 4]),
        ("3:11", [1, 2, 3, 4]),
        ("6.5", [1]),
        ("7,7.2,8.1,9.5", [2]),
        ("5, 14.19", [1, 4]),
        ("14", [4]),
        ("0:14", [0, 1, 2, 3, 4]),
        ("3:9.9", [1, 2, 3]),
        ("9.9:15", [3, 4, 5]),
        ("15", [5]),
        ("0", []),
    ]
    for input_str, expected_output in valid_cases:
        assert parse_time(ids_times, input_str) == expected_output

    # Assert invalid usage (raises UsageError)
    invalid_cases = [
        "",
        ",",
        ":",
        "a:b",
        "a,b,c",
        "3e",
        "5:3",
        "3::5",
        ":5",
        "5:",
    ]
    for input_str in invalid_cases:
        with pytest.raises(UsageError):
            parse_time(ids_times, input_str)
