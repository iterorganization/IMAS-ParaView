import subprocess
from pathlib import Path

import pytest

# Load XML scripts for integration testing
ROOT_DIRECTORY = Path(__file__).resolve().parents[3]
TEST_DIRECTORY = (
    ROOT_DIRECTORY
    / "imas_paraview"
    / "tests"
    / "integration_tests"
    / "xml_integration_tests"
)
TEST_SCRIPTS = list(TEST_DIRECTORY.glob("*.xml"))


def run_test(script, timeout=120):
    """Run a single test script with ParaView and a timeout."""
    print(f"Running test script: {script}")

    try:
        result = subprocess.run(
            [
                "xvfb-run",
                "--error-file=/dev/stdout",
                "-a",
                "paraview",
                "--dr",
                # "-v",
                # "9",
                "--no-mpi",
                "--test-script",
                script,
                "--exit",
            ],
            timeout=timeout,
        )

        if result.returncode != 0:
            print(f"Test failed: {script} with return code {result.returncode}")
            return False

        return True

    except subprocess.TimeoutExpired:
        print(
            f"Test {script} timed out after {timeout} seconds. "
            "This indicates that ParaView encountered an issue and failed to complete "
            "the test, or the timeout value is too short for the test to finish."
        )
        return False
    except Exception as e:
        print(f"An error occurred: {e}")
        return False


@pytest.mark.integration
def test_xvfb_pass():
    """Test a passing integration test, by testing an empty XML."""
    test_passed = run_test(
        ROOT_DIRECTORY
        / "imas_paraview"
        / "tests"
        / "integration_tests"
        / "xml_xvfb_tests"
        / "passing_test.xml"
    )
    assert test_passed


@pytest.mark.integration
def test_xvfb_fail():
    """Test failing integration test."""
    test_passed = run_test(
        ROOT_DIRECTORY
        / "imas_paraview"
        / "tests"
        / "integration_tests"
        / "xml_xvfb_tests"
        / "failing_test.xml"
    )
    assert not test_passed


@pytest.mark.integration
@pytest.mark.parametrize(
    "test_script", [pytest.param(script, id=script.name) for script in TEST_SCRIPTS]
)
def test_integration_scripts(test_script):
    """Parameterized test function for running integration tests."""
    if "hdf5" in str(test_script) and not Path("hdf5_testdb").exists():
        pytest.skip("HDF5 testdb does not exist")
    if "mdsplus" in str(test_script) and not Path("mdsplus_testdb").exists():
        pytest.skip("MDSplus testdb does not exist")
    if "profiles_1d_mapper" in str(test_script):
        pytest.skip("Requires ITER data, not available on GH Actions")
    test_passed = run_test(test_script)
    assert test_passed
