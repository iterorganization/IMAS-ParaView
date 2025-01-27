import subprocess
from pathlib import Path

import pytest

# Load XML scripts for integration testing
TEST_DIRECTORY = "./vtkggdtools/tests/integration_tests/xml_integration_tests"
TEST_SCRIPTS = list(Path(TEST_DIRECTORY).glob("*.xml"))


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


def test_xvfb_pass():
    """Test a passing integration test, by testing an empty XML."""
    test_passed = run_test(
        "./vtkggdtools/tests/integration_tests/xml_xvfb_tests/passing_test.xml"
    )
    assert test_passed


def test_xvfb_fail():
    """Test failing integration test."""
    test_passed = run_test(
        "./vtkggdtools/tests/integration_tests/xml_xvfb_tests/failing_test.xml"
    )
    assert not test_passed


@pytest.mark.parametrize("test_script", TEST_SCRIPTS)
def test_integration_scripts(test_script):
    """Parameterized test function for running integration tests."""
    test_passed = run_test(test_script)
    assert test_passed
