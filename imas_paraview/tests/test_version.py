import pytest
from packaging.version import Version

import imas_paraview


@pytest.mark.skip(reason="version checking doesnt work from docker container")
def test_version():
    version = imas_paraview.__version__
    assert Version(version) > Version("0")
    assert "unknown" not in version
