from packaging.version import Version

import vtkggdtools


def test_version():
    version = vtkggdtools.__version__
    assert Version(version) > Version("0")
    assert "unknown" not in version
