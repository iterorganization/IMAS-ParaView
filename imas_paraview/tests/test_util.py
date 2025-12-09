import imas
import numpy as np
import pytest

from imas_paraview.util import find_closest_indices, format_units


@pytest.mark.parametrize(
    "request_times, expected_indices",
    [
        ([2.2, 3.3], [1, 2]),
        ([4.4, 5.5], [3, 4]),
        ([1.5, 3.8, 5.9], [0, 2, 4]),
        ([10], [4]),
        ([0.1, 0.5, 0.9, 1.01], []),
        ([0.5, 1.01, 2.5, 4.7], [1, 3]),
        ([2.2, 2.2, 4.4, 4.4], [1, 1, 3, 3]),
    ],
)
def test_find_closest_indices(request_times, expected_indices):
    time_array = np.array([1.1, 2.2, 3.3, 4.4, 5.5])
    time_indices = find_closest_indices(request_times, time_array)
    assert time_indices == expected_indices


def test_format_units():
    cs = imas.IDSFactory("4.0.0").core_sources()
    cs.source.resize(1)
    cs.source[0].profiles_1d.resize(1)
    p1d = cs.source[0].profiles_1d[0]

    assert format_units(p1d.electrons.particles) == "[m⁻³·s⁻¹]"
    assert format_units(p1d.electrons.energy) == "[W·m⁻³]"
    assert format_units(p1d.momentum_phi) == "[kg·m⁻¹·s⁻²]"
