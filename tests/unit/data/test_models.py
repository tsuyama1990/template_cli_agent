import numpy as np
import pytest
from pydantic import ValidationError

from mlip_autopipec.data.models import DFTResults


def test_dft_results_success():
    """Test successful creation of a DFTResults object."""
    energy = -100.0
    forces = np.array([[0.1, 0.2, 0.3]])
    stress = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])

    results = DFTResults(energy=energy, forces=forces, stress=stress)

    assert results.energy == energy
    assert np.array_equal(results.forces, forces)
    assert np.array_equal(results.stress, stress)

def test_dft_results_list_conversion():
    """Test that lists are correctly converted to numpy arrays."""
    energy = -100.0
    forces_list = [[0.1, 0.2, 0.3]]
    stress_list = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]

    results = DFTResults(energy=energy, forces=forces_list, stress=stress_list)

    assert isinstance(results.forces, np.ndarray)
    assert isinstance(results.stress, np.ndarray)
    assert np.array_equal(results.forces, np.array(forces_list))
    assert np.array_equal(results.stress, np.array(stress_list))

def test_dft_results_invalid_type():
    """Test that a ValidationError is raised for invalid input types."""
    with pytest.raises(ValidationError):
        DFTResults(energy=-100.0, forces="not-an-array", stress=[1.0, 2.0])

    with pytest.raises(ValidationError):
        DFTResults(energy=-100.0, forces=[[0.1]], stress="not-an-array")
