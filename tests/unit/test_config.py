import pytest
import numpy as np
from pydantic import ValidationError
from mlip_autopipec.config import DFTResult, DFTInputConfig, MLIPTrainingConfig

def test_dft_result_valid():
    """Tests successful creation and serialization of DFTResult."""
    forces = np.array([[0.0, 0.0, 0.1], [-0.0, -0.0, -0.1]])
    stress = np.random.rand(3, 3)

    result = DFTResult(energy=-10.5, forces=forces, stress=stress)

    assert result.energy == -10.5
    np.testing.assert_array_equal(result.forces, forces)

    # Test JSON serialization and deserialization (round-trip)
    json_data = result.model_dump_json()
    reloaded_result = DFTResult.model_validate_json(json_data)

    assert reloaded_result.energy == result.energy
    np.testing.assert_array_equal(reloaded_result.forces, result.forces)
    np.testing.assert_array_equal(reloaded_result.stress, result.stress)

def test_dft_result_invalid_forces():
    """Tests that a validation error is raised for incorrect forces type."""
    with pytest.raises(ValidationError) as excinfo:
        DFTResult(energy=1.0, forces="not_an_array", stress=np.zeros((3, 3)))
    assert "Value must be a list or a numpy array" in str(excinfo.value)

def test_dft_input_config_valid():
    """Tests successful creation of DFTInputConfig."""
    config = DFTInputConfig(
        pseudopotentials={"H": "H.UPF", "O": "O.UPF"},
        kpoints=(1, 1, 1),
        ecutwfc=50,
        control={"calculation": "scf"}
    )
    assert config.ecutwfc == 50
    assert config.control["calculation"] == "scf"

def test_dft_input_config_missing_field():
    """Tests validation error for missing required field."""
    with pytest.raises(ValidationError):
        DFTInputConfig(
            pseudopotentials={"H": "H.UPF"},
            kpoints=(1, 1, 1),
            # ecutwfc is missing
        )

def test_mlip_training_config_defaults():
    """Tests the default values in MLIPTrainingConfig."""
    config = MLIPTrainingConfig(r_cut=5.0, loss_weights={"energy": 1.0})
    assert config.model_type == "ACE"
