
import numpy as np
import pytest
from pydantic import ValidationError

from mlip_autopipec.config import DFTInputConfig, DFTResult, MLIPTrainingConfig, ModelType


def test_dft_input_config_valid():
    """Tests that a valid DFTInputConfig model is parsed correctly."""
    config = {
        "pseudopotentials": {"H": "H.UPF", "O": "O.UPF"},
        "kpoints": (1, 1, 1),
        "ecutwfc": 60,
        "control": {"calculation": "scf"},
    }
    model = DFTInputConfig(**config)
    assert model.pseudopotentials == config["pseudopotentials"]
    assert model.kpoints == config["kpoints"]


def test_dft_input_config_invalid():
    """Tests that a DFTInputConfig with missing fields raises a validation error."""
    with pytest.raises(ValidationError):
        DFTInputConfig(
            pseudopotentials={"H": "H.UPF"}, kpoints=(1, 1, 1), control={}
        )


def test_dft_result_valid():
    """Tests that a valid DFTResult model is parsed correctly."""
    result = {
        "energy": -1.0,
        "forces": np.array([[0.0, 0.0, 0.0]]),
        "stress": np.zeros((3, 3)),
    }
    model = DFTResult(**result)
    assert model.energy == result["energy"]
    np.testing.assert_array_equal(model.forces, result["forces"])


def test_dft_result_list_conversion():
    """Tests that lists are correctly converted to numpy arrays."""
    result = {
        "energy": -1.0,
        "forces": [[0.0, 0.0, 0.0]],
        "stress": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
    }
    model = DFTResult(**result)
    assert isinstance(model.forces, np.ndarray)
    assert isinstance(model.stress, np.ndarray)


def test_dft_result_invalid_type():
    """Tests that an invalid type for forces raises a TypeError."""
    with pytest.raises(ValidationError):
        DFTResult(energy=-1.0, forces="not-an-array", stress=np.zeros((3, 3)))


def test_dft_result_json_serialization():
    """Tests the JSON serialization and deserialization of DFTResult."""
    result = {
        "energy": -1.0,
        "forces": np.array([[0.1, 0.2, 0.3]]),
        "stress": np.ones((3, 3)),
    }
    model = DFTResult(**result)

    # Serialize to JSON
    json_data = model.model_dump_json()

    # Deserialize back to Pydantic model
    deserialized_model = DFTResult.model_validate_json(json_data)

    assert model.energy == deserialized_model.energy
    np.testing.assert_array_almost_equal(model.forces, deserialized_model.forces)
    np.testing.assert_array_almost_equal(model.stress, deserialized_model.stress)


def test_mlip_training_config_valid():
    """Tests that a valid MLIPTrainingConfig model is parsed correctly."""
    config = {
        "model_type": ModelType.ACE,
        "r_cut": 5.0,
        "loss_weights": {"energy": 1.0, "forces": 100.0},
    }
    model = MLIPTrainingConfig(**config)
    assert model.model_type == config["model_type"]


def test_mlip_training_config_invalid_enum():
    """Tests that an invalid model_type raises a validation error."""
    with pytest.raises(ValidationError):
        MLIPTrainingConfig(
            model_type="INVALID_MODEL",
            r_cut=5.0,
            loss_weights={"energy": 1.0},
        )
