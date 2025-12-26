import pytest
from pydantic import ValidationError

from mlip_autopipec.data.models import DFTResult, TrainingConfig

def test_dft_result_successful():
    """Tests successful instantiation of a valid DFTResult."""
    result = DFTResult(
        total_energy_ev=-100.5,
        forces=[[0.0, 0.0, 0.1]],
        stress=[[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]],
        was_successful=True,
    )
    assert result.total_energy_ev == -100.5
    assert result.was_successful is True
    assert result.error_message is None

def test_dft_result_failed():
    """Tests successful instantiation of a failed DFTResult."""
    result = DFTResult(
        total_energy_ev=0.0,
        forces=[],
        stress=[],
        was_successful=False,
        error_message="SCF failed to converge",
    )
    assert result.was_successful is False
    assert result.error_message == "SCF failed to converge"

def test_dft_result_missing_field_raises_error():
    """Tests that missing a required field raises a ValidationError."""
    with pytest.raises(ValidationError):
        DFTResult(
            forces=[[0.0, 0.0, 0.1]],
            stress=[],
            was_successful=True,
        )

def test_training_config_instantiation():
    """Tests successful instantiation of a valid TrainingConfig."""
    config = TrainingConfig(
        model_type="mace",
        learning_rate=0.01,
        num_epochs=100,
        r_cut=5.0,
        delta_learn=True,
        baseline_potential="zbl",
    )
    assert config.model_type == "mace"
    assert config.delta_learn is True

def test_training_config_invalid_type_raises_error():
    """Tests that providing an incorrect type raises a ValidationError."""
    with pytest.raises(ValidationError):
        TrainingConfig(
            model_type="mace",
            learning_rate="fast",  # Should be a float
            num_epochs=100,
            r_cut=5.0,
            delta_learn=True,
            baseline_potential="zbl",
        )
