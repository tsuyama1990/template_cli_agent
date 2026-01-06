import pytest
from pydantic import ValidationError

from mlip_autopipec.common.pydantic_models import FullConfig, MDConfig, SamplingConfig, SystemConfig


@pytest.fixture
def valid_config_dict():
    """Provides a valid configuration dictionary."""
    return {
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "supercell_size": [3, 3, 3],
        },
        "exploration": {
            "temperature_k": 300.0,
            "pressure_gpa": 1.0,
            "time_step_fs": 1.0,
            "total_steps": 1000,
        },
        "sampling": {
            "method": "random",
            "num_samples": 10,
        },
    }

def test_full_config_valid(valid_config_dict):
    """Tests that a valid configuration dictionary parses correctly."""
    try:
        FullConfig(**valid_config_dict)
    except ValidationError as e:
        pytest.fail(f"Valid configuration failed to parse: {e}")

def test_system_config_missing_elements():
    """Tests validation error for missing 'elements' in SystemConfig."""
    with pytest.raises(ValidationError):
        SystemConfig(composition={"Fe": 1.0}, supercell_size=[1, 1, 1])

def test_md_config_invalid_temperature():
    """Tests validation error for non-positive temperature in MDConfig."""
    with pytest.raises(ValidationError):
        MDConfig(
            temperature_k=0,
            pressure_gpa=1.0,
            time_step_fs=1.0,
            total_steps=100
        )

def test_sampling_config_invalid_method():
    """Tests validation error for an invalid sampling method."""
    with pytest.raises(ValidationError):
        SamplingConfig(method="invalid_method", num_samples=10)

def test_full_config_extra_field(valid_config_dict):
    """Tests that an extra field in the config raises a ValidationError."""
    valid_config_dict["system"]["extra_parameter"] = "should_fail"
    with pytest.raises(ValidationError) as excinfo:
        FullConfig(**valid_config_dict)
    assert "extra_parameter" in str(excinfo.value)
