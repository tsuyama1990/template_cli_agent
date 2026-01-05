"""Unit tests for the Pydantic configuration models."""
import pytest
from pydantic import ValidationError

from mlip_autopipec.config.models import FullConfig


@pytest.fixture
def valid_config_dict() -> dict:
    """Return a dictionary representing a valid configuration."""
    return {
        "project_name": "test_project",
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 0.8},
    }


def test_valid_config(valid_config_dict: dict) -> None:
    """Test that a valid configuration dictionary is parsed correctly."""
    try:
        config = FullConfig(**valid_config_dict)
        assert config.project_name == "test_project"
        assert config.system.lattice == "fcc"
        assert config.sampling.fraction == 0.8
    except ValidationError as e:
        pytest.fail(f"Valid configuration failed validation: {e}")


def test_invalid_composition_sum(valid_config_dict: dict) -> None:
    """Test that a validation error is raised for invalid composition sum."""
    valid_config_dict["system"]["composition"] = {"Fe": 0.6, "Pt": 0.5}
    with pytest.raises(ValidationError, match="sum of compositions must be 1.0"):
        FullConfig(**valid_config_dict)


def test_mismatched_elements(valid_config_dict: dict) -> None:
    """Test that a validation error is raised for mismatched elements."""
    valid_config_dict["system"]["elements"] = ["Fe", "Au"]
    with pytest.raises(ValidationError, match="must exactly match the elements list"):
        FullConfig(**valid_config_dict)


@pytest.mark.parametrize(
    ("field", "value", "match"),
    [
        ("temperature", -10.0, "Input should be greater than 0"),
        ("num_structures", 0, "Input should be greater than 0"),
        ("fraction", 1.1, "Input should be less than or equal to 1"),
        ("fraction", 0.0, "Input should be greater than 0"),
    ],
)
def test_invalid_numeric_fields(
    valid_config_dict: dict, field: str, value: float | int, match: str
) -> None:
    """Test validation for out-of-bounds numeric fields."""
    if field in valid_config_dict["system"]:
        valid_config_dict["system"][field] = value
    elif field in valid_config_dict["exploration"]:
        valid_config_dict["exploration"][field] = value
    else:
        valid_config_dict["sampling"][field] = value

    with pytest.raises(ValidationError, match=match):
        FullConfig(**valid_config_dict)


def test_invalid_type(valid_config_dict: dict) -> None:
    """Test that a validation error is raised for incorrect data types."""
    valid_config_dict["system"]["num_structures"] = "ten"
    with pytest.raises(ValidationError, match="Input should be a valid integer"):
        FullConfig(**valid_config_dict)


def test_extra_field_forbidden(valid_config_dict: dict) -> None:
    """Test that extra fields are not allowed."""
    valid_config_dict["system"]["extra_field"] = "some_value"
    with pytest.raises(ValidationError, match="Extra inputs are not permitted"):
        FullConfig(**valid_config_dict)
