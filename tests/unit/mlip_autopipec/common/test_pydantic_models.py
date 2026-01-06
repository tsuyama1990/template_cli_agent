# mypy: ignore-errors
from typing import Any, Dict

import pytest
from pydantic import ValidationError

from mlip_autopipec.common.pydantic_models import (
    FullConfig,
    MDConfig,
    SystemConfig,
)


def test_system_config_valid() -> None:
    """Tests that a valid SystemConfig model passes validation."""
    config = {
        "elements": ["Fe", "Pt"],
        "composition": {"Fe": 0.75, "Pt": 0.25},
        "supercell_size": [3, 3, 3],
    }
    model = SystemConfig(**config)
    assert model.elements == ["Fe", "Pt"]
    assert model.composition == {"Fe": 0.75, "Pt": 0.25}
    assert model.supercell_size == [3, 3, 3]


@pytest.mark.parametrize(
    ("invalid_data", "error_message"),
    [
        (
            {"elements": [], "composition": {"Fe": 1.0}, "supercell_size": [3, 3, 3]},
            "List should have at least 1 item",
        ),
        (
            {
                "elements": ["Fe"],
                "composition": {"Fe": 1.0},
                "supercell_size": [3, 3],
            },
            "List should have at least 3 items",
        ),
        (
            {
                "elements": ["Fe"],
                "composition": {"Fe": 1.0},
                "supercell_size": [3, 3, 3, 3],
            },
            "List should have at most 3 items",
        ),
    ],
)
def test_system_config_invalid(
    invalid_data: Dict[str, Any], error_message: str
) -> None:
    """Tests that invalid SystemConfig models raise ValidationError."""
    with pytest.raises(ValidationError) as excinfo:
        SystemConfig(**invalid_data)
    assert error_message in str(excinfo.value)


def test_md_config_valid() -> None:
    """Tests that a valid MDConfig model passes validation."""
    config = {"temperature_k": 300.0, "pressure_gpa": 1.0}
    model = MDConfig(**config)
    assert model.temperature_k == 300.0
    assert model.pressure_gpa == 1.0


def test_md_config_invalid_temperature() -> None:
    """Tests that a non-positive temperature raises a ValidationError."""
    config = {"temperature_k": 0, "pressure_gpa": 1.0}
    with pytest.raises(ValidationError) as excinfo:
        MDConfig(**config)
    assert "Input should be greater than 0" in str(excinfo.value)


def test_full_config_valid() -> None:
    """Tests that a valid nested FullConfig model passes validation."""
    config = {
        "system": {
            "elements": ["Si", "O"],
            "composition": {"Si": 0.33, "O": 0.67},
            "supercell_size": [4, 4, 4],
        },
        "exploration": {"temperature_k": 1000.0, "pressure_gpa": 0.0},
    }
    model = FullConfig(**config)
    assert model.system.elements == ["Si", "O"]
    assert model.exploration.temperature_k == 1000.0


def test_full_config_invalid_nested() -> None:
    """Tests that an invalid value in a nested model raises ValidationError."""
    config = {
        "system": {
            "elements": ["Si", "O"],
            "composition": {"Si": 0.33, "O": 0.67},
            "supercell_size": [4, 4, 4],
        },
        "exploration": {"temperature_k": -100.0, "pressure_gpa": 0.0},
    }
    with pytest.raises(ValidationError):
        FullConfig(**config)


def test_full_config_extra_fields_forbidden() -> None:
    """Tests that extra fields in the config raise a ValidationError."""
    config = {
        "system": {
            "elements": ["Si", "O"],
            "composition": {"Si": 0.33, "O": 0.67},
            "supercell_size": [4, 4, 4],
        },
        "exploration": {"temperature_k": 1000.0, "pressure_gpa": 0.0},
        "extra_field": "should not be allowed",
    }
    with pytest.raises(ValidationError) as excinfo:
        FullConfig(**config)
    assert "Extra inputs are not permitted" in str(excinfo.value)
