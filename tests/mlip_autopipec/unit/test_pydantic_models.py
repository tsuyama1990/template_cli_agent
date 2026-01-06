# tests/mlip_autopipec/unit/test_pydantic_models.py
# mypy: ignore-errors
"""Unit tests for the Pydantic configuration models."""

import pytest
from pydantic import ValidationError

from mlip_autopipec.common.pydantic_models import FullConfig, MDConfig, SystemConfig


def test_system_config_valid() -> None:
    """Test that a valid SystemConfig model parses correctly."""
    config = {
        "elements": ["Fe", "Pt"],
        "composition": {"Fe": 0.75, "Pt": 0.25},
        "supercell_size": [3, 3, 3],
    }
    system = SystemConfig(**config)
    assert system.elements == ["Fe", "Pt"]
    assert system.supercell_size == [3, 3, 3]


@pytest.mark.parametrize(
    "invalid_config",
    [
        {"elements": [], "composition": {"Fe": 1.0}, "supercell_size": [3, 3, 3]},  # elements empty
        {
            "elements": ["Fe"],
            "composition": {},
            "supercell_size": [3, 3, 3, 3],
        },  # supercell too long
        {"elements": ["Fe"], "composition": {}, "supercell_size": [3, 3]},  # supercell too short
    ],
)
def test_system_config_invalid(invalid_config) -> None:
    """Test that invalid SystemConfig models raise ValidationError."""
    with pytest.raises(ValidationError):
        SystemConfig(**invalid_config)


def test_md_config_valid() -> None:
    """Test that a valid MDConfig model parses correctly."""
    config = {"temperature_k": 300.0, "pressure_gpa": 1.0, "n_steps": 1000}
    md = MDConfig(**config)
    assert md.temperature_k == 300.0
    assert md.n_steps == 1000


@pytest.mark.parametrize(
    "invalid_config",
    [
        {"temperature_k": 0, "n_steps": 1000},  # temp must be > 0
        {"temperature_k": 300, "n_steps": 0},  # steps must be > 0
    ],
)
def test_md_config_invalid(invalid_config) -> None:
    """Test that invalid MDConfig models raise ValidationError."""
    with pytest.raises(ValidationError):
        MDConfig(**invalid_config)


def test_full_config_valid() -> None:
    """Test that a full, valid configuration is parsed correctly."""
    config = {
        "system": {
            "elements": ["Si"],
            "composition": {"Si": 1.0},
            "supercell_size": [2, 2, 2],
        },
        "exploration": {"temperature_k": 500, "n_steps": 5000},
        "sampling": {"method": "Random", "n_samples": 100},
        "db_path": "test.db",
    }
    full = FullConfig(**config)
    assert full.system.elements == ["Si"]
    assert full.sampling.n_samples == 100
    assert full.db_path == "test.db"


def test_full_config_extra_fields_invalid() -> None:
    """Test that extra fields in the config cause a validation error."""
    config = {
        "system": {
            "elements": ["Si"],
            "composition": {"Si": 1.0},
            "supercell_size": [2, 2, 2],
        },
        "exploration": {"temperature_k": 500, "n_steps": 5000},
        "sampling": {"n_samples": 100},
        "db_path": "test.db",
        "some_extra_field": "should_fail",
    }
    with pytest.raises(ValidationError):
        FullConfig(**config)
