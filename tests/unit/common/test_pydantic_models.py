# tests/unit/common/test_pydantic_models.py
import pytest
from pydantic import ValidationError

from mlip_autopipec.common.pydantic_models import (
    FullConfig,
    MDConfig,
    SystemConfig,
)


def test_system_config_valid() -> None:
    """Test that a valid SystemConfig passes validation."""
    config = {
        "elements": ["Fe", "Pt"],
        "composition": {"Fe": 0.75, "Pt": 0.25},
        "supercell_size": [3, 3, 3],
    }
    validated = SystemConfig(**config)  # type: ignore[arg-type]
    assert validated.elements == ["Fe", "Pt"]
    assert validated.composition == {"Fe": 0.75, "Pt": 0.25}


def test_system_config_invalid_composition_sum() -> None:
    """Test that validation fails if composition fractions do not sum to 1."""
    config = {
        "elements": ["Fe", "Pt"],
        "composition": {"Fe": 0.75, "Pt": 0.5},  # Sums to 1.25
        "supercell_size": [3, 3, 3],
    }
    with pytest.raises(ValidationError) as excinfo:
        SystemConfig(**config)  # type: ignore[arg-type]
    assert "Composition fractions must sum to 1.0" in str(excinfo.value)


def test_md_config_invalid_temperature() -> None:
    """Test that validation fails for a non-positive temperature."""
    config = {"temperature_k": -100, "pressure_gpa": 1.0, "timestep_fs": 1.0, "n_steps": 100}
    with pytest.raises(ValidationError):
        MDConfig(**config)  # type: ignore[arg-type]


def test_full_config_valid() -> None:
    """Test that a full, valid configuration can be instantiated."""
    config = {
        "system": {
            "elements": ["Si"],
            "composition": {"Si": 1.0},
            "supercell_size": [2, 2, 2],
        },
        "exploration": {
            "temperature_k": 300,
            "pressure_gpa": 0,
            "timestep_fs": 1.0,
            "n_steps": 100,
        },
        "sampling": {"method": "random", "n_samples": 50},
        "db_path": "test.db",
    }
    validated = FullConfig(**config)  # type: ignore[arg-type]
    assert validated.system.elements == ["Si"]
    assert validated.sampling.n_samples == 50
    assert validated.db_path == "test.db"


def test_full_config_forbids_extra_fields() -> None:
    """Test that extra fields are forbidden due to `extra='forbid'`."""
    config = {
        "system": {
            "elements": ["Si"],
            "composition": {"Si": 1.0},
            "supercell_size": [2, 2, 2],
        },
        "exploration": {
            "temperature_k": 300,
            "pressure_gpa": 0,
            "timestep_fs": 1.0,
            "n_steps": 100,
        },
        "sampling": {"method": "random", "n_samples": 50},
        "db_path": "test.db",
        "some_extra_field": "not_allowed",
    }
    with pytest.raises(ValidationError) as excinfo:
        FullConfig(**config)  # type: ignore[arg-type]
    assert "Extra inputs are not permitted" in str(excinfo.value)
