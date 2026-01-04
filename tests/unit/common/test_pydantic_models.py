import pytest
from pydantic import ValidationError

from mlip_autopipec.common.pydantic_models import (
    FullConfig,
    MDConfig,
    SamplingConfig,
    SystemConfig,
)


def test_system_config_valid():
    config = {
        "elements": ["Fe", "Pt"],
        "composition": {"Fe": 0.5, "Pt": 0.5},
        "supercell_size": [3, 3, 3],
    }
    validated = SystemConfig(**config)
    assert validated.elements == ["Fe", "Pt"]
    assert validated.composition == {"Fe": 0.5, "Pt": 0.5}
    assert validated.supercell_size == (3, 3, 3)


def test_system_config_invalid_elements():
    with pytest.raises(ValidationError):
        SystemConfig(
            elements=[],
            composition={"Fe": 0.5, "Pt": 0.5},
            supercell_size=(3, 3, 3),
        )


def test_md_config_valid():
    config = {
        "temperature_k": 300.0,
        "pressure_gpa": 1.0,
        "time_step_fs": 1.0,
        "total_steps": 1000,
    }
    validated = MDConfig(**config)
    assert validated.temperature_k == 300.0
    assert validated.pressure_gpa == 1.0


def test_md_config_invalid_temperature():
    with pytest.raises(ValidationError):
        MDConfig(
            temperature_k=0,
            pressure_gpa=1.0,
            time_step_fs=1.0,
            total_steps=1000,
        )
    with pytest.raises(ValidationError):
        MDConfig(
            temperature_k=300,
            pressure_gpa=1.0,
            time_step_fs=0,
            total_steps=1000,
        )
    with pytest.raises(ValidationError):
        MDConfig(
            temperature_k=300,
            pressure_gpa=1.0,
            time_step_fs=1.0,
            total_steps=0,
        )


def test_sampling_config_valid():
    config = {"method": "random", "number_of_samples": 10}
    validated = SamplingConfig(**config)
    assert validated.number_of_samples == 10


def test_sampling_config_invalid_samples():
    with pytest.raises(ValidationError):
        SamplingConfig(method="random", number_of_samples=0)


def test_full_config_valid():
    config = {
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
        "sampling": {"method": "random", "number_of_samples": 10},
    }
    validated = FullConfig(**config)
    assert validated.system.elements == ["Fe", "Pt"]
    assert validated.exploration.temperature_k == 300.0
    assert validated.sampling.number_of_samples == 10


def test_full_config_extra_fields():
    config = {
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
        "sampling": {"method": "random", "number_of_samples": 10},
        "extra_field": "should_fail",
    }
    with pytest.raises(ValidationError):
        FullConfig(**config)
