import pytest
from pydantic import ValidationError

from mlip_autopipec.config.models import DFTParams, FullConfig, TrainingParams


def test_dft_params_success():
    """Tests successful creation of DFTParams."""
    params = DFTParams(
        command="pw.x",
        pseudopotentials={"Si": "Si.upf"},
        ecutwfc=60.0,
    )
    assert params.command == "pw.x"
    assert params.pseudopotentials == {"Si": "Si.upf"}
    assert params.ecutwfc == 60.0
    assert params.kpoints_density == 5000  # Check default


def test_dft_params_missing_required_field():
    """Tests that ValidationError is raised for missing required fields."""
    with pytest.raises(ValidationError) as excinfo:
        DFTParams(
            pseudopotentials={"Si": "Si.upf"},
            ecutwfc=60.0,
        )
    assert "command" in str(excinfo.value)


def test_training_params_success():
    """Tests successful creation of TrainingParams."""
    params = TrainingParams(r_cut=5.0)
    assert params.model_type == "ace"
    assert params.r_cut == 5.0
    assert params.delta_learning is True


def test_training_params_override_defaults():
    """Tests that default values can be overridden."""
    params = TrainingParams(r_cut=5.0, model_type="mace", delta_learning=False)
    assert params.model_type == "mace"
    assert params.delta_learning is False


def test_full_config_success():
    """Tests successful creation of the full configuration."""
    config_dict = {
        "dft_compute": {
            "command": "pw.x -in PREFIX.pwi > PREFIX.pwo",
            "pseudopotentials": {"Si": "Si.UPF"},
            "ecutwfc": 50.0,
        },
        "training": {"r_cut": 4.5},
        "ase_db_path": "test_data.db",
    }
    config = FullConfig(**config_dict)
    assert config.dft_compute.command == "pw.x -in PREFIX.pwi > PREFIX.pwo"
    assert config.training.r_cut == 4.5
    assert config.ase_db_path == "test_data.db"


def test_full_config_validation_error():
    """Tests that nested validation errors are caught."""
    config_dict = {
        "dft_compute": {
            "command": "pw.x",
            # Missing pseudopotentials
            "ecutwfc": 50.0,
        },
        "training": {"r_cut": 4.5},
        "ase_db_path": "test_data.db",
    }
    with pytest.raises(ValidationError) as excinfo:
        FullConfig(**config_dict)
    assert "dft_compute" in str(excinfo.value)
    assert "pseudopotentials" in str(excinfo.value)


def test_type_coercion_and_validation():
    """Tests that Pydantic correctly handles types."""
    with pytest.raises(ValidationError):
        DFTParams(
            command="pw.x",
            pseudopotentials={"Si": "Si.upf"},
            ecutwfc="not-a-float",  # Invalid type
        )

    # Pydantic can coerce string numbers to float, which is acceptable
    params = DFTParams(
        command="pw.x",
        pseudopotentials={"Si": "Si.upf"},
        ecutwfc="60.0",
    )
    assert isinstance(params.ecutwfc, float)
    assert params.ecutwfc == 60.0
