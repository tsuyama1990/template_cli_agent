# tests/unit/test_config_models.py

import pytest
from pydantic import ValidationError

from mlip_autopipec.configs.models import (
    DFTComputeConfig,
    MainConfig,
    MLIPTrainingConfig,
)


def test_dft_compute_config_valid():
    """Tests that a valid DFTComputeConfig model is accepted."""
    config = {
        "code": "quantum_espresso",
        "command": "mpirun -np 4 pw.x",
        "pseudopotentials": "SSSP_1.3_PBE_precision",
        "ecutwfc": 60.0,
        "ecutrho": 240.0,
        "kpoints_density": 2.5,
        "smearing": "mv",
        "degauss": 0.01,
    }
    model = DFTComputeConfig(**config)
    assert model.code == "quantum_espresso"
    assert model.ecutwfc == 60.0


def test_dft_compute_config_extra_field_fails():
    """Tests that an extra field in DFTComputeConfig raises a ValidationError."""
    config = {
        "code": "quantum_espresso",
        "command": "mpirun -np 4 pw.x",
        "pseudopotentials": "SSSP_1.3_PBE_precision",
        "ecutwfc": 60.0,
        "ecutrho": 240.0,
        "kpoints_density": 2.5,
        "smearing": "mv",
        "degauss": 0.01,
        "extra_field": "should_fail",
    }
    with pytest.raises(ValidationError):
        DFTComputeConfig(**config)


def test_mlip_training_config_valid():
    """Tests that a valid MLIPTrainingConfig model is accepted."""
    config = {
        "model_type": "ace",
        "r_cut": 5.0,
        "delta_learning": True,
        "base_potential": "lj_auto",
        "loss_weights": {"energy": 1.0, "force": 100.0},
    }
    model = MLIPTrainingConfig(**config)
    assert model.model_type == "ace"
    assert model.delta_learning is True


def test_mlip_training_config_missing_field_fails():
    """Tests that a missing field in MLIPTrainingConfig raises a ValidationError."""
    config = {
        "model_type": "ace",
        "r_cut": 5.0,
        "delta_learning": True,
        # "base_potential" is missing
        "loss_weights": {"energy": 1.0, "force": 100.0},
    }
    with pytest.raises(ValidationError):
        MLIPTrainingConfig(**config)


def test_main_config_valid():
    """Tests that a valid MainConfig model with nested models is accepted."""
    config = {
        "system": {"elements": ["Si"]},
        "dft_compute": {
            "code": "quantum_espresso",
            "command": "mpirun -np 4 pw.x",
            "pseudopotentials": "SSSP_1.3_PBE_precision",
            "ecutwfc": 60.0,
            "ecutrho": 240.0,
            "kpoints_density": 2.5,
            "smearing": "mv",
            "degauss": 0.01,
        },
        "mlip_training": {
            "model_type": "ace",
            "r_cut": 5.0,
            "delta_learning": True,
            "base_potential": "lj_auto",
            "loss_weights": {"energy": 1.0, "force": 100.0},
        },
    }
    model = MainConfig(**config)
    assert model.system["elements"] == ["Si"]
    assert model.dft_compute.code == "quantum_espresso"
    assert model.mlip_training.model_type == "ace"


def test_main_config_invalid_nested_fails():
    """Tests that an invalid nested model in MainConfig raises a ValidationError."""
    config = {
        "system": {"elements": ["Si"]},
        "dft_compute": {
            "code": "quantum_espresso",
            "command": "mpirun -np 4 pw.x",
            "pseudopotentials": "SSSP_1.3_PBE_precision",
            "ecutwfc": 60.0,
            # ecutrho is missing
            "kpoints_density": 2.5,
            "smearing": "mv",
            "degauss": 0.01,
        },
        "mlip_training": {
            "model_type": "ace",
            "r_cut": 5.0,
            "delta_learning": True,
            "base_potential": "lj_auto",
            "loss_weights": {"energy": 1.0, "force": 100.0},
        },
    }
    with pytest.raises(ValidationError):
        MainConfig(**config)
