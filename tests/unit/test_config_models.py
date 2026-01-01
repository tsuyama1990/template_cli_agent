import pytest
from pydantic import ValidationError

from mlip_autopipec.configs.models import (
    DFTComputeConfig,
    MLIPTrainingConfig,
    MainConfig,
)


def test_dft_compute_config_valid():
    config = {
        "code": "quantum_espresso",
        "command": "mpirun -np 4 pw.x",
        "pseudopotentials": "SSSP_1.3_PBE_precision",
        "ecutwfc": 60.0,
        "ecutrho": 240.0,
        "kpoints_density": 3.5,
        "smearing": "mv",
        "degauss": 0.01,
    }
    dft_config = DFTComputeConfig(**config)
    assert dft_config.code == "quantum_espresso"
    assert dft_config.ecutwfc == 60.0


def test_dft_compute_config_invalid_extra_field():
    config = {
        "code": "quantum_espresso",
        "command": "mpirun -np 4 pw.x",
        "pseudopotentials": "SSSP_1.3_PBE_precision",
        "ecutwfc": 60.0,
        "ecutrho": 240.0,
        "kpoints_density": 3.5,
        "smearing": "mv",
        "degauss": 0.01,
        "extra_field": "should_fail",
    }
    with pytest.raises(ValidationError):
        DFTComputeConfig(**config)


def test_mlip_training_config_valid():
    config = {
        "model_type": "ace",
        "r_cut": 5.0,
        "delta_learning": True,
        "base_potential": "lj_auto",
        "loss_weights": {"energy": 1.0, "force": 100.0},
    }
    mlip_config = MLIPTrainingConfig(**config)
    assert mlip_config.model_type == "ace"
    assert mlip_config.delta_learning is True


def test_main_config_valid():
    config = {
        "system": {"elements": ["Si"]},
        "dft_compute": {
            "code": "quantum_espresso",
            "command": "mpirun -np 4 pw.x",
            "pseudopotentials": "SSSP_1.3_PBE_precision",
            "ecutwfc": 60.0,
            "ecutrho": 240.0,
            "kpoints_density": 3.5,
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
    main_config = MainConfig(**config)
    assert main_config.system["elements"] == ["Si"]
    assert main_config.dft_compute.code == "quantum_espresso"
    assert main_config.mlip_training.model_type == "ace"
