import numpy as np
import pytest
from pydantic import ValidationError

from mlip_autopipec.data.models import (
    Cycle01Config,
    DFTCompute,
    DFTResults,
    MLIPTraining,
)


# Tests for DFTCompute
def test_dft_compute_valid():
    """Tests successful validation of DFTCompute with valid data."""
    data = {
        "code": "quantum_espresso",
        "command": "pw.x",
        "pseudopotentials": {"Si": "Si.UPF"},
        "ecutwfc": 60.0,
        "ecutrho": 240.0,
        "kpoints_density": 3.0,
    }
    model = DFTCompute(**data)
    assert model.code == "quantum_espresso"
    assert model.ecutwfc == 60.0


def test_dft_compute_invalid_ecutrho():
    """Tests ValidationError when ecutrho is less than ecutwfc."""
    data = {
        "code": "quantum_espresso",
        "command": "pw.x",
        "pseudopotentials": {"Si": "Si.UPF"},
        "ecutwfc": 60.0,
        "ecutrho": 50.0,  # Invalid value
        "kpoints_density": 3.0,
    }
    with pytest.raises(ValidationError, match="ecutrho must be greater than or equal to ecutwfc"):
        DFTCompute(**data)


def test_dft_compute_extra_field():
    """Tests ValidationError on extra fields due to `extra='forbid'`."""
    data = {
        "code": "quantum_espresso",
        "command": "pw.x",
        "pseudopotentials": {"Si": "Si.UPF"},
        "ecutwfc": 60.0,
        "ecutrho": 240.0,
        "kpoints_density": 3.0,
        "extra_field": "should not be here",
    }
    with pytest.raises(ValidationError, match="Extra inputs are not permitted"):
        DFTCompute(**data)


# Tests for MLIPTraining
def test_mlip_training_valid():
    """Tests successful validation of MLIPTraining with valid data."""
    data = {
        "model_type": "ace",
        "r_cut": 5.0,
        "delta_learning": False,
        "loss_weights": {"energy": 1.0, "forces": 1.0},
    }
    model = MLIPTraining(**data)
    assert model.model_type == "ace"
    assert not model.delta_learning


def test_mlip_training_delta_learning_valid():
    """Tests successful validation of MLIPTraining with delta_learning enabled."""
    data = {
        "model_type": "ace",
        "r_cut": 5.0,
        "delta_learning": True,
        "base_potential": "some_potential",
        "loss_weights": {"energy": 1.0, "forces": 1.0},
    }
    model = MLIPTraining(**data)
    assert model.delta_learning
    assert model.base_potential == "some_potential"


def test_mlip_training_delta_learning_missing_base():
    """Tests ValidationError when delta_learning is True but base_potential is missing."""
    data = {
        "model_type": "ace",
        "r_cut": 5.0,
        "delta_learning": True,
        "loss_weights": {"energy": 1.0, "forces": 1.0},
    }
    with pytest.raises(
        ValidationError, match="base_potential must be specified for delta_learning"
    ):
        MLIPTraining(**data)


# Tests for DFTResults
def test_dft_results_valid():
    """Tests successful validation of DFTResults and numpy array conversion."""
    data = {
        "energy": -1.0,
        "forces": [[0.1, 0.2, 0.3]],
        "stress": [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
    }
    model = DFTResults(**data)
    assert model.energy == -1.0
    assert isinstance(model.forces, np.ndarray)
    assert isinstance(model.stress, np.ndarray)
    assert np.array_equal(model.forces, np.array([[0.1, 0.2, 0.3]]))


def test_dft_results_extra_field():
    """Tests ValidationError on extra fields for DFTResults."""
    data = {
        "energy": -1.0,
        "forces": [[0.0, 0.0, 0.0]],
        "stress": [[0.0] * 3] * 3,
        "another_field": 123,
    }
    with pytest.raises(ValidationError, match="Extra inputs are not permitted"):
        DFTResults(**data)


# Tests for Cycle01Config (integration of other models)
def test_cycle01_config_valid(tmp_path):
    """Tests successful validation of the top-level Cycle01Config."""
    db_file = tmp_path / "test.db"
    db_file.touch()

    data = {
        "dft_compute": {
            "code": "quantum_espresso",
            "command": "pw.x",
            "pseudopotentials": {"Si": "Si.UPF"},
            "ecutwfc": 60.0,
            "ecutrho": 240.0,
            "kpoints_density": 3.0,
        },
        "mlip_training": {
            "model_type": "ace",
            "r_cut": 5.0,
            "delta_learning": False,
            "loss_weights": {"energy": 1.0, "forces": 1.0},
        },
        "database_path": str(db_file),
    }
    model = Cycle01Config(**data)
    assert isinstance(model.dft_compute, DFTCompute)
    assert isinstance(model.mlip_training, MLIPTraining)


def test_cycle01_config_invalid_path():
    """Tests ValidationError for a non-existent database_path."""
    data = {
        "dft_compute": {
            "code": "quantum_espresso",
            "command": "pw.x",
            "pseudopotentials": {"Si": "Si.UPF"},
            "ecutwfc": 60.0,
            "ecutrho": 240.0,
            "kpoints_density": 3.0,
        },
        "mlip_training": {
            "model_type": "ace",
            "r_cut": 5.0,
            "delta_learning": False,
            "loss_weights": {"energy": 1.0, "forces": 1.0},
        },
        "database_path": "non_existent_file.db",
    }
    with pytest.raises(ValidationError, match="Path does not point to a file"):
        Cycle01Config(**data)
