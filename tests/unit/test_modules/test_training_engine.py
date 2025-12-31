"""Unit tests for the TrainingEngine."""

from unittest.mock import ANY, MagicMock, patch

import numpy as np
import pytest
from ase import Atoms

from mlip_autopipec.modules.training_engine import (
    RealAgnosticInteractionBlock,
    TrainingEngine,
)


@pytest.fixture
def mock_db():
    """Fixture for a mocked AseDB that returns a pre-calculated Atoms object."""
    db = MagicMock()
    ar_dimer = Atoms(
        "Ar2", positions=[[0, 0, 0], [0, 0, 3.8]], pbc=True, cell=np.diag([10, 10, 10])
    )
    ar_dimer.get_potential_energy = MagicMock(return_value=-0.02)
    ar_dimer.get_forces = MagicMock(return_value=np.array([[0, 0, 0.01], [0, 0, -0.01]]))
    db.get_atoms_by_state.return_value = [ar_dimer]
    return db


@pytest.fixture
def config(tmp_path):
    """Fixture for a sample configuration dictionary."""
    return {
        "cutoff": 5.0,
        "model_path": str(tmp_path / "test_model.pt"),
    }


@patch("mlip_autopipec.modules.training_engine.torch.save")
@patch("mlip_autopipec.modules.training_engine.MACE")
def test_training_engine_run(mock_mace, mock_torch_save, mock_db, config):
    """Test the full run method of the TrainingEngine."""
    engine = TrainingEngine(mock_db, config)
    engine.run()

    mock_db.get_atoms_by_state.assert_called_once_with("labelled")
    mock_mace.assert_called_once()
    init_kwargs = mock_mace.call_args[1]
    assert init_kwargs["interaction_cls"] == RealAgnosticInteractionBlock
    mock_torch_save.assert_called_once_with(ANY, config["model_path"])


def test_prepare_dataset(mock_db):
    """Test the data preparation step."""
    engine = TrainingEngine(mock_db, {})
    configs, dft_results = engine._prepare_dataset()

    assert len(configs) == 1
    assert len(dft_results) == 1
    assert hasattr(configs[0], "positions")
    assert "energy" in dft_results[0]
    assert dft_results[0]["energy"] == -0.02


def test_calculate_baseline(mock_db):
    """Test the baseline calculation step."""
    engine = TrainingEngine(mock_db, {})
    configs, _ = engine._prepare_dataset()
    baseline_results = engine._calculate_baseline(configs)

    assert len(baseline_results) == 1
    assert "energy" in baseline_results[0]
    assert "forces" in baseline_results[0]
    # For an Ar dimer near equilibrium, the LJ energy should be attractive (negative)
    assert baseline_results[0]["energy"] < 0


def test_calculate_delta(mock_db):
    """Test the delta calculation step."""
    engine = TrainingEngine(mock_db, {})
    configs, dft_results = engine._prepare_dataset()
    baseline_results = engine._calculate_baseline(configs)
    delta_results = engine._calculate_delta(dft_results, baseline_results)

    assert len(delta_results) == 1
    assert "energy" in delta_results[0]
    assert "forces" in delta_results[0]

    expected_delta_energy = dft_results[0]["energy"] - baseline_results[0]["energy"]
    assert np.isclose(delta_results[0]["energy"], expected_delta_energy)

    expected_delta_forces = dft_results[0]["forces"] - baseline_results[0]["forces"]
    np.testing.assert_allclose(delta_results[0]["forces"], expected_delta_forces)
