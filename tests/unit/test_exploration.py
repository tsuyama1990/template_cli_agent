from __future__ import annotations

from unittest.mock import patch

import pytest
from ase.build import bulk, fcc111

from mlip_autopipec.exploration.engine import MDExplorer


@pytest.fixture
def default_config():
    return {"md_steps": 100, "temperature": 300.0, "mlip_model": "emt"}


def test_explorer_initialization(default_config):
    """Test that the MDExplorer initializes correctly with a given config."""
    explorer = MDExplorer(config=default_config)
    assert explorer.md_steps == 100
    assert explorer.temperature == 300.0
    assert explorer.mlip_model == "emt"


@patch("mlip_autopipec.exploration.engine.detect_vacuum")
@patch("mlip_autopipec.exploration.engine.NPT")
def test_ensemble_selection_bulk(mock_npt, mock_detect_vacuum, default_config):
    """Verify that NPT ensemble is selected for bulk systems."""
    mock_detect_vacuum.return_value = "bulk"
    atoms = bulk("Cu", "fcc", a=3.6)
    explorer = MDExplorer(config=default_config)
    explorer.explore([atoms])
    mock_detect_vacuum.assert_called_once_with(atoms)
    mock_npt.assert_called_once()


@patch("mlip_autopipec.exploration.engine.detect_vacuum")
@patch("mlip_autopipec.exploration.engine.Langevin")
def test_ensemble_selection_slab(mock_langevin, mock_detect_vacuum, default_config):
    """Verify that NVT (Langevin) ensemble is selected for slab systems."""
    mock_detect_vacuum.return_value = "slab"
    atoms = fcc111("Au", size=(2, 2, 3), vacuum=10.0)
    explorer = MDExplorer(config=default_config)
    explorer.explore([atoms])
    mock_detect_vacuum.assert_called_once_with(atoms)
    mock_langevin.assert_called_once()


def test_graceful_error_handling(default_config, caplog):
    """Test that the explorer handles simulation errors gracefully."""
    atoms = bulk("Cu", "fcc", a=3.6)
    explorer = MDExplorer(config=default_config)

    with patch(
        "mlip_autopipec.exploration.engine.LBFGS.run",
        side_effect=Exception("Test Crash"),
    ):
        results = explorer.explore([atoms])
        assert not results  # No frames should be collected
        assert "MD simulation failed" in caplog.text
        assert "Test Crash" in caplog.text
