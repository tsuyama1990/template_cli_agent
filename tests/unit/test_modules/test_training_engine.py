from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.calculators.singlepoint import SinglePointCalculator

from mlip_autopipec.config.models import TrainingParams
from mlip_autopipec.modules.d_training_engine import TrainingEngine


@pytest.fixture
def sample_training_params():
    """Provides a sample TrainingParams object for tests."""
    return TrainingParams(r_cut=5.0, delta_learning=True)


@pytest.fixture
def sample_labelled_dataset():
    """
    Creates a sample dataset of two labelled Atoms objects (H2 dimer).
    The DFT values are completely fabricated for this test.
    """
    atoms1 = Atoms("H2", positions=[[0, 0, 0], [0.75, 0, 0]])
    atoms1.calc = SinglePointCalculator(
        atoms1,
        energy=-30.0,
        forces=np.array([[0.1, 0, 0], [-0.1, 0, 0]]),
        stress=np.zeros(6),
    )

    atoms2 = Atoms("H2", positions=[[0, 0, 0], [1.5, 0, 0]])
    atoms2.calc = SinglePointCalculator(
        atoms2,
        energy=-15.0,
        forces=np.array([[0.05, 0, 0], [-0.05, 0, 0]]),
        stress=np.zeros(6),
    )
    return [atoms1, atoms2]


def test_prepare_data_delta_learning(sample_training_params, sample_labelled_dataset):
    """
    Tests that the _prepare_data method correctly calculates the 'delta'
    between the DFT labels and a baseline potential (Lennard-Jones).
    """
    engine = TrainingEngine(train_config=sample_training_params)

    # For this test, we are pretending the output is just a list of the modified atoms.
    # The actual implementation will convert this to a MACE-specific format.
    modified_dataset = engine._prepare_data(sample_labelled_dataset)

    # Use the same baseline potential to verify the delta.
    baseline = LennardJones()

    # --- Verification for the first structure ---
    original_atoms1 = sample_labelled_dataset[0]
    baseline.calculate(original_atoms1)
    expected_energy1 = original_atoms1.get_potential_energy() - baseline.get_potential_energy()
    expected_forces1 = original_atoms1.get_forces() - baseline.get_forces()

    # The test assumes the _prepare_data returns atoms with modified calculators
    # We will adjust this test once the real implementation is in place.
    # For now, this drives the implementation of the calculation itself.
    assert modified_dataset[0].get_potential_energy() == pytest.approx(expected_energy1)
    assert np.allclose(modified_dataset[0].get_forces(), expected_forces1)


@patch("mlip_autopipec.modules.d_training_engine.TrainingEngine._save_model")
@patch("mlip_autopipec.modules.d_training_engine.TrainingEngine._train_model")
@patch("mlip_autopipec.modules.d_training_engine.TrainingEngine._prepare_data")
def test_run_workflow(mock_prepare, mock_train, mock_save_model, sample_training_params, sample_labelled_dataset):
    """
    Tests that the main `run` method calls the internal methods in the correct sequence.
    """
    # Setup mocks
    mock_prepare.return_value = "prepared_data"
    mock_model = MagicMock()
    mock_train.return_value = mock_model

    engine = TrainingEngine(train_config=sample_training_params)
    engine.run(sample_labelled_dataset)

    # Assertions
    mock_prepare.assert_called_once_with(sample_labelled_dataset)
    mock_train.assert_called_once_with("prepared_data")
    mock_save_model.assert_called_once_with(mock_model)
