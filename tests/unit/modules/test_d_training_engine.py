from unittest.mock import patch

import numpy as np
import pytest
import torch
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from mlip_autopipec.data.models import MLIPTraining
from mlip_autopipec.modules.d_training_engine import TrainingEngine


@pytest.fixture
def labeled_atoms_list():
    """Fixture for a list of labeled Atoms objects."""
    atoms1 = Atoms('H', positions=[(0, 0, 0)])
    calc1 = SinglePointCalculator(
        atoms1, energy=-1.0, forces=np.array([[0.1, 0.1, 0.1]])
    )
    atoms1.calc = calc1

    atoms2 = Atoms('H', positions=[(0, 0, 1)])
    calc2 = SinglePointCalculator(
        atoms2, energy=-2.0, forces=np.array([[0.2, 0.2, 0.2]])
    )
    atoms2.calc = calc2

    return [atoms1, atoms2]

@pytest.fixture
def training_config_no_delta():
    """Fixture for training config without delta learning."""
    return MLIPTraining(
        model_type="ace", r_cut=5.0, delta_learning=False,
        loss_weights={'energy': 1.0, 'forces': 1.0}
    )

@pytest.fixture
def training_config_with_delta():
    """Fixture for training config with delta learning."""
    return MLIPTraining(
        model_type="ace", r_cut=5.0, delta_learning=True,
        base_potential="lj_auto", loss_weights={'energy': 1.0, 'forces': 1.0}
    )

def test_training_engine_no_delta_data_preparation(
    training_config_no_delta, labeled_atoms_list
):
    """Test the data preparation step of the engine without delta learning."""
    engine = TrainingEngine(training_config_no_delta)
    prepared_data = engine._prepare_training_data(labeled_atoms_list)

    assert len(prepared_data) == 2
    assert prepared_data[0].get_potential_energy() == pytest.approx(-1.0)
    assert prepared_data[1].get_potential_energy() == pytest.approx(-2.0)

@patch('ase.calculators.lj.LennardJones')
def test_training_engine_with_delta_data_preparation(
    mock_lj, training_config_with_delta, labeled_atoms_list
):
    """Test the data preparation step of the engine with delta learning."""
    mock_lj_instance = mock_lj.return_value
    mock_lj_instance.get_potential_energy.side_effect = [-0.5, -0.8]
    mock_lj_instance.get_forces.side_effect = [
        np.array([[0.05, 0.05, 0.05]]),
        np.array([[0.1, 0.1, 0.1]])
    ]

    engine = TrainingEngine(training_config_with_delta)
    prepared_data = engine._prepare_training_data(labeled_atoms_list)

    assert prepared_data[0].get_potential_energy() == pytest.approx(-0.5)
    assert prepared_data[1].get_potential_energy() == pytest.approx(-1.2)
    expected_forces1 = np.array([[0.1, 0.1, 0.1]]) - np.array([[0.05, 0.05, 0.05]])
    expected_forces2 = np.array([[0.2, 0.2, 0.2]]) - np.array([[0.1, 0.1, 0.1]])
    assert np.allclose(prepared_data[0].get_forces(), expected_forces1)
    assert np.allclose(prepared_data[1].get_forces(), expected_forces2)

@patch('mlip_autopipec.modules.d_training_engine.MACE')
@patch('mlip_autopipec.modules.d_training_engine.torch.save')
@patch('mlip_autopipec.modules.d_training_engine.TrainingEngine._prepare_training_data')
def test_training_engine_execute_flow(
    mock_prepare_data, mock_torch_save, mock_mace,
    training_config_no_delta, labeled_atoms_list
):
    """Test the overall execution flow of the TrainingEngine."""
    mock_model_instance = mock_mace.return_value
    mock_model_instance.parameters.return_value = [
        torch.nn.Parameter(torch.randn(1))
    ]
    mock_model_instance.return_value = {'energy': torch.randn(1, requires_grad=True)}
    mock_prepare_data.return_value = labeled_atoms_list

    engine = TrainingEngine(training_config_no_delta)
    engine.execute(labeled_atoms_list)

    mock_prepare_data.assert_called_once_with(labeled_atoms_list)
    assert mock_mace.called
    mock_torch_save.assert_called_once()

@patch('mlip_autopipec.modules.d_training_engine.MACE')
@patch('mlip_autopipec.modules.d_training_engine.torch.save')
def test_training_engine_no_data(
    mock_torch_save, mock_mace, training_config_no_delta
):
    """Test that the engine handles the case where no labeled data is found."""
    engine = TrainingEngine(training_config_no_delta)
    engine.execute([])

    mock_mace.assert_not_called()
    mock_torch_save.assert_not_called()
