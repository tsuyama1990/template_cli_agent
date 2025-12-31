import unittest
from unittest.mock import MagicMock, patch

import numpy as np
import torch
from ase import Atoms

from mlip_autopipec.data.models import DFTResult
from mlip_autopipec.modules.training_engine import TrainingEngine

# Mock the MACE model to avoid actual training
mock_mace = MagicMock()
mock_mace_instance = MagicMock()
mock_mace_instance.return_value = {
    "energy": torch.randn(1, requires_grad=True),
    "forces": torch.randn(2, 3, requires_grad=True),
}
mock_mace.return_value = mock_mace_instance


class TestTrainingEngine(unittest.TestCase):
    def setUp(self):
        self.engine = TrainingEngine(model_type="mace")
        self.atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])
        self.dft_result = DFTResult(
            energy=-10.0,
            forces=np.random.rand(2, 3),
            stress=np.zeros(6),
            status="success",
        )
        self.training_data = [(self.atoms, self.dft_result)]

    def test_compute_delta(self):
        baseline = self.engine._get_baseline_potential(["H"])
        delta_energy, delta_forces = self.engine._compute_delta(
            self.atoms, self.dft_result, baseline
        )
        # Check that the delta is not equal to the original DFT result
        self.assertNotAlmostEqual(delta_energy, self.dft_result.energy)
        self.assertFalse(np.allclose(delta_forces, self.dft_result.forces))

    @patch("mace.modules.models.MACE", new=mock_mace)
    @patch("torch.optim.Adam")
    def test_train_success(self, mock_adam):
        model = self.engine.train(self.training_data)
        self.assertIsNotNone(model)
        # Check that the MACE model was instantiated
        mock_mace.assert_called()
        # Check that the optimizer was instantiated
        mock_adam.assert_called()
        # Check that the model's forward and backward passes were called
        self.assertTrue(mock_mace_instance.call_count > 0)

    def test_train_empty_data(self):
        with self.assertRaises(ValueError):
            self.engine.train([])


if __name__ == "__main__":
    unittest.main()
