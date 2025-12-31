import os
import unittest
from unittest.mock import MagicMock, patch

import numpy as np
from ase import Atoms
from click.testing import CliRunner

from mlip_autopipec.data.models import DFTResult
from mlip_autopipec.main import cli


class TestCycle01Workflow(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()
        self.db_path = "test_e2e.db"
        if os.path.exists(self.db_path):
            os.remove(self.db_path)

    def tearDown(self):
        if os.path.exists(self.db_path):
            os.remove(self.db_path)
        if os.path.exists("model_cycle01.pt"):
            os.remove("model_cycle01.pt")

    @patch("mlip_autopipec.orchestrator.LabellingEngine")
    @patch("mlip_autopipec.orchestrator.TrainingEngine")
    @patch("torch.save")
    def test_e2e_happy_path(self, mock_torch_save, mock_training_engine, mock_labelling_engine):
        # Mock LabellingEngine to return a successful DFTResult
        mock_dft_result = DFTResult(
            energy=-10.0,
            forces=np.zeros((2, 3)),
            stress=np.zeros(6),
            status="success",
        )
        mock_labelling_engine.return_value.run.return_value = mock_dft_result

        # Mock TrainingEngine to return a dummy model
        mock_training_engine.return_value.train.return_value = MagicMock()

        result = self.runner.invoke(
            cli, ["run-cycle01", "--db-path", self.db_path, "--qe-command", "mock_pw.x"]
        )

        self.assertEqual(result.exit_code, 0)
        self.assertIn("Cycle 01 workflow finished.", result.output)

        # Verify that the database was created and populated
        self.assertTrue(os.path.exists(self.db_path))

        # Verify that the model was "saved"
        mock_torch_save.assert_called_once()


if __name__ == "__main__":
    unittest.main()
