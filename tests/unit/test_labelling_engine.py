import unittest
from unittest.mock import MagicMock, patch

import numpy as np
from ase import Atoms

from mlip_autopipec.modules.labelling_engine import LabellingEngine

# Sample QE output for mocking
SUCCESS_QE_OUTPUT = """
!    total energy              =     -29.95595568 Ry
Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =     0.00000000    0.00000000   -0.00000000
     atom    2 type  1   force =    -0.00000000   -0.00000000    0.00000000

The total force is     0.00000000    Total SCF correction is     0.00000000
"""

FAILURE_QE_OUTPUT = """
     Error: SCF not converged in 100 steps
"""


class TestLabellingEngine(unittest.TestCase):
    def setUp(self):
        self.engine = LabellingEngine(
            qe_command="pw.x -in stdin", pseudo_dir="/path/to/pseudos"
        )
        self.atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])

    @patch("mlip_autopipec.modules.labelling_engine.subprocess.run")
    def test_run_success(self, mock_subprocess_run):
        # Mock a successful QE run
        mock_process = MagicMock()
        mock_process.stdout = SUCCESS_QE_OUTPUT
        mock_process.stderr = ""
        mock_process.returncode = 0
        mock_subprocess_run.return_value = mock_process

        result = self.engine.run(self.atoms)

        self.assertEqual(result.status, "success")
        self.assertAlmostEqual(result.energy, -29.95595568 * 13.605693122994, places=5)
        self.assertIsInstance(result.forces, np.ndarray)
        self.assertEqual(result.forces.shape, (2, 3))
        mock_subprocess_run.assert_called_once()

    @patch("mlip_autopipec.modules.labelling_engine.subprocess.run")
    def test_run_failure(self, mock_subprocess_run):
        # Mock a failed QE run
        mock_process = MagicMock()
        mock_process.stdout = FAILURE_QE_OUTPUT
        mock_process.stderr = "SCF not converged"
        mock_process.returncode = 1
        mock_subprocess_run.return_value = mock_process

        result = self.engine.run(self.atoms)

        self.assertEqual(result.status, "failed")
        self.assertIn("stderr", result.metadata)
        self.assertEqual(result.metadata["stderr"], "SCF not converged")

    def test_generate_input_file(self):
        input_file = self.engine._generate_input_file(self.atoms)
        self.assertIn("calculation      = 'scf'", input_file)
        self.assertIn("pseudo_dir       = '/path/to/pseudos'", input_file)
        self.assertIn("nat              = 2", input_file)
        self.assertIn("H 0.0000000000 0.0000000000 0.0000000000", input_file)


if __name__ == "__main__":
    unittest.main()
