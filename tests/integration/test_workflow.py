import subprocess
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest
from ase import Atoms
from click.testing import CliRunner

from mlip_autopipec.cli import app
from mlip_autopipec.database import AseDBWrapper
from mlip_autopipec.interfaces import IProcessRunner

# A more realistic dummy Quantum Espresso output file content
DUMMY_QE_OUTPUT = """
     Program PWSCF v.6.5 starts on 12Jan2021 at 12:00:00

     lattice parameter (alat)  =  18.897260  a.u.
     celldm(1)   18.897260

     number of atoms/cell      =           3
     number of atomic types    =           2

     crystal axes: (angstrom)
               a(1) = (   10.000   0.000   0.000 )
               a(2) = (    0.000  10.000   0.000 )
               a(3) = (    0.000   0.000  10.000 )

     site n.     atom                  positions (alat units)
         1           O      tau(   1) = (   0.00000   0.00000   0.01170  )
         2           H      tau(   2) = (   0.00000   0.07570  -0.04690  )
         3           H      tau(   3) = (   0.00000  -0.07570  -0.04690  )

     Forces acting on atoms (cartesian axes, Ry/au):

          atom    1   type  O   force =     0.001   0.002   0.003
          atom    2   type  H   force =    -0.001  -0.001  -0.001
          atom    3   type  H   force =     0.000   0.000  -0.002

!    total energy              =     -17.83123456 Ry

     total stress  (Ry/bohr**3)                (kbar)     P=       -0.01
      -0.00000011   0.00000000   0.00000000     -0.01      0.00      0.00
       0.00000000  -0.00000015   0.00000000      0.00     -0.02      0.00
       0.00000000   0.00000000  -0.00000020      0.00      0.00     -0.03
"""


@pytest.fixture
def mock_process_runner(mocker):
    """Mocks IProcessRunner to simulate a QE execution."""

    def mock_run(command, stdout_path):
        with open(stdout_path, "w") as f:
            f.write(DUMMY_QE_OUTPUT)
        return subprocess.CompletedProcess(args=command, returncode=0)

    mock = MagicMock(spec=IProcessRunner)
    mock.run.side_effect = mock_run
    mocker.patch(
        "mlip_autopipec.factories.SubprocessRunner", return_value=mock
    )
    return mock


def test_label_and_train_workflow(mock_process_runner, tmp_path):
    """An integration test for the full label-and-train workflow."""
    runner = CliRunner()
    with runner.isolated_filesystem() as temp_dir:
        db_path = Path(temp_dir) / "test.db"
        db_wrapper = AseDBWrapper(db_path=str(db_path))

        atoms = Atoms("H2O", positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0)])
        structure_id = db_wrapper.add_atoms(atoms)
        assert structure_id == 1

        result = runner.invoke(
            app,
            ["label", "--id", str(structure_id), "--db-path", str(db_path)],
            catch_exceptions=False,
        )
        assert result.exit_code == 0
        assert (
            f"Labeling complete for structure ID: {structure_id}"
            in result.output
        )

        mock_process_runner.run.assert_called_once()

        labeled_data = db_wrapper.get_all_labeled_atoms()
        assert len(labeled_data) == 1
        _, dft_result = labeled_data[0]
        assert dft_result.energy == pytest.approx(-242.62, abs=1e-1)
        np.testing.assert_array_equal(dft_result.stress, np.zeros((3, 3)))

        result = runner.invoke(
            app, ["train", "--db-path", str(db_path)], catch_exceptions=False
        )
        assert result.exit_code == 0
        assert "Training not implemented" in result.output
