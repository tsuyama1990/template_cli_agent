import subprocess
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from click.testing import CliRunner

from mlip_autopipec.cli import app
from mlip_autopipec.database import AseDBWrapper

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
def mock_subprocess_run(mocker):
    """Mocks subprocess.run to simulate a QE execution."""

    def mock_run(*args, **kwargs):
        command_list = args[0]
        # The output is now redirected to stdout, which is a file object
        output_file_handle = kwargs.get("stdout")
        if output_file_handle:
            output_file_handle.write(DUMMY_QE_OUTPUT)
            output_file_handle.flush()
        return subprocess.CompletedProcess(args=command_list, returncode=0)

    return mocker.patch("subprocess.run", side_effect=mock_run)


def test_label_and_train_workflow(mock_subprocess_run, tmp_path):
    """An integration test for the full label-and-train workflow."""
    # Use CliRunner's built-in temporary file system
    runner = CliRunner()
    with runner.isolated_filesystem() as temp_dir:
        db_path = Path(temp_dir) / "asedb.db"  # Match the hardcoded path in cli.py
        db_wrapper = AseDBWrapper(db_path=str(db_path))

        # 1. Add an initial unlabeled structure
        atoms = Atoms("H2O", positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0)])
        id = db_wrapper.add_atoms(atoms)
        assert id == 1

        # 2. Run the labeling CLI command
        result = runner.invoke(
            app, ["label", "--id", str(id)], catch_exceptions=False
        )
        assert result.exit_code == 0
        assert f"Labeling complete for structure ID: {id}" in result.output

        # 3. Verify the subprocess was called
        mock_subprocess_run.assert_called_once()

        # 4. Verify the database was updated
        labeled_data = db_wrapper.get_all_labeled_atoms()
        assert len(labeled_data) == 1
        retrieved_atoms, dft_result = labeled_data[0]
        assert dft_result.energy == pytest.approx(-242.62, abs=1e-1)
        # Verify that stress defaults to zero when parsing fails
        np.testing.assert_array_equal(dft_result.stress, np.zeros((3, 3)))

        # 5. Run the training command (which is expected to fail gracefully)
        result = runner.invoke(app, ["train"], catch_exceptions=False)
        assert result.exit_code == 0
        assert "Training not implemented" in result.output
