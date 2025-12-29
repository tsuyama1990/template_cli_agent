# ruff: noqa: D101, D102, D103, D104, D105, D107, S101
import subprocess
from unittest.mock import MagicMock, patch

import pytest
from ase.atoms import Atoms

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine


@pytest.fixture
def mock_db() -> MagicMock:
    """Provides a mock AseDB instance."""
    db = MagicMock(spec=AseDB)
    db.write.return_value = 1  # Simulate returning a db_id
    return db


@pytest.fixture
def atoms_h2() -> Atoms:
    """Provides a simple H2 Atoms object."""
    return Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]], cell=[5, 5, 5], pbc=True)


@patch("subprocess.run")
def test_labelling_engine_execute(
    mock_subprocess_run: MagicMock, mock_db: MagicMock, atoms_h2: Atoms
):
    """Test the main execution flow of the LabellingEngine."""
    # Mock the return value of the subprocess call
    mock_process_result = MagicMock(spec=subprocess.CompletedProcess)
    # A minimal successful output string for the parser to work with
    mock_process_result.stdout = """
!    total energy              =     -2.0 Ry
Forces acting on atoms (cartesian axes, Ry/au):
     atom    1 type  1   force =     0.0  0.0  0.0
     atom    2 type  1   force =     0.0  0.0  0.0
total stress  (Ry/bohr**3)                (kbar)     P=       0.0
      0.0   0.0   0.0    0.0   0.0   0.0
    """
    mock_subprocess_run.return_value = mock_process_result

    engine = LabellingEngine(
        qe_command="mpirun -np 4 pw.x",
        db=mock_db,
        pseudo_potentials={"H": "H.pbe-rrkjus.UPF"},
    )

    db_id = engine.execute(atoms_h2)

    # Verify that subprocess.run was called correctly
    assert mock_subprocess_run.call_count == 1
    call_args = mock_subprocess_run.call_args[0][0]
    assert "mpirun" in call_args
    assert "pw.x" in call_args
    assert "-in" in call_args

    # Verify that the database write method was called
    assert mock_db.write.call_count == 1
    write_args = mock_db.write.call_args[0]
    written_atoms = write_args[0]
    written_result = write_args[1]

    assert isinstance(written_atoms, Atoms)
    assert written_result.was_successful
    assert written_result.total_energy_ev < 0

    # Verify the correct db_id is returned
    assert db_id == 1
