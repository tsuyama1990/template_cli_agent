import pytest
from unittest.mock import patch, MagicMock
from ase import Atoms
import numpy as np
import subprocess
import unittest
import os

from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.data.database import AseDB

# Get the absolute path of the test directory
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(TEST_DIR, '..', '..', 'data')


@pytest.fixture
def mock_db(tmp_path):
    """Fixture for a mock AseDB."""
    db_path = tmp_path / "test.db"
    return AseDB(str(db_path))

@pytest.fixture
def labelling_engine(mock_db):
    """Fixture for a LabellingEngine with a mock DB."""
    return LabellingEngine(
        qe_command="mock_pw.x",
        db=mock_db,
        default_params={},
        pseudos={'H': 'H.upf'}
    )

@patch('subprocess.run')
def test_labelling_engine_success(mock_subprocess_run, labelling_engine, mock_db):
    """Tests a successful run of the LabellingEngine."""
    with open(os.path.join(DATA_DIR, 'qe_success.out'), 'r') as f:
        successful_output = f.read()

    # Configure the mock to simulate a successful run and to write the output to the file
    def side_effect(*args, **kwargs):
        with open(kwargs['stdout'].name, 'w') as f:
            f.write(successful_output)
        return MagicMock(returncode=0, stderr="")

    mock_subprocess_run.side_effect = side_effect

    atoms = Atoms('H', positions=[[0, 0, 0]])
    db_id = labelling_engine.execute(atoms)

    assert db_id is not None
    mock_subprocess_run.assert_called_once()

    retrieved_atoms, kvp = mock_db.get(db_id)
    assert kvp['was_successful']
    assert np.isclose(retrieved_atoms.get_potential_energy(), -11.23032512 * 13.605693122994)


@patch('subprocess.run')
def test_labelling_engine_failure(mock_subprocess_run, labelling_engine, mock_db):
    """Tests a failed run of the LabellingEngine."""
    with open(os.path.join(DATA_DIR, 'qe_fail.out'), 'r') as f:
        failed_output = f.read()

    def side_effect(*args, **kwargs):
        with open(kwargs['stdout'].name, 'w') as f:
            f.write(failed_output)
        raise subprocess.CalledProcessError(returncode=1, cmd="mock_pw.x", stderr="SCF NOT CONVERGED")

    mock_subprocess_run.side_effect = side_effect

    atoms = Atoms('H', positions=[[0, 0, 0]])
    db_id = labelling_engine.execute(atoms)

    assert db_id is not None

    _, kvp = mock_db.get(db_id)
    assert not kvp['was_successful']
    assert "SCF did not converge" in kvp['error_message']
