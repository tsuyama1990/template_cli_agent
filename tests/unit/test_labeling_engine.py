# tests/unit/test_labeling_engine.py

import subprocess
from unittest.mock import MagicMock

import pytest
from ase import Atoms

from mlip_autopipec.configs.models import DFTComputeConfig
from mlip_autopipec.modules.c_labeling_engine import LabelingEngine

# Reusing the sample QE output from the dft_utils test
SAMPLE_QE_OUTPUT = """
!    total energy              =     -15.85217439 Ry
Forces acting on atoms (cartesian axes, Ry/au):
     atom    1   force =     -0.00000014    -0.00000014    -0.00000014
total   stress  (Ry/bohr**3)                (kbar)     P=      -0.34
      -0.00000215    -0.00000000     0.00000000
      -0.00000000    -0.00000215     0.00000000
       0.00000000     0.00000000    -0.00000215
"""


@pytest.fixture
def mock_db_wrapper():
    """Fixture to create a mock AseDBWrapper."""
    db_wrapper = MagicMock()
    mock_row = MagicMock()
    mock_row.toatoms.return_value = Atoms("Si", cell=[1, 1, 1], pbc=True)
    mock_row.id = 1
    db_wrapper.get_unlabeled_rows.return_value = [mock_row]
    return db_wrapper


@pytest.fixture
def dft_config():
    """Fixture to provide a standard DFTComputeConfig."""
    return DFTComputeConfig(
        code="quantum_espresso",
        command="pw.x",
        pseudopotentials="SSSP",
        ecutwfc=60.0,
        ecutrho=240.0,
        kpoints_density=1.0,
        smearing="mv",
        degauss=0.01,
    )


def test_labeling_engine_run_success(
    mocker, mock_db_wrapper, dft_config, tmp_path
):
    """Tests the successful run of the LabelingEngine."""
    # Mock subprocess.run to return a successful result
    mock_subprocess_run = mocker.patch(
        "subprocess.run",
        return_value=subprocess.CompletedProcess(
            args=[], returncode=0, stdout=SAMPLE_QE_OUTPUT, stderr=""
        ),
    )

    engine = LabelingEngine(
        config=dft_config, db_wrapper=mock_db_wrapper, calculation_dir=tmp_path
    )
    engine.run()

    # Verify database interactions
    mock_db_wrapper.get_unlabeled_rows.assert_called_once()
    # Check that update_row was called with parsed data
    args, kwargs = mock_db_wrapper.update_row.call_args
    assert args[0] == 1  # row_id
    assert "energy" in args[1]
    assert "forces" in args[1]
    assert "stress" in args[1]
    assert args[2] == {"state": "labeled"}  # Corrected assertion

    # Verify that the QE command was called correctly
    expected_command = ["pw.x", "-in", str(tmp_path / "id_1/pw.in")]
    mock_subprocess_run.assert_called_once_with(
        expected_command, capture_output=True, text=True, check=True, cwd=mocker.ANY
    )


def test_labeling_engine_run_subprocess_error(
    mocker, mock_db_wrapper, dft_config, tmp_path
):
    """Tests that the engine correctly handles a CalledProcessError."""
    # Mock subprocess.run to raise an error
    mocker.patch(
        "subprocess.run",
        side_effect=subprocess.CalledProcessError(
            returncode=1,
            cmd=[],
            output="QE run failed",  # Corrected from 'stdout' to 'output'
            stderr="Error details",
        ),
    )

    engine = LabelingEngine(
        config=dft_config, db_wrapper=mock_db_wrapper, calculation_dir=tmp_path
    )
    engine.run()

    # Verify that the row is updated with an 'error' state
    mock_db_wrapper.update_row.assert_called_once_with(1, {}, {"state": "error"})


def test_labeling_engine_run_parsing_error(
    mocker, mock_db_wrapper, dft_config, tmp_path
):
    """Tests that the engine handles a failure in parsing the QE output."""
    # Mock a successful run but with invalid output that can't be parsed
    mocker.patch(
        "subprocess.run",
        return_value=subprocess.CompletedProcess(
            args=[], returncode=0, stdout="Invalid output", stderr=""
        ),
    )

    engine = LabelingEngine(
        config=dft_config, db_wrapper=mock_db_wrapper, calculation_dir=tmp_path
    )
    engine.run()

    # Verify that the row is updated with an 'error' state
    mock_db_wrapper.update_row.assert_called_once_with(1, {}, {"state": "error"})
