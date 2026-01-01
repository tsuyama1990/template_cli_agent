from unittest.mock import MagicMock, patch, ANY

import pytest
from ase import Atoms

# Under test
from mlip_autopipec.modules.c_labeling_engine import LabelingEngine
from mlip_autopipec.configs.models import DFTComputeConfig

# Sample data that the mocked parser will return
PARSED_QE_DATA = {
    "energy": -150.0,
    "forces": [[0.1, 0.1, 0.1]],
    "stress": [[0.0] * 3] * 3,
}

SAMPLE_QE_OUTPUT_STR = "This is a fake QE output"


@pytest.fixture
def mock_db_wrapper():
    """Fixture for a mocked AseDBWrapper."""
    db = MagicMock()
    # Simulate get_unlabeled_rows returning one row to be processed
    mock_row = MagicMock()
    mock_row.id = 1
    mock_row.toatoms.return_value = Atoms("Si")
    db.get_unlabeled_rows.return_value = [mock_row]
    return db


@pytest.fixture
def dft_config():
    """Fixture for a simple DFTComputeConfig."""
    return DFTComputeConfig(
        code="quantum_espresso",
        command="pw.x",
        pseudopotentials="SSSP",
        ecutwfc=50.0,
        ecutrho=200.0,
        kpoints_density=3.0,
        smearing="mv",
        degauss=0.01,
    )

# Patch the dependencies of the LabelingEngine module
@patch("mlip_autopipec.modules.c_labeling_engine.subprocess.run")
@patch("mlip_autopipec.modules.c_labeling_engine.parse_qe_output")
@patch("mlip_autopipec.modules.c_labeling_engine.ase.io.write")
@patch("mlip_autopipec.modules.c_labeling_engine.Espresso")
def test_labeling_engine_run_correct_ase_usage(
    mock_espresso_calc, mock_ase_write, mock_parse_qe, mock_subprocess_run, mock_db_wrapper, dft_config
):
    """
    Test the main run loop of the LabelingEngine, verifying the correct
    ASE calculator pattern is used.
    """
    # --- Setup Mocks ---
    # Mock for subprocess.run
    mock_subprocess_run.return_value = MagicMock(
        returncode=0, stdout=SAMPLE_QE_OUTPUT_STR, stderr=""
    )
    # Mock for the parser
    mock_parse_qe.return_value = PARSED_QE_DATA
    # Mock for the Espresso calculator instance
    mock_calc_instance = MagicMock()
    mock_espresso_calc.return_value = mock_calc_instance

    # --- Instantiate and Run ---
    engine = LabelingEngine(config=dft_config, db_wrapper=mock_db_wrapper)

    # Temporarily mock open to simulate reading the output file
    with patch("builtins.open", MagicMock()) as mock_open:
        mock_open.return_value.__enter__.return_value.read.return_value = SAMPLE_QE_OUTPUT_STR
        engine.run()


    # --- Assertions ---
    # 1. Verify it fetched unlabeled data
    mock_db_wrapper.get_unlabeled_rows.assert_called_once()

    # 2. Verify an Espresso calculator was instantiated with the correct params
    mock_espresso_calc.assert_called_once()
    calc_args = mock_espresso_calc.call_args[1]
    assert calc_args['input_data']['system']['ecutwfc'] == dft_config.ecutwfc

    # 3. Verify ase.io.write was called. We don't care about the atoms object itself (ANY)
    # but we care that it was called with the correct format.
    mock_ase_write.assert_called_once()
    write_args = mock_ase_write.call_args[0]
    # Check that the second argument (the atoms object) has our mocked calculator attached
    assert write_args[1].calc is mock_calc_instance

    # 4. Verify the external command was run
    mock_subprocess_run.assert_called_once()
    command_list = mock_subprocess_run.call_args[0][0]
    assert dft_config.command in command_list
    assert "-in" in command_list

    # 5. Verify the output file was opened and read
    mock_open.assert_called_with(ANY, "r")

    # 6. Verify the output was parsed
    mock_parse_qe.assert_called_once_with(SAMPLE_QE_OUTPUT_STR)

    # 7. Verify the database was updated
    mock_db_wrapper.update_row.assert_called_once_with(
        row_id=1,
        data=PARSED_QE_DATA,
        key_value_pairs={"state": "labeled"}
    )
