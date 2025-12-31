from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from ase import Atoms

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.data.database import AseDB
from mlip_autopipec.orchestrator import Orchestrator


@pytest.fixture
def temp_db_path(tmp_path):
    """Fixture to create a temporary database file path for the test."""
    return str(tmp_path / "test_e2e.db")


@pytest.fixture
def test_config(temp_db_path):
    """Provides a full configuration for the e2e test."""
    return FullConfig(
        ase_db_path=temp_db_path,
        dft_compute={
            "command": "pw.x -in fake.pwi",
            "pseudopotentials": {"H": "H.upf"},
            "ecutwfc": 30.0,
        },
        training={"r_cut": 4.0, "delta_learning": False},
    )


def load_qe_output(filename: str) -> str:
    """Helper to load a QE output file from the test data directory."""
    path = Path(__file__).parent / ".." / "unit" / "test_data" / "qe_outputs" / filename
    return path.read_text()


@patch("torch.save")
@patch("subprocess.run")
def test_cycle01_e2e_workflow(mock_subprocess_run, mock_torch_save, test_config, temp_db_path):
    """
    Tests the full end-to-end workflow for Cycle 01, mocking the external
    Quantum Espresso execution.
    """
    # --- 1. Setup ---
    # Mock the external process call
    mock_process = MagicMock()
    mock_process.returncode = 0
    # Provide a realistic output to be parsed
    mock_process.stdout = load_qe_output("successful_run.out").replace("Si", "H ")
    mock_process.stderr = ""
    mock_subprocess_run.return_value = mock_process

    # Create a database and add unlabelled structures
    db = AseDB(db_path=temp_db_path)
    structures_to_label = [
        Atoms("H2", positions=[[0, 0, 0], [0.8, 0, 0]]),
        Atoms("H2", positions=[[0, 0, 0], [1.5, 0, 0]]),
    ]
    db.add_structures(structures_to_label, status="needs_labelling")

    # --- 2. Execution ---
    # Instantiate and run the orchestrator
    orchestrator = Orchestrator(config=test_config)
    orchestrator.run_cycle_01()

    # --- 3. Verification ---
    # Verify that the external command was called for each structure
    assert mock_subprocess_run.call_count == 2

    # Verify database state changes
    labelled = db.get_structures_by_status("labelled")
    needs_labelling = db.get_structures_by_status("needs_labelling")
    assert len(labelled) == 2
    assert len(needs_labelling) == 0

    # Verify that the labelled atoms have calculator results attached
    atoms = labelled[0]
    assert atoms.calc is not None
    # Check for a realistic energy value from the parsed output
    assert atoms.get_potential_energy() == pytest.approx(-15.85311181 * 13.605693122994)

    # Verify that the model training was run and saved
    mock_torch_save.assert_called_once()
    # Check that the saved object is a state dictionary (a dict)
    assert isinstance(mock_torch_save.call_args[0][0], dict)
