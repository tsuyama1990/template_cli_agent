from unittest.mock import MagicMock

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from mlip_autopipec.config import DFTComputeConfig
from mlip_autopipec.modules.labeling_engine import LabelingEngine


@pytest.fixture
def mock_process_runner():
    """Fixture for a mocked IProcessRunner."""
    runner = MagicMock()
    completed_process = MagicMock()
    completed_process.returncode = 0
    runner.run.return_value = completed_process
    return runner


@pytest.fixture
def mock_ase_io(mocker):
    """Fixture to mock ase.io.read and ase.io.write."""
    mock_write = mocker.patch("mlip_autopipec.modules.labeling_engine.write")

    mock_atoms = Atoms("H", positions=[(0, 0, 0)], cell=[10, 10, 10], pbc=True)
    mock_atoms.calc = SinglePointCalculator(
        atoms=mock_atoms,
        energy=-1.0,
        forces=np.array([[0.1, 0.2, 0.3]]),
        stress=np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6]),
    )
    mock_read = mocker.patch("mlip_autopipec.modules.labeling_engine.read", return_value=mock_atoms)

    return mock_write, mock_read


def test_labeling_engine_uses_dft_config(mock_process_runner, mock_ase_io):
    """
    Tests that the LabelingEngine correctly uses the DFTComputeConfig
    to generate the Quantum Espresso input file.
    """
    mock_write, _ = mock_ase_io

    dft_config = DFTComputeConfig(
        ecutwfc=50.0,
        ecutrho=400.0,
        kpoints_density=3.5,
        magnetism=None,
        control={"calculation": "vc-relax"},
        pseudopotentials={"H": "H.UPF"},
    )
    labeling_engine = LabelingEngine(
        dft_compute_config=dft_config,
        process_runner=mock_process_runner,
        qe_command="pw.x",
    )
    atoms_to_label = Atoms("H", positions=[(0, 0, 0)], cell=[10, 10, 10], pbc=True)
    result = labeling_engine.label_structure(atoms_to_label)

    # Assert that ase.io.write was called with the correct parameters
    mock_write.assert_called_once()
    _, kwargs = mock_write.call_args
    assert kwargs["format"] == "espresso-in"
    assert kwargs["ecutwfc"] == 50.0
    assert kwargs["ecutrho"] == 400.0
    assert kwargs["kpts"] == (4, 4, 4)  # Based on density and cell
    assert kwargs["control"]["calculation"] == "vc-relax"

    # Assert that the result is parsed correctly
    assert result.energy == -1.0
    assert np.allclose(result.forces, [[0.1, 0.2, 0.3]])


def test_labeling_engine_handles_qe_failure(mock_ase_io):
    """
    Tests that the LabelingEngine raises a RuntimeError if the
    Quantum Espresso calculation fails.
    """
    mock_process_runner = MagicMock()
    completed_process = MagicMock()
    completed_process.returncode = 1  # Simulate a failure
    mock_process_runner.run.return_value = completed_process

    dft_config = DFTComputeConfig(
        ecutwfc=50.0,
        ecutrho=400.0,
        kpoints_density=3.5,
        magnetism=None,
        control={"calculation": "vc-relax"},
        pseudopotentials={"H": "H.UPF"},
    )
    labeling_engine = LabelingEngine(
        dft_compute_config=dft_config,
        process_runner=mock_process_runner,
        qe_command="pw.x",
    )
    atoms_to_label = Atoms("H", positions=[(0, 0, 0)], cell=[10, 10, 10], pbc=True)

    with pytest.raises(RuntimeError) as excinfo:
        labeling_engine.label_structure(atoms_to_label)

    assert "Quantum Espresso calculation failed" in str(excinfo.value)
