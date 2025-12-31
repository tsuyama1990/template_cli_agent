from unittest.mock import MagicMock

import pytest
from ase import Atoms

from mlip_autopipec.data.models import DFTCompute
from mlip_autopipec.modules.c_labeling_engine import LabelingEngine


@pytest.fixture
def mock_db_wrapper():
    return MagicMock()


@pytest.fixture
def mock_dft_config():
    return DFTCompute(
        code="quantum_espresso",
        command="pw.x -np 4",
        pseudopotentials="SSSP",
        ecutwfc=60.0,
        ecutrho=240.0,
        kpoints_density=2.5,
    )


def test_labeling_engine_execute(mocker, mock_db_wrapper, mock_dft_config):
    # Arrange
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]], cell=[10, 10, 10])
    mock_row = MagicMock()
    mock_row.id = 1
    mock_row.toatoms.return_value = atoms

    mock_db_wrapper.get_rows_to_label.return_value = [mock_row]

    mock_subprocess_run = mocker.patch("subprocess.run")

    # Mock file I/O but allow directories to be created
    mocker.patch("pathlib.Path.write_text")
    mocker.patch(
        "pathlib.Path.read_text",
        return_value="!    total energy = -1.23 Ry\nForces acting on atoms (cartesian axes, Ry/au):\n\n     atom    1 type  1   force =     0.0   0.0    0.0\n     atom    2 type  1   force =     0.0   0.0    0.0",
    )

    # Act
    engine = LabelingEngine(mock_db_wrapper, mock_dft_config)
    engine.execute()

    # Assert
    mock_db_wrapper.get_rows_to_label.assert_called_once()

    # Check subprocess call
    mock_subprocess_run.assert_called_once()

    # Check database update
    mock_db_wrapper.update_row_with_dft_results.assert_called_once()
    updated_atoms = mock_db_wrapper.update_row_with_dft_results.call_args[0][1]
    assert updated_atoms.get_potential_energy() is not None
