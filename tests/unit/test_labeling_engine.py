import subprocess
from unittest.mock import MagicMock, patch

import pytest
from ase import Atoms

from mlip_autopipec.config import DFTInputConfig, DFTResult
from mlip_autopipec.interfaces import IProcessRunner
from mlip_autopipec.modules.labeling_engine import LabelingEngine


@pytest.fixture
def mock_process_runner():
    """Fixture for a mocked IProcessRunner."""
    return MagicMock(spec=IProcessRunner)


@pytest.fixture
def dft_config():
    """Fixture for a sample DFTInputConfig."""
    return DFTInputConfig(
        pseudopotentials={"H": "H.UPF", "O": "O.UPF"},
        kpoints=(1, 1, 1),
        ecutwfc=60,
        control={"calculation": "scf"},
    )


def test_label_structure_success(mock_process_runner, dft_config):
    """Tests the successful execution of the label_structure method."""
    # Arrange
    atoms_to_label = Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])

    def mock_run(command, stdout_path):
        with open(stdout_path, "w") as f:
            f.write(
                """
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
!    total energy              =     -17.83123456 Ry
     Forces acting on atoms (cartesian axes, Ry/au):
          atom    1   type  O   force =     0.001   0.002   0.003
          atom    2   type  H   force =    -0.001  -0.001  -0.001
          atom    3   type  H   force =     0.000   0.000  -0.002
"""
            )

    mock_process_runner.run.side_effect = mock_run

    engine = LabelingEngine(
        dft_input_configuration=dft_config,
        process_runner=mock_process_runner,
        qe_command="pw.x",
    )

    # Act
    result = engine.label_structure(atoms_to_label)

    # Assert
    assert isinstance(result, DFTResult)
    assert result.energy == pytest.approx(-242.606, abs=1e-3)
    mock_process_runner.run.assert_called_once()


def test_label_structure_qe_failure(mock_process_runner, dft_config):
    """Tests the error handling when Quantum Espresso execution fails."""
    # Arrange
    atoms_to_label = Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    mock_process_runner.run.side_effect = subprocess.CalledProcessError(returncode=1, cmd="pw.x")
    engine = LabelingEngine(
        dft_input_configuration=dft_config,
        process_runner=mock_process_runner,
        qe_command="pw.x",
    )

    # Act & Assert
    with pytest.raises(subprocess.CalledProcessError):
        engine.label_structure(atoms_to_label)


@patch("mlip_autopipec.modules.labeling_engine.read")
def test_label_structure_malformed_output(mock_ase_read, mock_process_runner, dft_config):
    """Tests the error handling when the QE output is malformed."""
    # Arrange
    atoms_to_label = Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    mock_ase_read.side_effect = Exception("Malformed output")

    def mock_run(command, stdout_path):
        # The content of the output file is irrelevant for this test,
        # as the read function is mocked.
        with open(stdout_path, "w") as f:
            f.write("dummy content")

    mock_process_runner.run.side_effect = mock_run

    engine = LabelingEngine(
        dft_input_configuration=dft_config,
        process_runner=mock_process_runner,
        qe_command="pw.x",
    )

    # Act & Assert
    with pytest.raises(Exception, match="Malformed output"):
        engine.label_structure(atoms_to_label)
