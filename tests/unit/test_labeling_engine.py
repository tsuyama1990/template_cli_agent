import subprocess
from unittest.mock import MagicMock, patch

import pytest
from ase import Atoms

from mlip_autopipec.config import DFTInputConfig
from mlip_autopipec.database import AseDBWrapper
from mlip_autopipec.modules.labeling_engine import LabelingEngine


@pytest.fixture
def mock_db_wrapper():
    """Fixture for a mocked AseDBWrapper."""
    mock = MagicMock(spec=AseDBWrapper)
    mock.get_atoms_by_id.return_value = Atoms(
        "H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    )
    return mock


@pytest.fixture
def dft_config():
    """Fixture for a sample DFTInputConfig."""
    return DFTInputConfig(
        pseudopotentials={"H": "H.UPF", "O": "O.UPF"},
        kpoints=(1, 1, 1),
        ecutwfc=60,
        control={"calculation": "scf"},
    )


@patch("subprocess.run")
def test_label_structure_success(mock_subprocess_run, mock_db_wrapper, dft_config):
    """Tests the successful execution of the label_structure method."""
    # Arrange
    def mock_run(*args, **kwargs):
        output_file_handle = kwargs.get("stdout")
        if output_file_handle:
            output_file_handle.write(
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
            output_file_handle.flush()
        return subprocess.CompletedProcess(args=args[0], returncode=0)

    mock_subprocess_run.side_effect = mock_run

    engine = LabelingEngine(
        dft_config=dft_config,
        db_wrapper=mock_db_wrapper,
        qe_command="pw.x",
    )

    # Act
    engine.label_structure(1)

    # Assert
    mock_db_wrapper.get_atoms_by_id.assert_called_once_with(1)
    mock_subprocess_run.assert_called_once()
    mock_db_wrapper.update_labels.assert_called_once()


@patch("subprocess.run")
def test_label_structure_qe_failure(
    mock_subprocess_run, mock_db_wrapper, dft_config
):
    """Tests the error handling when Quantum Espresso execution fails."""
    # Arrange
    mock_subprocess_run.side_effect = subprocess.CalledProcessError(
        returncode=1, cmd="pw.x"
    )
    engine = LabelingEngine(
        dft_config=dft_config,
        db_wrapper=mock_db_wrapper,
        qe_command="pw.x",
    )

    # Act
    engine.label_structure(1)

    # Assert
    mock_db_wrapper.get_atoms_by_id.assert_called_once_with(1)
    mock_subprocess_run.assert_called_once()
    mock_db_wrapper.update_state.assert_called_once_with(
        1, "labeling_failed"
    )
    mock_db_wrapper.update_labels.assert_not_called()
